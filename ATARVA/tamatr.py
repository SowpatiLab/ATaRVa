import sys, os
import pysam
import threading
import polars as pl
from functools import reduce

def processor(process_df, outfile, tidx, each_thread, total_samples):
    # print('starting process = ', each_thread)
    #print(process_df.columns)
    out = open(f'{outfile}_reader{tidx}_processor{each_thread}.vcf', 'w')
    for row_dict in process_df.iter_rows(named=True):
        genotyped_samples = 0
        sample_wise_full_gt = []
        ALT = []
        alt_seq_lens = []
        alt_seq_count = {}
        for file_id in range(total_samples):
            current_sample = row_dict[f's{file_id}']
            if current_sample:
                splited_sample = current_sample.split(':')
                individual_sample_gt = splited_sample[1]
            else:
                sample_wise_full_gt.append('.:.:.:.:.:.:.:.:.')
                continue
            if individual_sample_gt=='.':
                sample_wise_full_gt.append('.:.:.:.:.:.:.:.:.')
            else:
                genotyped_samples += 1
                GT = []
                alt_seqs = splited_sample[0].split(',') if splited_sample[0]!='.' else ""
                seq_lens = [0 if i=='<DEL>' else len(i) for i in alt_seqs]
                for idx,lens in enumerate(seq_lens):
                    if lens in alt_seq_lens:
                        alt_seq_count[lens] += 1 # count of that alt allele
                        GT.append(str(alt_seq_lens.index(lens) + 1))
                    else:
                        ALT.append(alt_seqs[idx])
                        alt_seq_lens.append(lens)
                        alt_seq_count[lens] = 1 # initialize count of that alt allele
                        GT.append(str(len(alt_seq_lens)))
                alt_count = len(GT)
                if len(individual_sample_gt) > 1: # autosomes
                    phaser = individual_sample_gt[1] # either '/' or '|'
                    sep_gt = individual_sample_gt.split(phaser) # separated genotype
                            
                    if alt_count==2: # if there are two alt alleles
                        new_GT = phaser.join(GT)
                    elif alt_count == 1: # if there is only one alt alleles
                        the_single_gt = GT[0]
                        if len(set(sep_gt)) == 2: # if it is heterozyous
                            new_GT = phaser.join(['0', the_single_gt])
                        else: # if it is homozygous
                            new_GT = phaser.join([the_single_gt, the_single_gt])
                    else:
                        new_GT = '0'+phaser+'0' 
                else: # Sex chromosomes
                    if alt_count==1:
                        new_GT = str(GT[0])
                    else:
                        new_GT = '0'                        

                splited_sample[1] = new_GT
                sample_wise_full_gt.append(':'.join(splited_sample[1:]))
        if genotyped_samples:
            pass
        else:
            continue
        if alt_seq_lens:
            AC = []
            for i in alt_seq_lens:
                AC.append(str(alt_seq_count[i]))
            AC = ','.join(AC)
        else:
            AC = '0'
        AN = str(genotyped_samples * 2)
        info = 'AC='+AC+';AN='+AN+';' + row_dict['i']
        ref_seq = row_dict['r']
        start = row_dict['s']
        chrom = row_dict['c']
        filter = '.'
        id = '.'
        q = '.'
        alt = ','.join(ALT) if ALT else '.'
        format = 'GT:AL:AR:SD:DP:SN:SQ:MM:MR'
        #sample = '\t'.join(sample_wise_full_gt)

        repeat_info = [chrom, start, id, ref_seq, alt, q, filter, info, format, *sample_wise_full_gt]
        del sample_wise_full_gt
        tot_tabs = len(repeat_info)
        chunk_size = 100
        for i in range(0, tot_tabs, chunk_size):
            chunk = repeat_info[i:i + chunk_size]
            out.write("\t".join(map(str, chunk)))
            if i<tot_tabs-1:
                out.write("\t")
        out.write("\n")
        #out.write("\t".join(map(str, repeat_info)) + "\n")
        del repeat_info
    out.close()
    # print('DONE Processing....')

def reader(outfile, bedfile, ref, vcfs, contigs, tidx, process_thread):

    total_samples = len(vcfs)
    #print("total_samples = ", total_samples)
    tbx = pysam.TabixFile(bedfile)
    ref_file = pysam.FastaFile(ref)
    vcf_instance = []
    for each_vcf in vcfs:
        vcf_instance.append(pysam.TabixFile(each_vcf))
    #print(len(vcf_instance))
    if tidx!=-1: # multi thread
        if tidx==0: # first process
            vcf_names = [file_path.split("/")[-1].split('.')[0] for file_path in vcfs]
            out = open(f'{outfile}.vcf', 'w')
            vcf_writer(out, vcf_names, vcfs[0])
        else:
            out = open(f'{outfile}_thread_{tidx}.vcf', 'w')
    else: # single thread
        vcf_names = [file_path.split("/")[-1].split('.')[0] for file_path in vcfs]
        out = open(f'{outfile}.vcf', 'w')
        vcf_writer(out, vcf_names, vcfs[0])


    thread_pool = list()
    print('Reader thread = ', tidx)
    print(f'Inside reader{tidx} = length of contig = {len(contigs)}')
    for contig in contigs:
        
        Chrom, Start, End = contig
        # print("\nReading new block............")
        frames = []
        base_frame = pl.DataFrame().lazy()
        parquet_batch = 0
        file_count = 0
        for file_id,file in enumerate(vcf_instance):
            #print(file_id)
            #if Chrom not in file.contigs: break
            file_data_dict = {}
        
            # dictionary for sample file
            file_data_dict['s'] = []
            file_data_dict['e'] = []
        
            # variable for each column
            file_start = file_data_dict['s']
            file_end = file_data_dict['e']
        
            if file_id == 0:

                schema = {"s": pl.Int32,
                         "e": pl.Int32,
                         "c": pl.Categorical,
                         "r": pl.Categorical,
                         "i": pl.Categorical,
                         "s0": pl.Categorical}
                
                # dictionary for sample file
                file_data_dict['c'] = [] # chrom
                file_data_dict['r'] = [] # ref
                file_data_dict['i'] = [] # info
                file_data_dict['s0'] = [] # sample
                
                # variable for each column
                file_ref = file_data_dict['r']
                file_info = file_data_dict['i']
                file_sample = file_data_dict['s0']
             

                for line in tbx.fetch(Chrom, Start[0], End[1]):
                    line = line.strip().split('\t')
                    chrom = line[0]
                    start = int(line[1])
                    end = int(line[2])
                    
                    if (start>=Start[0]) and (end<=End[1]):
                        if start==Start[0]:
                            if end==Start[1]: pass
                            else: continue
                        pass
                    elif start<Start[0]:
                        continue
                    elif start>=End[0]: break
                    
                    motif_value = line[3]
                    ref_value = end-start # +1)//period_value
                    del line
                    
                    ref_string = True
                    has_region = False

                    if Chrom in file.contigs:
                        for entry in file.fetch(chrom, start+1, end):
                            entry = entry.strip().split('\t')
                            st = int(entry[1])
                            if (st-1)!=start: # -1 to match with 0-based coord
                                continue

                            info = entry[7].split(';', 5)[:5]
                            en = int(info[4].split('=')[1])
                            if ((st-1)==start) & (en==end): # -1 to match with 0-based coord
                                has_region = True
                                file_start.append(st) 
                                file_end.append(en)
                                file_ref.append(entry[3])
                                file_info.append(f"MOTIF={motif_value};START={start};END={end}")
                                sample = entry[9]
                                if sample[0]=='.':
                                    file_sample.append(None)
                                else:
                                    file_sample.append(entry[4]+':'+':'.join(sample.split(':', 9)[:9]))
                                del entry
                                del sample
                            break
                            
                    if not has_region:
                        file_start.append(start+1)
                        file_end.append(end)
                        file_ref.append(ref_file.fetch(chrom, start, end))
                        file_info.append(f"MOTIF={motif_value};START={start};END={end}")
                        file_sample.append(None)

                file_count += 1
                
                file_data_dict['c'].extend([Chrom]*len(file_start))
                df = pl.DataFrame(file_data_dict, schema=schema).lazy()
                df = df.unique(subset=['s', 'e'], keep='first', maintain_order=True)
                frames.append(df)
                base_frame = df.collect().select(['s', 'e']).lazy()
                del df
                #print('base_frame shape = ', base_frame.collect().shape)
            else:
                file_count += 1
                
                schema = {"s": pl.Int32,
                         "e": pl.Int32,
                         f's{file_id}': pl.Categorical}

                file_data_dict[f's{file_id}'] = [] # sample
                file_sample = file_data_dict[f's{file_id}']

                if file_count >= 200:
                    joiner(frames, parquet_batch, tidx, outfile)
                    parquet_batch += 1
                    del frames
                    frames = []
                    frames.append(base_frame)
                    file_count = 0

                if Chrom not in file.contigs:
                    df = pl.DataFrame(file_data_dict, schema=schema).lazy()
                    frames.append(df)
                    print(f'Continuing due to no chr {Chrom} in {file_id}')
                    continue

                for entry in file.fetch(Chrom, Start[0], End[1]):
                    entry = entry.strip().split('\t')
                    
                    sample = entry[9]
                    if sample[0] == '.':
                        continue
                    else:
                        st = int(entry[1])
                        info = entry[7].split(';', 5)[:5]
                        en = int(info[4].split('=')[1])
                        
        
                        file_start.append(st)
                        file_end.append(en)
                        file_sample.append(entry[4]+':'+':'.join(sample.split(':', 9)[:9]))
                        
                df = pl.DataFrame(file_data_dict, schema=schema).lazy()
                df = df.unique(subset=['s', 'e'], keep='first', maintain_order=True)
                frames.append(df)
                del df

        if frames:
            joiner(frames, parquet_batch, tidx, outfile)
            parquet_batch += 1
            del frames

        # print('Done reading & joining!!!!!!!!!')
        if thread_pool:
            # joining previous threads - waiting for previous threads to be over
            for thread_x in thread_pool:
                # print('waiting for ', thread_x)
                thread_x.join()
            thread_pool.clear()

            # print('Concatenating processor files..............')
            for each_thread in range(process_thread):
                thread_out = f'{outfile}_reader{tidx}_processor{each_thread}.vcf'
                # print('opening ', thread_out)
                with open(thread_out, 'r') as fh:
                    for line in fh:
                        repeat_info = line.strip().split('\t')
                        tot_tabs = len(repeat_info)
                        chunk_size = 100
                        for i in range(0, tot_tabs, chunk_size):
                            chunk = repeat_info[i:i + chunk_size]
                            out.write("\t".join(map(str, chunk)))
                            if i<tot_tabs-1:
                                out.write("\t")
                        out.write("\n")
                        del repeat_info
                        #out.write("\t".join(map(str, repeat_info)) + "\n")
                # print('Removing ', thread_out)
                #del repeat_info
                os.remove(thread_out)

        batch_files = [f"{outfile}_reader{tidx}_batch{batch_val}.parquet" for batch_val in range(parquet_batch)]
        parquet_frames = [pl.read_parquet(f).lazy() for f in batch_files]
        for p_files in batch_files:
            os.remove(p_files)
        
        if parquet_frames:
            merged = reduce(lambda l, r: l.join(r, on=['s','e'], how='left'), parquet_frames)
            whole_df = merged.collect(engine="streaming")
            print("Shape = ", whole_df.shape)
            
            if process_thread > 0:
                loci_count = whole_df.shape[0]
                split_count = loci_count // process_thread
                if split_count == 0:
                    split_count = 1
                initial = 0
                track = split_count
        
                # initializing threads
                for each_thread in range(process_thread):
                    if each_thread+1 == process_thread:
                        process_df = whole_df[initial : ]
                    else:
                        process_df = whole_df[initial : track]
                        
                    thread_x = threading.Thread(target = processor, args = (process_df, outfile, tidx, each_thread, total_samples))
                    thread_x.start()
                    thread_pool.append(thread_x)
                    
                    initial = track
                    track += split_count
    
            else:
                processor(whole_df, outfile, tidx, 0, total_samples)
                thread_out = f'{outfile}_reader{tidx}_processor{0}.vcf'
                with open(thread_out, 'r') as fh:
                    for line in fh:
                        repeat_info = line.strip().split('\t')
                        out.write("\t".join(map(str, repeat_info)) + "\n")
                os.remove(thread_out)
    
            del whole_df

    if thread_pool:
        # joining previous threads - waiting for previous threads to be over
        for thread_x in thread_pool:
            # print('waiting for ', thread_x)
            thread_x.join()
        thread_pool.clear()
    
        # print('Concatenating processor files..............')
        for each_thread in range(process_thread):
            thread_out = f'{outfile}_reader{tidx}_processor{each_thread}.vcf'
            # print('opening ', thread_out)
            with open(thread_out, 'r') as fh:
                for line in fh:
                    repeat_info = line.strip().split('\t')
                    out.write("\t".join(map(str, repeat_info)) + "\n")
            # print('Removing ', thread_out)
            os.remove(thread_out)
    
    for i in vcf_instance:
        i.close()
    ref_file.close()
    tbx.close()
    out.close()


def joiner(frames, parquet_batch, tidx, outfile):
    base = reduce(lambda l, r: l.join(r, on=['s', 'e'], how='left'), frames)
    df = base.collect(engine="streaming")
    df.write_parquet(f"{outfile}_reader{tidx}_batch{parquet_batch}.parquet", compression="zstd")


def vcf_writer(out, bam_name, source_vcf_path):

    source_vcf = pysam.VariantFile(source_vcf_path)

    vcf_header = pysam.VariantHeader()

    # command
    vcf_header.add_line(f"##command=Tamatr {' '.join(sys.argv)}")

    # print(source_vcf)

    for contig, metadata in source_vcf.header.contigs.items():
        vcf_header.contigs.add(contig, length=metadata.length)
    source_vcf.close()
    
    #sample_name
    for each_sample in bam_name:
        vcf_header.add_sample(each_sample)
    # FILTER
    vcf_header.filters.add('LESS_READS', number=None, type=None, description="Read depth below threshold")
    # INFO
    vcf_header.info.add("AC", number='A', type="Integer", description="Number of alternate alleles in called genotypes")
    vcf_header.info.add("AN", number=1, type="Integer", description="Number of alleles in called genotypes")
    vcf_header.info.add("MOTIF", number=1, type="String", description="Repeat motif")
    vcf_header.info.add("START", number=1, type="Integer", description="Start position of the repeat region in 0-based coordinate system")
    vcf_header.info.add("END", number=1, type="Integer", description="End position of the repeat region")
    vcf_header.info.add("ID", number=1, type="String", description="Locus identifier tag")
    vcf_header.info.add("CT", number=1, type="String", description="Cluster type")
    vcf_header.info.add("EAC", number=1, type="String", description="Each Allele Count")
    # FORMAT
    vcf_header.formats.add("GT", number=1, type="String", description="Genotype")
    vcf_header.formats.add("AL", number=2, type="Integer", description="Allele length in base pairs")
    vcf_header.formats.add("AR", number='.', type="String", description="Allele length range")
    vcf_header.formats.add("SD", number='.', type="Integer", description="Number of reads supporting for the alleles")
    vcf_header.formats.add("PC", number=2, type="Integer", description="Number of reads in the phased cluster for each allele")
    vcf_header.formats.add("DP", number=1, type="Integer", description="Number of the supporting reads for the repeat locus")
    vcf_header.formats.add("SN", number='.', type="Integer", description="Number of SNPs used for phasing")
    vcf_header.formats.add("SQ", number='.', type="Float", description="Phred-scale qualities of the SNPs used for phasing")
    vcf_header.formats.add("MM", number='.', type="Float", description="Mean methylation level for each allele")
    vcf_header.formats.add("MR", number='.', type="Integer", description="Number of reads providing methylation info for each allele")
    vcf_header.formats.add("DS", number='A', type="String", description="Motif decomposed sequence")
    vcf_header.formats.add("MV", number='.', type="String", description="Visual methylation encodings for the alleles")

    out.write(str(vcf_header))