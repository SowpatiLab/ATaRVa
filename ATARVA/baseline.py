from ATARVA.locus_utils import process_locus
from ATARVA.cstag_utils import parse_cstag
from ATARVA.cigar_utils import parse_cigar_tag
from ATARVA.operation_utils import update_homopolymer_coords
from ATARVA.genotype_utils import analyse_genotype
from ATARVA.vcf_writer import *
from ATARVA.consensus import consensus_seq_poa

from tqdm import tqdm
import pysam
import sys
import logging

def locus_processor(global_loci_keys, global_loci_ends, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minR, ref, Chrom, global_loci_info, out, snpQ, snpC, snpD, snpR, phasingR, tbx, flank, sorted_global_ins_rpos_set, log_bool, logger, male, prev_locus_end, decomp, hp_code):

    genotyped_loci = 0
    popped    = global_loci_ends.pop(0)
    locus_key = global_loci_keys.pop(0)
    lstart = int(locus_key[locus_key.index(':')+1 : locus_key.index('-')])
    lend = int(locus_key[locus_key.index('-')+1:])
    near_by_loci = []
    for row in tbx.fetch(Chrom, lstart-flank, lend+flank):
        row = row.split('\t')
        near_by_loci.append( ( int(row[1]), int(row[2]) ) )

    if locus_key in global_loci_variations:
        

        if not sorted_global_snp_list:
            sorted_global_snp_list = sorted(global_snp_positions.keys())
        

        prev_reads, category, homozygous_allele, reads_of_homozygous, hallele_counter, skip_point, max_limit, haplotypes = process_locus(locus_key, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minR, global_loci_info, near_by_loci, sorted_global_ins_rpos_set, Chrom, lstart, lend, ref, log_bool, logger, snpD, prev_locus_end, hp_code)

        read_seqs = global_loci_variations[locus_key]['read_sequence']
        if category == 1:
            ref_allele_length = lend - lstart
            if homozygous_allele != ref_allele_length:
                seqs = [seq for seq in [read_seqs[read_id][0] for read_id in reads_of_homozygous] if seq!='']
                if len(seqs)>0:
                    ALT = consensus_seq_poa(seqs)
                    homozygous_allele = len(ALT)
                else:
                    ALT = '<DEL>'
                    homozygous_allele = 0
            else:
                ALT = '.'
            vcf_homozygous_writer(ref, Chrom, locus_key, global_loci_info, homozygous_allele, global_loci_variations, len(reads_of_homozygous), out, ALT, log_bool, '.', decomp)
            genotyped_loci += 1
        elif category == 2:
            state, skip_point = analyse_genotype(Chrom, locus_key, global_loci_info, global_loci_variations, global_read_variations, global_snp_positions, hallele_counter, ref, out, sorted_global_snp_list, snpQ, snpC, snpD, snpR, phasingR, maxR, max_limit, male, log_bool, decomp)
            if state: genotyped_loci += 1
            else:
                skip_messages = {
                    0: 'Locus skipped due to insignificant snps at the level of read split.',
                    1: 'Locus skipped due to less read contribution of Significant snps.',
                    2: 'Locus skipped due to less read contribution in the phased clusters.',
                    6: f'Locus {locus_key} skipped due to wide distribution of alleles with one read supporting to it.',
                }
                tqdm.write(skip_messages.get(skip_point, 'Locus skipped due to less number of significant snps based on user\'s parameter.'))
        elif category == 3:
            genotypes = []
            allele_count = {}
            ALT_seqs = []
            phased_read = []
            for hap_reads in haplotypes:
                phased_read.append(len(hap_reads))
                seqs = [seq for seq in [read_seqs[read_id][0] for read_id in hap_reads] if seq!='']
                if len(seqs)>0:
                    ALT = consensus_seq_poa(seqs)
                    allele_length = len(ALT)
                else:
                    ALT = '<DEL>'
                    allele_length = 0

                ALT_seqs.append(ALT)
                genotypes.append(allele_length)

                if allele_length not in allele_count:
                    allele_count[allele_length] = len(hap_reads)
                else:
                    allele_count[str(allele_length)] = len(hap_reads)
            vcf_heterozygous_writer(Chrom, genotypes, lstart, global_loci_variations, lend, allele_count, len(reads_of_homozygous), global_loci_info, ref, out, '.', phased_read, 0, ALT_seqs, log_bool, 'HP', decomp)
            genotyped_loci += 1
        else:
            if skip_point == 0:
                # print('Locus skipped due to minimum number of supporting reads')
                vcf_fail_writer(Chrom, locus_key, global_loci_info, ref, out, len(prev_reads), skip_point)
            # else:
            #     print('Locus skipped due to maximum number of supporting reads')
                
        del global_loci_variations[locus_key]
        prev_locus_end = popped
        
    return genotyped_loci, prev_locus_end


def cooper(bam_file, tbx_file, ref_file, aln_format, contigs, mapq_threshold, outfile, snpQ, snpC, snpD, maxR, minR, snpR, phasingR, tidx, flank, log_bool, karyotype, decomp, hp_code):
    # this function iterates through each contig and processes the genotypes for each locus
    tbx  = pysam.Tabixfile(tbx_file)
    bam  = pysam.AlignmentFile(bam_file, aln_format)
    ref  = pysam.FastaFile(ref_file)

    logger = None

    # Determine the output file and log file names based on tidx
    if tidx == -1 or tidx == 0:
        out_filename = f'{outfile}.vcf'
        log_name = f'{outfile}_debug.log'
    else:
        idx = outfile.rfind('/')
        hid_outfile = outfile[:idx+1] + '.' + outfile[idx+1:]
        out_filename = f'{hid_outfile}_thread_{tidx}.vcf'
        log_name = f'{hid_outfile}_debug_{tidx}.log'
    
    # Open the output file
    out = open(out_filename, 'w')
    
    # Write the VCF header if tidx is -1 or 0
    if tidx == -1 or tidx == 0:
        vcf_writer(out, bam, bam_file.split("/")[-1].split('.')[0])
    
    # Initialize the logger if log_bool is True
    if log_bool:
        with open(log_name, 'w'):
            pass
        logging.basicConfig(
            filename=log_name,
            level=logging.DEBUG,
            format='%(levelname)s - %(message)s'
        )
        logger = logging.getLogger("MyLogger")
    
    dwrite = tidx!=-1

    whole_loci = 0
    tot_loci_list = []
    for contig in contigs:
        Chrom, Start, End = contig

        total_loci = 0
        for row in tbx.fetch(Chrom, Start[0], End[1]):
            row = row.split('\t')
            if (total_loci == 0) and (int(row[2]) != Start[1]): continue
            total_loci += 1
            if End[0] == int(row[1]):
                tot_loci_list.append(total_loci)
                whole_loci += total_loci
                break
                    
    progress_bar = tqdm(total= whole_loci, disable= dwrite, desc="Processing ", ascii="_>", ncols=75, bar_format="{l_bar}{bar}{n_fmt}/{total_fmt}")

    for cidx,contig in enumerate(contigs):

        Chrom, Start, End = contig
        male = (Chrom in {'chrX', 'chrY', 'X', 'Y'}) and karyotype
        end_coord = End[1]
        prev_locus_end = 0

        genotyped_loci_count = 0

        if not dwrite: tqdm.write(f"> {Chrom} {Start} {End} Total loci =  {tot_loci_list[cidx]}")


        global_snp_positions = dict()       # tracking the encountered SNPs
        global_read_variations = {}         # tracking the variations on each read
        global_loci_variations = {}         # tracking the variation for each locus
        global_loci_info = {}               # saving the information of each loci

        # tracking the loci
        global_loci_ends = []; global_loci_keys    = []        
        global_read_ends = []; global_read_indices = []

        prev_reads = set()
        sorted_global_snp_list = []
        sorted_global_ins_rpos_set = set()

        read_index = 0

        for read in bam.fetch(Chrom, Start[0], End[1]):
        

            # skip read with low mapping quality
            if read.mapping_quality < mapq_threshold:
                continue

            read_chrom = read.reference_name
            read_start = read.reference_start
            read_end   = read.reference_end

            while global_loci_ends and read_start > global_loci_ends[0]:

                genotype_status = locus_processor(global_loci_keys, global_loci_ends, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minR, ref, Chrom, global_loci_info, out, snpQ, snpC, snpD, snpR, phasingR, tbx, flank, sorted_global_ins_rpos_set, log_bool, logger, male, prev_locus_end, decomp, hp_code)
                genotyped_loci_count += genotype_status[0]
                prev_locus_end = genotype_status[1]
                progress_bar.update(1)
                

            while global_read_ends and read_start > global_read_ends[0]:
                # if the read is beyond the end of the first read that was tracked
                if global_loci_ends and global_read_ends[0] > global_loci_ends[0]:
                    # if the initial read useful for the first locus being tracked then it is retained
                    break
                else:

                    # remove the read information if the current read is beyond the first read and the locus
                    popped = global_read_ends.pop(0)
                    rindex = global_read_indices.pop(0)
                    if rindex in global_read_variations:
                        for pos in global_read_variations[rindex]['snps']:
                            if pos in global_snp_positions:
                                global_snp_positions[pos]['cov'] -= 1
                                
                        del_snps = [pos for pos in global_snp_positions if global_snp_positions[pos]['cov'] == 0]
                        for snp in del_snps:
                            del global_snp_positions[snp]
                            sorted_global_snp_list.remove(snp)
                        del global_read_variations[rindex]
                        del del_snps

                        if rindex in prev_reads: prev_reads.remove(rindex)
                            
                    del_ins_pos_idx = 0
                    list_rpos = sorted(sorted_global_ins_rpos_set)
                    for i in list_rpos:
                        del_ins_pos_idx+=1
                        if i > popped: break
                    del list_rpos[:del_ins_pos_idx]
                    sorted_global_ins_rpos_set = set(list_rpos)
                    


            # if the read is beyond the last locus in the bed file the loop stops
            if read_start > end_coord:
                while global_loci_ends:
                    genotype_status = locus_processor(global_loci_keys, global_loci_ends, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minR, ref, Chrom, global_loci_info, out, snpQ, snpC, snpD, snpR, phasingR, tbx, flank, sorted_global_ins_rpos_set, log_bool, logger, male, prev_locus_end, decomp, hp_code)
                    genotyped_loci_count += genotype_status[0]
                    prev_locus_end = genotype_status[1]
                    progress_bar.update(1)
                # process the loci left in global_loci_variation
                break

            # information locally saved for each read and all the loci it covers
            read_loci_variations = {}
            # set of homopolymer positions within the reference part that is covered by the read
            homopoly_positions = {}

            # repeat loci covered by the read
            loci_coords = []; loci_keys = []
            left_flank_list = []; right_flank_list = []

            for row in tbx.fetch(read_chrom, read_start, read_end):
                
                # adjust read start and end based on soft and hard clippings
                # soft and hard clippings do not consume the reference bases

                row = row.split('\t')
                locus_start = int(row[1]);  locus_end = int(row[2]); locus_len = locus_end-locus_start

                if (locus_start>=Start[0]) and (locus_end<=End[1]):
                    if locus_start==Start[0]:
                        if locus_end==Start[1]: pass
                        else: continue
                    pass
                elif locus_start<Start[0]:
                    continue
                elif locus_start>=End[0]: break
                

                # if only the read completely covers the repeat
                if ( locus_start >= read_start ) & ( locus_end <= read_end ):
                    left_flank = min(flank, locus_start - read_start)
                    right_flank = min(flank, read_end - locus_end)
                    left_flank_list.append(left_flank)
                    right_flank_list.append(right_flank)

                    loci_coords.append((locus_start - left_flank, locus_end + right_flank))

                    locus_key = f'{read_chrom}:{locus_start}-{locus_end}'
                    loci_keys.append(locus_key)
                    read_loci_variations[locus_key] = {'halen': locus_len, 'alen': locus_len, 'rlen': locus_len, 'seq': []}

                    if locus_key not in global_loci_variations:
                        global_loci_variations[locus_key] = {'rlen': locus_len, 'reads': [], 'read_allele': {}, 'read_sequence': {}, 'read_tag':[]}
                        global_loci_info[locus_key] = row

                        # adding the locus key when it is first encountered
                        global_loci_ends.append(locus_end)
                        global_loci_keys.append(locus_key)

            # if no repeats are covered by the read
            if len(loci_coords) == 0: continue


            read_index += 1
            read_quality = read.query_qualities
            cigar_tuples = read.cigartuples
            read_sequence = read.query_sequence

            cigar_one = cigar_tuples[0]

            global_read_ends.append(read_end)
            global_read_indices.append(read_index)
            global_read_variations[read_index] = {'s': read_start, 'e': read_end, 'snps': set(), 'dels': []}


            if hp_code: hp = read.has_tag(hp_code)
            else: hp = False
            if hp: hp_tag = read.get_tag(hp_code)
            else: hp_tag = None
                
            if read.has_tag('cs'):
                del cigar_tuples
                cs_tag = read.get_tag('cs')
                parse_cstag(read_index, cs_tag, read_start, loci_keys, loci_coords, read_loci_variations, homopoly_positions, global_read_variations, global_snp_positions, read_sequence, read_quality, cigar_one, sorted_global_snp_list, left_flank_list, right_flank_list, male, hp)
            else :
                parse_cigar_tag(read_index, cigar_tuples, read_start, loci_keys, loci_coords, read_loci_variations,
                                homopoly_positions, global_read_variations, global_snp_positions, read_sequence, read, ref, read_quality, sorted_global_snp_list, left_flank_list, right_flank_list, male, hp)

            for locus_key in read_loci_variations:
                global_loci_variations[locus_key]['reads'].append(read_index)
                global_loci_variations[locus_key]['read_allele'][read_index] = [read_loci_variations[locus_key]['halen'], read_loci_variations[locus_key]['alen']]
                global_loci_variations[locus_key]['read_sequence'][read_index] = read_loci_variations[locus_key]['seq']
                global_loci_variations[locus_key]['read_tag'].append(hp_tag)

        while global_loci_ends:
            genotype_status = locus_processor(global_loci_keys, global_loci_ends, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minR, ref, Chrom, global_loci_info, out, snpQ, snpC, snpD, snpR, phasingR, tbx, flank, sorted_global_ins_rpos_set, log_bool, logger, male, prev_locus_end, decomp, hp_code)
            genotyped_loci_count += genotype_status[0]
            prev_locus_end = genotype_status[1]
            progress_bar.update(1)
                
        if not dwrite: tqdm.write(f'\nTotal genotyped loci = {genotyped_loci_count} out of {tot_loci_list[cidx]} in {Chrom} {Start[0]}-{End[1]}\n')
        del global_snp_positions
        del global_read_variations
        del global_loci_variations
        del global_loci_info
        del global_loci_ends
        del global_loci_keys
        del global_read_ends
        del global_read_indices
        del prev_reads
        del sorted_global_snp_list
        del sorted_global_ins_rpos_set
        
    bam.close()
    ref.close()
    tbx.close()
    out.close()


def mini_cooper(bam_file, tbx_file, ref_file, aln_format, contigs, mapq_threshold, outfile, snpQ, snpC, snpD, maxR, minR, snpR, phasingR, tidx, flank, log_bool, karyotype, decomp, hp_code):
    # this function iterates through each contig and processes the genotypes for each locus
    tbx  = pysam.Tabixfile(tbx_file)
    tbx2  = pysam.Tabixfile(tbx_file)
    bam  = pysam.AlignmentFile(bam_file, aln_format)
    ref  = pysam.FastaFile(ref_file)


    logger = None

    # Determine the output file and log file names based on tidx
    if tidx == -1 or tidx == 0:
        out_filename = f'{outfile}.vcf'
        log_name = f'{outfile}_debug.log'
    else:
        idx = outfile.rfind('/')
        hid_outfile = outfile[:idx+1] + '.' + outfile[idx+1:]
        out_filename = f'{hid_outfile}_thread_{tidx}.vcf'
        log_name = f'{hid_outfile}_debug_{tidx}.log'
    
    # Open the output file
    out = open(out_filename, 'w')
    
    # Write the VCF header if tidx is -1 or 0
    if tidx == -1 or tidx == 0:
        vcf_writer(out, bam, bam_file.split("/")[-1].split('.')[0])
    
    # Initialize the logger if log_bool is True
    if log_bool:
        with open(log_name, 'w'):
            pass
        logging.basicConfig(
            filename=log_name,
            level=logging.DEBUG,
            format='%(levelname)s - %(message)s'
        )
        logger = logging.getLogger("MyLogger")

    dwrite = tidx!=-1

    whole_loci = 0
    tot_loci_list = []
    for contig in contigs:
        Chrom, Start, End = contig

        total_loci = 0
        for row in tbx.fetch(Chrom, Start[0], End[1]):
            row = row.split('\t')
            if (total_loci == 0) and (int(row[2]) != Start[1]): continue
            total_loci += 1
            if End[0] == int(row[1]):
                tot_loci_list.append(total_loci)
                whole_loci += total_loci
                break
                    
    progress_bar = tqdm(total= whole_loci, disable= dwrite, desc="Processing ", ascii="_>", ncols=75, bar_format="{l_bar}{bar}{n_fmt}/{total_fmt}")

    for cidx, contig in enumerate(contigs):

        Chrom, Start, End = contig
        male = (Chrom in {'chrX', 'chrY', 'X', 'Y'}) and karyotype
        prev_locus_end = 0

        genotyped_loci_count = 0

        if not dwrite: tqdm.write(f"> {Chrom} {Start} {End} Total loci =  {tot_loci_list[cidx]}")


        global_snp_positions = dict()       # tracking the encountered SNPs
        global_read_variations = {}         # tracking the variations on each read
        global_loci_variations = {}         # tracking the variation for each locus
        global_loci_info = {}               # saving the information of each loci

        # tracking the loci
        global_loci_ends = []; global_loci_keys    = []        
        global_read_ends = []; global_read_indices = []

        prev_reads = set()
        sorted_global_snp_list = []
        sorted_global_ins_rpos_set = set()

        read_index = 0

        for row in tbx.fetch(Chrom, Start[0], End[1]):

            row = row.split('\t')
            locus_start = int(row[1]);  locus_end = int(row[2]); locus_len = locus_end-locus_start


            if (locus_start>=Start[0]) and (locus_end<=End[1]):
                if locus_start==Start[0]:
                    if locus_end==Start[1]: pass
                    else: continue
                pass
            elif locus_start<Start[0]:
                continue
            elif locus_start>=End[0]: break
            
           
            for read in bam.fetch(Chrom, int(row[1]), int(row[2])):
                if read.mapping_quality < mapq_threshold:
                    continue

                read_chrom = read.reference_name
                read_start = read.reference_start
                read_end   = read.reference_end
                
                # if only the read completely covers the repeat
                if not (( locus_start >= read_start ) & ( locus_end <= read_end )):
                    continue

                
                # information locally saved for each read and all the loci it covers
                read_loci_variations = {}
                # set of homopolymer positions within the reference part that is covered by the read
                homopoly_positions = {}
    
                # repeat loci covered by the read
                loci_coords = []; loci_keys = []
 
                left_flank_list = []; right_flank_list = []

                left_flank = min(flank, locus_start - read_start)
                right_flank = min(flank, read_end - locus_end)
                left_flank_list.append(left_flank)
                right_flank_list.append(right_flank)

                loci_coords.append((locus_start - left_flank, locus_end + right_flank))

                locus_key = f'{read_chrom}:{locus_start}-{locus_end}'
                loci_keys.append(locus_key)
                read_loci_variations[locus_key] = {'halen': locus_len, 'alen': locus_len, 'rlen': locus_len, 'seq': []}

                if locus_key not in global_loci_variations:
                    global_loci_variations[locus_key] = {'rlen': locus_len, 'reads': [], 'read_allele': {}, 'read_sequence': {}, 'read_tag':[]}
                    global_loci_info[locus_key] = row
                    global_loci_ends.append(locus_end)
                    global_loci_keys.append(locus_key)

        
                read_index += 1
                read_quality = read.query_qualities
                cigar_tuples = read.cigartuples
                read_sequence = read.query_sequence

                cigar_one = cigar_tuples[0]

                global_read_ends.append(read_end)
                global_read_indices.append(read_index)
                global_read_variations[read_index] = {'s': read_start, 'e': read_end, 'snps': set(), 'dels': []}
    
                if hp_code: hp = read.has_tag(hp_code)
                else: hp = False
                if hp: hp_tag = read.get_tag(hp_code)
                else: hp_tag = None
    
                if read.has_tag('cs'):
                    del cigar_tuples
                    cs_tag = read.get_tag('cs')
                    parse_cstag(read_index, cs_tag, read_start, loci_keys, loci_coords, read_loci_variations, homopoly_positions, global_read_variations, global_snp_positions, read_sequence, read_quality, cigar_one, sorted_global_snp_list, left_flank_list, right_flank_list, male, hp)
                else :
                    parse_cigar_tag(read_index, cigar_tuples, read_start, loci_keys, loci_coords, read_loci_variations,
                                homopoly_positions, global_read_variations, global_snp_positions, read_sequence, read, ref, read_quality, sorted_global_snp_list, left_flank_list, right_flank_list, male, hp)
    
                for locus_key in read_loci_variations:
                    global_loci_variations[locus_key]['reads'].append(read_index)
                    global_loci_variations[locus_key]['read_allele'][read_index] = [read_loci_variations[locus_key]['halen'], read_loci_variations[locus_key]['alen']]
                    global_loci_variations[locus_key]['read_sequence'][read_index] = read_loci_variations[locus_key]['seq']
                    global_loci_variations[locus_key]['read_tag'].append(hp_tag)
            
            while global_loci_ends:

                genotype_status = locus_processor(global_loci_keys, global_loci_ends, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minR, ref, Chrom, global_loci_info, out, snpQ, snpC, snpD, snpR, phasingR, tbx2, flank, sorted_global_ins_rpos_set, log_bool, logger, male, prev_locus_end, decomp, hp_code)
                genotyped_loci_count += genotype_status[0]
                prev_locus_end = genotype_status[1]
                progress_bar.update(1)

            while global_read_ends:

                popped = global_read_ends.pop(0)
                rindex = global_read_indices.pop(0)
                if rindex in global_read_variations:
                    for pos in global_read_variations[rindex]['snps']:
                        if pos in global_snp_positions:
                            global_snp_positions[pos]['cov'] -= 1
                            
                    del_snps = [pos for pos in global_snp_positions if global_snp_positions[pos]['cov'] == 0]
                    for snp in del_snps:
                        del global_snp_positions[snp]
                        sorted_global_snp_list.remove(snp)
                    del global_read_variations[rindex]
                    del del_snps

                    if rindex in prev_reads: prev_reads.remove(rindex)
                        
                del_ins_pos_idx = 0
                list_rpos = sorted(sorted_global_ins_rpos_set)
                for i in list_rpos:
                    del_ins_pos_idx+=1
                    if i > popped: break
                del list_rpos[:del_ins_pos_idx]
                sorted_global_ins_rpos_set = set(list_rpos)

        if not dwrite: tqdm.write(f'\nTotal genotyped loci = {genotyped_loci_count} out of {tot_loci_list[cidx]} in {Chrom} {Start[0]}-{End[1]}\n')
        del global_snp_positions
        del global_read_variations
        del global_loci_variations
        del global_loci_info
        del global_loci_ends
        del global_loci_keys
        del global_read_ends
        del global_read_indices
        del prev_reads
        del sorted_global_snp_list
        del sorted_global_ins_rpos_set
        
    bam.close()
    ref.close()
    tbx.close()
    tbx2.close()
    out.close()