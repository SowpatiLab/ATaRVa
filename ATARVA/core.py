#!/usr/bin/env python

"""
    ATaRVa (pronunced as atharva) is a tool designed to analyse tandem repeat variation
    from long/short read whole genome sequencing data.
"""

import sys, os
import pysam
import timeit as ti
import argparse as ap
from multiprocessing import Process

from version import __version__
from baseline import *

def parse_args():
    parser = ap.ArgumentParser()
    parser._action_groups.pop()

    print("ATaRVa (atharva) - Analysis of Tandem Repeat Variants\nSowpati Lab\n")

    required = parser.add_argument_group('Required arguments')
    required.add_argument('-fi',  '--fasta',   required=True, metavar='<FILE>', help='input reference fasta file')
    required.add_argument('--bams', nargs='+', required=True, metavar='<FILE>', help='samples alignment files. allowed formats: SAM, BAM, CRAM')
    required.add_argument('-bed', '--regions', required=True, metavar='<FILE>', help='input regions file. the regions file should be strictly in bgzipped tabix format. \
                                                                  If the regions input file is in bed format. First sort it using bedtools. Compress it using bgzip. \
                                                                  Index the bgzipped file with tabix command from samtools package.')

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--format', type=str, metavar='<STR>', default='bam', help='format of input alignment file. allowed options: [cram, bam, sam]. default: [bam]')
    optional.add_argument('-q', '--map-qual', type=int, metavar='<INT>', default=5, help='minimum mapping quality of the reads to be considered. [default: 5]')
    optional.add_argument('--contigs', nargs='+', help='contigs to get genotyped [chr1 chr12 chr22 ..]. If not mentioned every contigs in the region file will be genotyped.')
    optional.add_argument('--min-reads', type=int, metavar='<INT>', default=10, help='minimum read coverage after quality cutoff at a locus to be genotyped. [default: 10]')
    optional.add_argument('--max-reads', type=int, metavar='<INT>', default=300, help='maximum number of reads to be used for genotyping a locus. [default: 100]')
    optional.add_argument('--snp-dist', type=int, metavar='<INT>', default=5000, help='maximum distance of the SNP from repeat region to be considered for phasing. [default: 5000]')
    optional.add_argument('--snp-count', type=int, metavar='<INT>', default=3, help='number of SNPs to be considered for phasing (minimum value = 1). [default: 3]')
    optional.add_argument('--snp-qual', type=int, metavar='<INT>', default=13, help='minimum basecall quality at the SNP position to be considered for phasing. [default: 13]')
    # optional.add_argument('--level-split', type=int, metavar='<INT>', default=2, help='a positive integer(0 to 2, where 0 : 30 to 70%% ; 1 : 25 to 75%% ; 2 : 20 to 80%%) as the percentage level of read split of snps to be used for phasing. [default: 2]')
    optional.add_argument('--snp-read', type=float, metavar='<FLOAT>', default=0.2, help='a positive float as the minimum fraction of snp\'s read contribution to be used for phasing. [default: 0.25]')
    optional.add_argument('--phasing-read', type=float, metavar='<FLOAT>', default=0.4, help='a positive float as the minimum fraction of total read contribution from the phased read clusters. [default: 0.4]')
    optional.add_argument('-o',  '--vcf', type=str, metavar='<FILE>', default='', help='name of the output file, output is in vcf format. [default: sys.stdout]')
    optional.add_argument('--platform', type=str, metavar='<STR>', default='simplex', help='sequencing platform used for generating the data. changing this will have an affect \
                                                                           on phasing which is happening on SNPs. allowed options: [hifi, duplex, simplex-hq, simplex]. \
                                                                           default: [simplex]')
    optional.add_argument('-p',  '--processor', type=int, metavar='<INT>', default=1, help='number of processor. [default: 1]')
    optional.add_argument('-v', '--version', action='version', version=f'ATaRVa version {__version__}')
    optional.add_argument('-log', '--debug_mode', action='store_true', help="write the debug messages to log file. [default: False]")

    

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    return parser.parse_args()


def f_check(path):
    try:
        f = pysam.FastaFile(path)
        f.close()
    except (FileNotFoundError, ValueError, OSError) as e:
        print(f"Error: {path} is not a valid FASTA file. {str(e)}")
        sys.exit()
    except Exception as e:
        print("An unexpected error occurred:", str(e))
        sys.exit()

def b_check(path, aln_format):
    try:
        b = pysam.AlignmentFile(path, aln_format)
        header = b.header
        if 'HD' in header and 'SO' in header['HD']:
            sort_order = header['HD']['SO']
            if sort_order == 'coordinate':
                pass
                # print(f"Alignment file sort order: {sort_order}")
            else:
                print(f"Alignment file sort order: {sort_order}. It should be sorted by \'coordinate\'!!")
                print(f"Use: samtools sort sorted_{path.split('/')[-1]} {path.split('/')[-1]}")
                sys.exit()
        else:
            print("No sort order specified in the header.")
            print(f"Use: samtools sort sorted_{path.split('/')[-1]} {path.split('/')[-1]}")
            sys.exit()

        b.close()
    except (FileNotFoundError, ValueError, OSError) as e:
        print(f"Error: {path} is not a valid alignment file. {str(e)}")
        sys.exit()
    except Exception as e:
        print("An unexpected error occurred:", str(e))
        sys.exit()

def t_check(path):
    try:
        t = pysam.TabixFile(path)
        t.close()
    except (FileNotFoundError, ValueError, OSError) as e:
        print(f"Error: {path} is not a valid tabix file. {str(e)}")
        sys.exit()
    except Exception as e:
        print("An unexpected error occurred:", str(e))
        sys.exit()

def main():

    start_time = ti.default_timer()
    args = parse_args()

    for arg in vars(args):
        print (arg, getattr(args, arg))
    print('\n')

    f_check(args.fasta)

    aln_format = ''         # format of the alignment file
    if   args.format == 'cram': aln_format = 'rc'
    elif args.format == 'sam':  aln_format = 'r'
    else:            aln_format = 'rb'

    for each_bam in args.bams:
        b_check(each_bam, aln_format)
    t_check(args.regions)
    

    seq_platform = args.platform     # seq_tech of the alignment file


    out_file = sys.stdout
    if args.vcf:
        if '.vcf' == args.vcf[-4:]:
            out_file = f'{args.vcf}'[:-4]
        elif args.vcf[-1]=='/':
            out_file = args.vcf + "atarva_tr"
        else:
            out_file = f'{args.vcf}'
    # else: out_file = f'{".".join(args.bams.split(".")[:-1])}'
    external_name = out_file

    tbx  = pysam.Tabixfile(args.regions)
    total_loci = 0
    if not args.contigs:
        contigs = sorted(tbx.contigs)
        for row in tbx.fetch():
            total_loci += 1
    else:
        contigs = sorted(args.contigs)
        for each_contig in args.contigs:
            for row in tbx.fetch(each_contig):
                total_loci += 1

    threads = args.processor
    split_point = total_loci // threads
    if split_point == 0:
        split_point = 1

    
    fetcher = []
    line_count = 0
    current_split = []
    for each_contig in contigs:
        init = 0
        for row in tbx.fetch(each_contig):
            line_count += 1
            if init == 0:
                Row = row.split('\t')
                chrom = Row[0]
                start_coord = (int(Row[1]), int(Row[2]))
                init=1
            if len(fetcher) < threads-1:
                if line_count % split_point == 0:
                    end_coord = (int(row.split('\t')[1]), int(row.split('\t')[2]))
                    current_split.append([chrom, start_coord, end_coord])
                    fetcher.append(tuple(current_split))
                    line_count = 0
                    current_split = []
                    init = 0
        if init != 0:
            end_coord = (int(row.split('\t')[1]), int(row.split('\t')[2]))
            current_split.append([chrom, start_coord, end_coord])
    fetcher.append(tuple(current_split))
    tbx.close()

    mbso = 0
    if (len(args.bams)>1) and (args.vcf):
        mbso = 1
    
    for each_bam in args.bams:
        out_file = external_name
        print(f"Processing sample {each_bam.split('/')[-1]}\n")

        count = 0
        aln_file = pysam.AlignmentFile(each_bam, aln_format)
        length = 0
        for read in aln_file.fetch():
            if (read.flag & 0X400) or (read.flag & 0X100): continue 
            count+=1
            string = read.cigarstring
            length += read.query_length
            if read.has_tag('cs'):
                print("CS tag detected. Processing using CS tag...\n")
                break
            elif (string!=None) and (('X' in string) or ('=' in string)):
                print("CIGAR(X/=) tag detected. Processing using CIGAR(X/=) tag...\n")
                break
            elif read.has_tag('MD'):
                print("MD tag detected. Processing using MD tag...")
                print("Include CS tag or CIGAR tag with 'X/=' for faster processing.\n")
                break
            if count>100:
                print(f"No tags detected in {each_bam.split('/')[-1]}. Processing without tags...")
                print("Include the CS tag, MD tag, or CIGAR tag with 'X/=' for faster processing.\n")
                # sys.exit()
        aln_file.close()
        srs = False
        if length/count < 350:
            print('Short reads detected... Processing in short-read mode.')
            srs = True
        else: print('Long reads detected... Processing in long-read mode.')
        

        if not args.vcf:
            out_file = f'{".".join(each_bam.split("/")[-1].split(".")[:-1])}'
        elif mbso or (out_file[-1]=='/'):
            out_file = out_file + '_' + ".".join(each_bam.split("/")[-1].split('.')[:-1])

        

        if threads > 1:
            thread_pool = list()
            # initializing threads
            for tidx in range(threads):
                contig = fetcher[tidx]
                if srs:
                    thread_x = Process(
                        target = mini_cooper,
                        args = (each_bam, args.regions, args.fasta, aln_format, contig, args.map_qual, out_file, seq_platform, args.snp_qual, args.snp_count, args.snp_dist, args.max_reads, args.min_reads, args.snp_read, args.phasing_read, tidx, args.debug_mode))
                else:
                    thread_x = Process(
                        target = cooper,
                        args = (each_bam, args.regions, args.fasta, aln_format, contig, args.map_qual, out_file, seq_platform, args.snp_qual, args.snp_count, args.snp_dist, args.max_reads, args.min_reads, args.snp_read, args.phasing_read, tidx, args.debug_mode))
                thread_x.start()
                thread_pool.append(thread_x)
            # joining Threads 
            for tidx, thread_x in enumerate(thread_pool):
                thread_x.join()
            # emptying thread_pool
            thread_pool.clear()
        
            out = open(f'{out_file}.vcf', 'a')
            print('Concatenating thread outputs!', file=sys.stderr)
            for tidx in range(threads)[1:]:
                thread_out = f'{out_file}_thread_{tidx}.vcf'
                with open(thread_out, 'r') as fh:
                    # if tidx!=0: next(fh)
                    for line in fh:
                        repeat_info = line.strip().split('\t')
                        print(*repeat_info, file=out, sep='\t')
                os.remove(thread_out)
            out.close()
            print('Concatenation completed!! ^_^', file=sys.stderr)

            if args.debug_mode:
                out_log = open(f'{out_file}_debug.log', 'a')
                for tidx in range(threads)[1:]:
                    thread_log_out = f'{out_file}_debug_{tidx}.log'
                    with open(thread_log_out, 'r') as fh:
                        for line in fh:
                            log_info = line.strip()
                            print(log_info, file=out_log)
                    os.remove(thread_log_out)
                out_log.close()
        else:
            if srs:
                mini_cooper(each_bam, args.regions, args.fasta, aln_format, fetcher[0], args.map_qual, out_file, seq_platform, args.snp_qual, args.snp_count, args.snp_dist, args.max_reads, args.min_reads, args.snp_read, args.phasing_read, -1, args.debug_mode)
            else:
                cooper(each_bam, args.regions, args.fasta, aln_format, fetcher[0], args.map_qual, out_file, seq_platform, args.snp_qual, args.snp_count, args.snp_dist, args.max_reads, args.min_reads, args.snp_read, args.phasing_read, -1, args.debug_mode)

    time_now = ti.default_timer()
    sys.stderr.write('CPU time: {} seconds\n'.format(time_now - start_time))

if __name__ == '__main__':
    main()