#!/usr/bin/env python

"""
    ATaRVa (pronunced as atharva) is a tool designed to analyse tandem repeat variation
    from long/short read whole genome sequencing data.
"""

import sys, os
import pysam
import timeit as ti
import statistics
import argparse as ap
from multiprocessing import Process
from version import __version__

from baseline import cooper

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
    optional.add_argument('--level-split', type=int, metavar='<INT>', default=2, help='a positive integer(0 to 2, where 0 : 30 to 70%% ; 1 : 25 to 75%% ; 2 : 20 to 80%%) as the percentage level of read split of snps to be used for phasing. [default: 2]')
    optional.add_argument('--snp-read', type=float, metavar='<FLOAT>', default=0.2, help='a positive float as the minimum fraction of snp\'s read contribution to be used for phasing. [default: 0.25]')
    optional.add_argument('--phasing-read', type=float, metavar='<FLOAT>', default=0.4, help='a positive float as the minimum fraction of total read contribution from the phased read clusters. [default: 0.4]')
    optional.add_argument('-o',  '--vcf', type=str, metavar='<FILE>', default='', help='name of the output file, output is in vcf format. [default: sys.stdout]')
    optional.add_argument('--platform', type=str, metavar='<STR>', default='simplex', help='sequencing platform used for generating the data. changing this will have an affect \
                                                                           on phasing which is happening on SNPs. allowed options: [hifi, duplex, simplex-hq, simplex]. \
                                                                           default: [simplex]')
    optional.add_argument('-p',  '--processor', type=int, metavar='<INT>', default=1, help='number of processor. [default: 1]')
    optional.add_argument('-v', '--version', action='version', version=f'ATaRVa version {__version__}')

    

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    return parser.parse_args()



def main():

    start_time = ti.default_timer()
    args = parse_args()
    
    for arg in vars(args):
        print (arg, getattr(args, arg))

    aln_format = ''         # format of the alignment file
    if   args.format == 'cram': aln_format = 'rc'
    elif args.format == 'sam':  aln_format = 'r'
    else:            aln_format = 'rb'

    seq_platform = args.platform     # seq_tech of the alignment file


    out_file = sys.stdout
    if args.vcf:
        if '.vcf' == args.vcf[-4:]:
            out_file = f'{args.vcf}'[:-4]
        else:
            out_file = f'{args.vcf}'
    # else: out_file = f'{".".join(args.bams.split(".")[:-1])}'

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

        if not args.vcf:
            out_file = f'{".".join(each_bam.split(".")[:-1])}'
        elif mbso:
            out_file = ".".join(each_bam.split(".")[:-1]) + out_file

        if threads > 1:
            thread_pool = list()
            # initializing threads
            for tidx in range(threads):
                contig = fetcher[tidx]
                thread_x = Process(
                    target = cooper,
                    args = (each_bam, args.regions, args.fasta, aln_format, contig, args.map_qual, out_file, seq_platform, args.level_split, args.snp_qual, args.snp_count, args.snp_dist, args.max_reads, args.min_reads, args.snp_read, args.phasing_read, tidx))
                thread_x.start()
                thread_pool.append(thread_x)
            # joining Threads 
            for tidx, thread_x in enumerate(thread_pool):
                thread_x.join()
            # emptying thread_pool
            thread_pool.clear()
        
            out = open(f'{out_file}.vcf', 'w')
            print('Concatenating thread outputs!', file=sys.stderr)
            for tidx in range(threads):
                thread_out = f'{out_file}_thread_{tidx}.vcf'
                with open(thread_out, 'r') as fh:
                    if tidx!=0: next(fh)
                    for line in fh:
                        repeat_info = line.strip().split('\t')
                        print(*repeat_info, file=out, sep='\t')
                os.remove(thread_out)
            out.close()
            print('Concatenation completed!! ^_^', file=sys.stderr)
        else:
            cooper(each_bam, args.regions, args.fasta, aln_format, fetcher[0], args.map_qual, out_file, seq_platform, args.level_split, args.snp_qual, args.snp_count, args.snp_dist, args.max_reads, args.min_reads, args.snp_read, args.phasing_read, -1)
                    

    time_now = ti.default_timer()
    sys.stderr.write('CPU time: {} seconds\n'.format(time_now - start_time))

if __name__ == '__main__':
    main()