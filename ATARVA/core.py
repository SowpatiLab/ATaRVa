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

from ATARVA.version import __version__
from ATARVA.baseline import *

def parse_args():
    """
    Parse command line arguments.
    """
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
    
    optional.add_argument('-o',  '--vcf', type=str, metavar='<FILE>', default='', help='name of the output file, output is in vcf format.')
    
    optional.add_argument('-q', '--map-qual', type=int, metavar='<INT>', default=5, help='minimum mapping quality of the reads to be considered. [default: 5]')
    optional.add_argument('--contigs', nargs='+', help='contigs to get genotyped [chr1 chr12 chr22 ..]. If not mentioned every contigs in the region file will be genotyped.')
    optional.add_argument('--min-reads', type=int, metavar='<INT>', default=10, help='minimum read coverage after quality cutoff at a locus to be genotyped. [default: 10]')
    optional.add_argument('--max-reads', type=int, metavar='<INT>', default=300, help='maximum number of reads to be used for genotyping a locus. [default: 100]')
    optional.add_argument('--flank', type=int, metavar='<INT>', default=10, help='length of the flanking region (in base pairs) to search for insertion with a repeat in it. [default: 10]')
    
    optional.add_argument('--snp-dist', type=int, metavar='<INT>', default=5000, help='maximum distance of the SNP from repeat region to be considered for phasing. [default: 5000]')
    optional.add_argument('--snp-count', type=int, metavar='<INT>', default=3, help='number of SNPs to be considered for phasing (minimum value = 1). [default: 3]')
    optional.add_argument('--snp-qual', type=int, metavar='<INT>', default=13, help='minimum basecall quality at the SNP position to be considered for phasing. [default: 13]')
    optional.add_argument('--snp-read', type=float, metavar='<FLOAT>', default=0.2, help='a positive float as the minimum fraction of snp\'s read contribution to be used for phasing. [default: 0.25]')
    
    optional.add_argument('--phasing-read', type=float, metavar='<FLOAT>', default=0.4, help='a positive float as the minimum fraction of total read contribution from the phased read clusters. [default: 0.4]')
    optional.add_argument('--platform', type=str, metavar='<STR>', default='simplex', help='sequencing platform used for generating the data. changing this will have an affect \
                                                                           on phasing which is happening on SNPs. allowed options: [hifi, duplex, simplex-hq, simplex]. \
                                                                           default: [simplex]')
    optional.add_argument('-p',  '--num-processes', type=int, metavar='<INT>', default=1, help='number of processor. [default: 1]')
    optional.add_argument('--karyotype', nargs='+', help='karyotype of the samples [XY XX]')

    optional.add_argument('-log', '--debug-mode', action='store_true', help="write the debug messages to log file. [default: False]")
    optional.add_argument('-v', '--version', action='version', version=f'ATaRVa version {__version__}')
    

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    return parser.parse_args()


class SNPParams():

    def __init__(self, snp_dist, snp_count, snp_qual, snp_read):
        self.snp_dist = snp_dist
        self.snp_count = snp_count
        self.snp_qual = snp_qual
        self.snp_read = snp_read

def fasta_check(path):
    """
    Check if the provided FASTA file is valid.
    
    Args:
        path (str): Path to the FASTA file.
    
    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file is not a valid FASTA file.
        OSError: If there is an error reading the file.
    """
    try:
        f = pysam.FastaFile(path)
        f.close()
    except (FileNotFoundError, ValueError, OSError) as e:
        print(f"Error: {path} is not a valid FASTA file. {str(e)}")
        sys.exit()
    except Exception as e:
        print("An unexpected error occurred:", str(e))
        sys.exit()


def bam_check(path, aln_format):
    """
    Check if the provided BAM file is valid and sorted by coordinate.
    
    Args:
        path (str): Path to the BAM file.
        aln_format (str): Format of the alignment file.
    
    Raises:
        FileNotFoundError: If the file does not exist.
        SortOrderError: If the BAM file is not sorted by coordinate.
        ValueError: If the file is not a valid BAM file.
        OSError: If there is an error reading the file.
    """

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


def tabix_check(path):
    """
    Check if the provided regions file is valid and indexed.
    
    Args:
        path (str): Path to the regions file.
        
    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file is not a valid tabix file.
        OSError: If there is an error reading the file.
    """

    try:
        t = pysam.TabixFile(path)
        t.close()
    except (FileNotFoundError, ValueError, OSError) as e:
        print(f"Error: {path} is not a valid tabix file. {str(e)}")
        sys.exit()
    except Exception as e:
        print("An unexpected error occurred:", str(e))
        sys.exit()


def splitfile_threads(tbx, contigs, total_loci, threads):
    split_point = total_loci // threads
    # split_point is 0 when the total_loci is less than the number of threads
    # this is a rare case with a bed file with very few loci; all these loci will be handled by a single thread
    if split_point == 0: split_point = total_loci
    
    fetcher = []
    line_count = 0
    current_split = []
    for each_contig in contigs:
        init = 0
        for row in tbx.fetch(each_contig):
            line_count += 1
            if init == 0:
                fields = row.split('\t')
                chrom  = fields[0]
                start_coord = (int(fields[1]), int(fields[2]))
                init = 1
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

    return fetcher


def main():

    start_time = ti.default_timer()
    parameters = parse_args()

    for arg in vars(parameters):
        print (arg, getattr(parameters, arg))
    print('\n')

    # checking the input alignment format. default is bam
    if    parameters.format == 'cram': parameters.format = 'rc'
    elif  parameters.format == 'sam':  parameters.format = 'r'
    else: parameters.format = 'rb'

    # checking the input files
    fasta_check(parameters.fasta)
    for each_bam in parameters.bams:
        bam_check(each_bam, parameters.format)
    tabix_check(parameters.regions)

    parameters.vcf = sys.std.out
    if parameters.vcf:
        # cleaning the output file name
        if '.vcf' == parameters.vcf[-4:]:
            parameters.vcf = f'{parameters.vcf}'[:-4]
        elif parameters.vcf[-1]=='/':             # ?? I don't understand this condition
            parameters.vcf = parameters.vcf + "atarva"

    tbx  = pysam.Tabixfile(parameters.regions)

    # the contigs to be genotyped
    total_loci = 0
    if not parameters.contigs:
        parameters.contigs = sorted(tbx.contigs)
        for row in tbx.fetch(): total_loci += 1
    else:
        parameters.contigs = sorted(parameters.contigs)
        for each_contig in parameters.contigs:
            for row in tbx.fetch(each_contig): total_loci += 1

    if not parameters.karyotype:
        karyotype_list = [False]*len(parameters.bams)
    else:
        karyotype_list = [i=='XY' for i in parameters.karyotype]

    threads = parameters.num_processes
    fetcher = splitfile_threads(tbx, parameters.contigs, threads)

    mbso = False    # mbso - multiple bams single output
    if (len(parameters.bams)>1) and (parameters.vcf): mbso = True
    
    # each bam file is processed in separately
    for kidx, each_bam in enumerate(parameters.bams):
        print(f"Processing sample {each_bam.split('/')[-1]}\n")

        reads_sampled = 0 # sample set of reads to look for tags
        aln_file = pysam.AlignmentFile(each_bam, parameters.format)
        read_length = 0
        for read in aln_file.fetch():
            # ?? could you add a comment exmplaing the flags here
            # 0x400 - read is PCR or optical duplicate
            # 0x100 - not primary alignment
            if (read.flag & 0X400) or (read.flag & 0X100): continue 
            reads_sampled +=1
            cigar = read.cigarstring
            # ?? I think the read length and also if the data is SRS or LRS should be provided by the user
            # ?? Also should read_length be summed?
            read_length += read.query_length
            
            if read.has_tag('cs'):
                print("CS tag detected. Processing using CS tag...\n")
                break
            
            elif (cigar!=None) and (('X' in cigar) or ('=' in cigar)):
                print("CIGAR(X/=) tag detected. Processing using CIGAR(X/=) tag...\n")
                break
            
            elif read.has_tag('MD'):
                print("MD tag detected. Processing using MD tag...")
                print("Include CS tag or CIGAR tag with 'X/=' for faster processing.\n")
                break
            
            if reads_sampled > 100:
                print(f"No tags detected in {each_bam.split('/')[-1]}. Processing without tags...")
                print("Include the CS tag, MD tag, or CIGAR tag with 'X/=' for faster processing.\n")
                break
                # sys.exit()
        aln_file.close()
        
        srs = False
        # if average read length is less than 350bp consider it as short read data
        # ?? also this needs to be fixed as error model we might consider for short read and long read data is different
        if read_length/reads_sampled < 350:
            print('Short reads detected... Processing in short-read mode.')
            srs = True
        else: print('Long reads detected... Processing in long-read mode.')

        if not parameters.vcf:
            parameters.vcf = f'{".".join(each_bam.split("/")[-1].split(".")[:-1])}'
        elif mbso or (parameters.vcf[-1]=='/'):
            parameters.vcf = parameters.vcf + '_' + ".".join(each_bam.split("/")[-1].split('.')[:-1])

        if threads > 1:
            thread_pool = list()
            # initializing threads
            for tidx in range(threads):
                contig = fetcher[tidx]
                target = ''
                if srs: target = mini_cooper
                else: target = cooper
                
                thread = Process(target = target, args = (parameters, each_bam, contig, tidx, karyotype_list[kidx]))
                thread.start()
                thread_pool.append(thread)
            
            # joining Threads 
            for tidx, thread in enumerate(thread_pool):
                thread.join()
            # emptying thread_pool
            thread_pool.clear()

            out = open(f'{parameters.vcf}.vcf', 'a')
            print('Concatenating thread outputs!', file=sys.stderr)
            for tidx in range(threads)[1:]:
                thread_out = f'{parameters.vcf}_thread_{tidx}.vcf'
                with open(thread_out, 'r') as fh:
                    # if tidx!=0: next(fh)
                    for line in fh:
                        repeat_info = line.strip().split('\t')
                        print(*repeat_info, file=out, sep='\t')
                os.remove(thread_out)
            out.close()
            print('Concatenation completed!! ^_^', file=sys.stderr)

            if parameters.debug_mode:
                log_file = open(f'{parameters.vcf}_debug.log', 'a')
                for tidx in range(threads)[1:]:
                    thread_log_out = f'{parameters.vcf}_debug_{tidx}.log'
                    with open(thread_log_out, 'r') as fh:
                        for line in fh:
                            log_info = line.strip()
                            print(log_info, file=log_file)
                    os.remove(thread_log_out)
                log_file.close()
        else:
            if srs:
                mini_cooper(parameters, each_bam, contig, tidx, karyotype_list[kidx])
            else:
                cooper(parameters) # ?? what is purpose of this?
                cooper(parameters, each_bam, contig, tidx, karyotype_list[kidx])

    time_now = ti.default_timer()
    sys.stderr.write('CPU time: {} seconds\n'.format(time_now - start_time))

if __name__ == '__main__':
    main()
