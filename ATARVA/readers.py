import pysam
import sys

def fasta_check(path):
    """
    checks file existence and format of reference fasta

    :param path: path to reference fasta file
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


def bam_check(path, aln_format, fasta):
    """
    checks file existence and format of alignment file

    :param path: path to alignment file
    :param aln_format: format of the alignment file
    :param fasta: path to reference fasta file (required to decode CRAM).
    """
    try:
        b = pysam.AlignmentFile(path, aln_format, reference_filename=fasta)
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
    checks file existence and format of tabix file

    :param path: path to tabix file
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


def check_alnfile(bamfile, args):
    """
    check the alignment file for tags and determine the genotyping mode (SRS or LRS)

    :param bamfile: path to the alignment file
    :param args: command line arguments
    :returns: srs (boolean) - True if SRS mode, False if LRS mode
    """
    n_reads   = 0
    aln_file  = pysam.AlignmentFile(bamfile, args.aln_format, reference_filename=args.fasta)
    bases     = 0
    for read in aln_file.fetch():
        if (read.flag & 0X400) or (read.flag & 0X100): continue 

        n_reads += 1
        bases   += read.query_length
        if read.has_tag('cs'):
            print("CS tag detected. Processing using CS tag...\n")
            break
        elif (read.cigarstring != None) and (('X' in read.cigarstring) or ('=' in read.cigarstring)):
            print("CIGAR(X/=) tag detected. Processing using CIGAR(X/=) tag...\n")
            break
        elif read.has_tag('MD'):
            print("MD tag detected. Processing using MD tag...")
            print("Include CS tag or CIGAR tag with 'X/=' for faster processing.\n")
            break
        if n_reads > 1000:
            print(f"No tags detected in {bamfile.split('/')[-1]}. Processing without tags...")
            print("Include the CS tag, MD tag, or CIGAR tag with 'X/=' for faster processing.\n")
            break

    aln_file.close()
    return bases/n_reads if n_reads > 0 else None
