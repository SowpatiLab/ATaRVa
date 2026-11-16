import subprocess
import sys
from pathlib import Path

def extract_vcf_fields(vcf_file):
    """
    Extracts the CHROM, POS, REF, ALT, Allele, Meth and Meth_tag   fields from a VCF file.

    Parameters:
        vcf_file (str): Path to the VCF file.
    """
    extracted_data = {}
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue  # Skip header lines
            columns = line.strip().split('\t')
            if len(columns) < 8:
                continue  # Skip malformed lines
            chrom = columns[0]
            ref = columns[3]
            alt = columns[4]
            info_field = columns[7].split(';')
            start = info_field[3].split('=')[1]
            end = info_field[4].split('=')[1] 
            format_field = columns[8].split(':')
            sample_fields = columns[9].split(':')

            # Extract Allele, Meth and Meth_tag from the SAMPLE field
            sample_dict = dict(zip(format_field, sample_fields))
            allele = [int(i) for i in sample_dict.get('AL').split(',')] if '.' not in sample_dict.get('AL') else []
            meth = [float(i) for i in sample_dict.get('MA').split(',')] if '0' in sample_dict.get('MA') else []
            meth_tag = sample_dict.get('MV').split(',') if '.' not in sample_dict.get('MV') else []

            extracted_data[f'{chrom}:{start}-{end}'] = (allele, meth, meth_tag)
    return extracted_data

expected_data = {'chr4:1000-1093' : ([108, 87], [0.0, 0.01], ['AAAAAAAA', 'AEEAAAAA']), # CAG - HD_HTT HG002 sample
                 'chrX:1000-1062' : ([95, 95], [0.01, 0.01], ['AAAAADADAAAAAAADAADAAAADAADAA', 'AAAAADADAAAAAAADAADAAAADAADAA']), # CGG - FMR1 HG002 sample
                 }
def assert_vcf_fields(vcf_file, expected_data):
    """
    Asserts that the extracted fields from a VCF file match the expected data.

    Parameters:
        vcf_file (str): Path to the VCF file.
        expected_data (dict): Expected data in the format {chrom_pos_end: (allele, meth, meth_tag)}.
    """
    extracted_data = extract_vcf_fields(vcf_file)
    assert extracted_data.keys() == expected_data.keys(), (
    f"Missing loci: {expected_data.keys() - extracted_data.keys()}\n"
    f"Unexpected loci: {extracted_data.keys() - expected_data.keys()}"
    )
    for locus_key in extracted_data:
        extracted_allele, extracted_meth, extracted_meth_tag = extracted_data[locus_key]
        expected_allele, expected_meth, expected_meth_tag = expected_data[locus_key]

        assert sorted(extracted_allele) == sorted(expected_allele), f"Allele mismatch for {locus_key}: expected {expected_allele}, got {extracted_allele}"
        assert sorted(extracted_meth) == sorted(expected_meth), f"Meth mismatch for {locus_key}: expected {expected_meth}, got {extracted_meth}"
        assert sorted(extracted_meth_tag) == sorted(expected_meth_tag), f"Meth_tag mismatch for {locus_key}: expected {expected_meth_tag}, got {extracted_meth_tag}"

def test_complete_pipeline(tmp_path):

    test_dir = Path(__file__).parent.parent / "tests/test_data"

    # Input BAM file
    input_bam = test_dir / "test_reads.bam"
    input_bed = test_dir / "test_region.bed.gz"
    input_fasta = test_dir / "test_ref.fasta"

    # Temporary output file
    output_vcf = tmp_path / "tmp_test_output.vcf"

    # Execute the complete tool
    result = subprocess.run(
        [
            sys.executable,
            "ATARVA/core.py",
            "genotype",
            "-b", str(input_bam), "-r", str(input_bed),
            "-f", str(input_fasta), "-o", str(output_vcf),
            "--decompose", "--methviz"
        ],
        capture_output=True,
        text=True
    )

    # Check successful execution
    assert result.returncode == 0, (
        f"Pipeline failed!\nSTDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )

    # Check output file
    assert output_vcf.exists(), "Output VCF was not generated"
    assert output_vcf.stat().st_size > 0, "Output VCF is empty"

    # Validate extracted VCF fields
    assert_vcf_fields(output_vcf, expected_data)