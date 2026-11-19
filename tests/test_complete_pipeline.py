import subprocess
import sys
from pathlib import Path

def split_values(text, cast=str):
    """'0.0,.' -> [0.0]   (a '.' means 'no value for this allele' and is skipped)."""
    return [cast(v) for v in text.split(',') if v != '.']

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
            if len(columns) < 10:
                continue  # Skip malformed lines
            chrom = columns[0]
            ref = columns[3]
            alt = columns[4]
            info_field = columns[7].split(';')
            start = info_field[3].split('=')[1]
            end = info_field[4].split('=')[1] 
            format_field = columns[8].split(':')
            sample_fields = columns[9].split(':')

            # Extract Allele, Meth, DS and Meth_tag from the SAMPLE field
            sample_dict = dict(zip(format_field, sample_fields))
            allele = split_values(sample_dict['AL'], int)
            meth = split_values(sample_dict['MA'], float)
            ds = split_values(sample_dict['DS'])
            meth_tag = split_values(sample_dict['MV'])

            extracted_data[f'{chrom}:{start}-{end}'] = (allele, meth, ds, meth_tag)
    return extracted_data

expected_data = {
'chr1:1000-1046': ([71, 71],
                    [],
                    ['TTTA-(TTTC)8-(TTTA)8-TTT', 'TTTA-(TTTC)9-(TTTA)7-TTT'],
                    []),  # TTAT
'chr10:1000-1007': ([7, 7],
                    [],
                    ['(GGA)2-A'],
                    []),  # GGA
'chr16:1000-1025': ([25, 25],
                    [],
                    ['TA-(GATG)5-GAT'],
                    []),  # TGGA
'chr2:1000-1053': ([67, 45],
                    [],
                    ['(TCTA)10-(TA)4-TC-(TA)6-TCTAT', '(TC)2-(TCTA)7-(TA)4-TCTAT'],
                    []),  # TCTA
'chr3:1000-1036': ([51, 48],
                    [0.0],
                    ['(AAGG)8-GGCGAGCAAGGAATGAAGG', '(AAGG)12'],
                    ['A']),  # CT
'chr4:1000-1057': ([72, 51],
                    [],
                    ['(CAG)24', '(CAG)17'],
                    []),  # CAG
'chr9:1000-1026': ([26, 26],
                    [],
                    [],
                    []),  # AAAAT
'chrX:1000-1060': ([93],
                    [0.01],
                    ['(CGG)10-AGG-(CGG)9-AGG-(CGG)10'],
                    ['AAAAAEAAAAAAAAAEAAAAAAAEAAEAA', 'AAAAAEAAAAAAAAAEAAAAAAAEAAEAA']),  # CGG
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
        extracted_allele, extracted_meth, extracted_ds, extracted_meth_tag = extracted_data[locus_key]
        expected_allele, expected_meth, expected_ds, expected_meth_tag = expected_data[locus_key]

        assert sorted(extracted_allele) == sorted(expected_allele), f"Allele mismatch for {locus_key}: expected {expected_allele}, got {extracted_allele}"
        assert sorted(extracted_ds) == sorted(expected_ds), f"DS mismatch for {locus_key}: expected {expected_ds}, got {extracted_ds}"
        assert sorted(extracted_meth) == sorted(expected_meth), f"Meth mismatch for {locus_key}: expected {expected_meth}, got {extracted_meth}"
        assert sorted(extracted_meth_tag) == sorted(expected_meth_tag), f"Meth_tag mismatch for {locus_key}: expected {expected_meth_tag}, got {extracted_meth_tag}"

def test_complete_pipeline(tmp_path):

 
    test_dir = Path(__file__).parent.parent / "tests/test_data"

    # Input BAM file
    input_bam = test_dir / "test_reads.bam"
    input_bed = test_dir / "test_region.bed.gz"
    input_fasta = test_dir / "test_ref.fasta"

    # Temporary output file
    output_vcf = test_dir / "tmp_test_output.vcf"
    
    # Execute the complete tool
    result = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).resolve().parent.parent / "ATARVA" / "core.py"),
            "genotype",
            "-b", str(input_bam), "-r", str(input_bed),
            "-f", str(input_fasta), "-o", str(output_vcf),
            "--decompose", "--methviz", "--karyotype", "XY"
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