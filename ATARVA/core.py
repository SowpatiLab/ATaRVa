#!/usr/bin/env python

"""
    ATaRVa (pronunced as atharva) is a tool designed to analyse tandem repeat variation
    from genome sequencing data. Besides genotyping the locus ATaRVa provides functionality
    to decompose motifs, analyse base modifications along with interactive visualizations.
"""

import sys
import argparse as ap

from ATARVA.genotype import genotype_parser
from ATARVA.merge import merge_parser
from ATARVA.version import __version__

contact = """\nFor queries or suggestions, please contact:
Divya Tej Sowpati - tej at csirccmb dot org
Abishek Kumar S - abishekks at csirccmb dot org
Akshay Kumar Avvaru - avvaruakshay at gmail dot com"""

def main():
    parser = ap.ArgumentParser(prog="atarva",
                               add_help=False,
                               formatter_class=ap.RawTextHelpFormatter)

    parser._action_groups.pop()

    print(f"ATaRVa - Analysis of Tandem Repeat Variants\nSowpati Lab\n")

    parser.add_argument('-h', '--help', action='store_true', help="Print help")
    parser.add_argument('-v', '--version', action='version', version=f'ATaRVa version {__version__}', help="Print version")

    subparsers = parser.add_subparsers(dest="command")

    genotype_parser(subparsers)
    merge_parser(subparsers)

    args = parser.parse_args()

    if args.command is None:
        print("Usage:")
        print("    atarva [OPTIONS] <COMMAND>\n")
        print("Commands:")
        for name, sp in subparsers.choices.items():
            print(f"  {name:<9} {sp.description}")
        print("\nOptions:")
        print("  -h, --help     Print help")
        print("  -v, --version  Print version")
        print(contact)
        sys.exit()

    args.func(args)

if __name__ == '__main__':
    main()
