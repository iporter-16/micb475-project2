#!/usr/bin/env python

import argparse
from importlib.metadata import version
from picrust2.util import check_files_exist, shuffle_predictions

parser = argparse.ArgumentParser(

    description='Shuffles sequence ids of prediction table for a specified '
                'number of replicates. Note that this means that the '
                'prediction table is not totally scrambled across all '
                'colunmns: each predicted genome is still the same, but the '
                'sequence ids (e.g. the ASV id) linked with each predicted '
                'genome is randomized.',

    epilog='''

Usage example:

shuffle_predictions.py -i EC_predicted.tsv.gz -o EC_predicted_shuffled --rep 5''',

    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input', metavar='OUTPUT', required=True,
                    type=str,
                    help='Path to input directory.')

parser.add_argument('-o', '--outdir', metavar='OUTPUT', required=True,
                    type=str,
                    help='Path to output directory.')

parser.add_argument('-r', '--rep', default=1, metavar='INT', required=False,
                    type=int, help="Number of shuffled replicates to create "
                                   "(default: %(default)d).")

parser.add_argument('-s', '--seed', default=None, metavar='INT',
                    required=False, type=int,
                    help="Random seed: set this if you want reproducible "
                         "shufflings (default: None).")

parser.add_argument('-v', '--version', default=False, action='version',
                    version="PICRUSt2 " + version('PICRUSt2'))


def main():

    args = parser.parse_args()

    # Check that input file exists.
    check_files_exist([args.input])

    shuffle_predictions(input=args.input, outdir=args.outdir, rep=args.rep,
                        seed=args.seed)


if __name__ == '__main__':
    main()
