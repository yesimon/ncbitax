import argparse
import sys
from ncbitax import subset


def main():
    parser = argparse.ArgumentParser(description='Utilities for NCBI taxonomy.')
    subparsers = parser.add_subparsers(dest='subparser_name')
    subset.add_parser_subset_taxonomy(subparsers)
    args = parser.parse_args()
    if not args.subparser_name:
        parser.print_help()
        return

    args.func(args)
