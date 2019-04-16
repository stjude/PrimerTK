#!/usr/bin/env python3
""" Main entry point for primer_tk """

import sys
import argparse

from primer_tk import genome_iterator, primer_cross_hyb, analyze_pcr_output, \
    __version__

from primer_tk import core

def get_args(argv):
    """ Get arguments

    Args:
        argv (list): List of strings as argv.

    Returns:
        args (Namespace): Argparse object.
    """

    parser = argparse.ArgumentParser(prog="primer_tk")
    parser.add_argument("-v", "--version", action="version",
                        version=__version__)

    subparser = parser.add_subparsers(help="Actions")
    genome_iterator.add_iterator_subparser(subparser)
    primer_cross_hyb.add_pre_subparser(subparser)
    analyze_pcr_output.add_post_subparser(subparser)

    args = parser.parse_args(argv)

    return args

def main():
    """ Main rountine """

    args = get_args(sys.argv[1:])

    action = sys.argv[1]

    if action == "iterator":
        genome_iterator.main(args)
    elif action == "pre":
        core.pre(args)
    elif action == "post":
        core.post(args)

    print(action)
    print(args)

main()
