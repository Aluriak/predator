"""entry point for predator.

"""

import argparse


def existant_file(filepath:str) -> str:
    """Argparse type, raising an error if given file does not exists"""
    if not os.path.exists(filepath):
        raise argparse.ArgumentTypeError("file {} doesn't exists".format(filepath))
    return filepath


def cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    # positionnals
    parser.add_argument('infile', type=existant_file,
                        help="SBML or ASP file containing the graph data")
    # arguments
    parser.add_argument('--targets-file', type=str, default=None,
                        help="file containing one target per line")
    parser.add_argument('--targets', nargs='*', default=[],
                        help="targets to activate in the graph")
    # flags
    parser.add_argument('--union', action='store_true',
                        help="Print the union of all solutions")
    parser.add_argument('--intersection', action='store_true',
                        help="Print the union of all solutions")

    return parser


def parse_args(args: iter = None) -> dict:
    return cli_parser().parse_args(args)


if __name__ == '__main__':
    cliargs = parse_args()
