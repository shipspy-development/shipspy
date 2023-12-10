import logging

from .sections import configure_sections_parser
from .dship import configure_dship_parser
from .rename import configure_rename_parser
from ._version import __version__


def get_parser():
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-v",
        "--verbose",
        metavar="DEBUG",
        help="Set the level of verbosity [DEBUG, INFO, WARNING, ERROR]",
        required=False,
        default="INFO",
    )
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=__version__,
    )

    parser.set_defaults(func=None)
    subparsers = parser.add_subparsers(dest="subcommand")
    configure_sections_parser(subparsers.add_parser("sections"))
    configure_dship_parser(subparsers.add_parser("dship"))
    configure_rename_parser(subparsers.add_parser("rename"))

    return parser


def main():
    args = get_parser().parse_args()

    logging.basicConfig(level=logging.getLevelName(args.verbose))

    return args.func(args)


if __name__ == "__main__":
    exit(main())
