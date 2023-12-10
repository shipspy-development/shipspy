import argparse


def get_args(configure_parser):
    parser = argparse.ArgumentParser()

    configure_parser(parser)

    parsed_args = vars(parser.parse_args())

    return parsed_args
