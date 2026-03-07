#!/usr/bin/env python
import src.args
from src.config import RuntimeConfig
from src.pipeline import run_prediction


def main():
    args = src.args.parse_arguments()
    config = RuntimeConfig.from_args(args)
    run_prediction(config)


if __name__ == '__main__':
    main()
