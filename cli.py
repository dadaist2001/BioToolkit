import argparse
from BioToolkit import filter_fastq

# super miniimalistic
def parse_args():
    parser = argparse.ArgumentParser(description="FastQ filter")

    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)

    parser.add_argument("--min_len", type=int, default=0)
    parser.add_argument("--max_len", type=int, default=float("inf"))

    parser.add_argument("--min_gc", type=float, default=0)
    parser.add_argument("--max_gc", type=float, default=100)

    return parser.parse_args()


def main():
    args = parse_args()

    filter_fastq(
        args.input,
        args.output,
        (args.min_len, args.max_len),
        (args.min_gc, args.max_gc),
    )


if __name__ == "__main__":
    main()
