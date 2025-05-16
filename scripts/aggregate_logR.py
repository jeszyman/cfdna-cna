#!/usr/bin/env python3

import os
import glob
import pandas as pd
from functools import reduce
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate logR values from ichorCNA .cna.seg files.")
    parser.add_argument("--seg_dir", required=True, help="Path to parent directory containing per-sample .cna.seg files.")
    parser.add_argument("--output_tsv", required=True, help="Path to write aggregated TSV output.")
    return parser.parse_args()


def main():
    args = parse_args()

    seg_files = glob.glob(os.path.join(args.seg_dir, "*/", "*.cna.seg"))
    if not seg_files:
        raise FileNotFoundError(f"No .cna.seg files found in {args.seg_dir}")

    dfs = []
    for f in seg_files:
        sample_id = os.path.basename(f).split(".")[0]
        try:
            df = pd.read_csv(f, sep="\t", usecols=["chr", "start", "end", f"{sample_id}.logR"])
        except ValueError as e:
            raise ValueError(f"Expected column '{sample_id}.logR' not found in {f}") from e
        df = df.rename(columns={f"{sample_id}.logR": sample_id})
        dfs.append(df)

    merged = reduce(lambda left, right: pd.merge(left, right, on=["chr", "start", "end"], how="outer"), dfs)
    merged = merged.sort_values(by=["chr", "start", "end"])
    merged.to_csv(args.output_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
