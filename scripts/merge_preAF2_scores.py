#!/usr/bin/env python3
"""
Merge scfvtools scores and sort designs.

Outputs: merged_scores.csv
"""

import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scfv", required=True, help="scfvtools CSV file")
    parser.add_argument("--out", required=True, help="Merged output CSV")
    parser.add_argument(
        "--sort_by",
        default="scfvtools_blosum_diff_score",
        help="Column to sort by (default: scfvtools_blosum_diff_score)"
    )
    args = parser.parse_args()

    df = pd.read_csv(args.scfv)

    if args.sort_by in df.columns:
        df = df.sort_values(by=args.sort_by, ascending=False, na_position="last")
    else:
        print(f"[WARNING] sort column not found: {args.sort_by}")

    df.to_csv(args.out, index=False)
    print(f"[merge] Wrote merged scores → {args.out}")

if __name__ == "__main__":
    main()
