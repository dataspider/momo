import argparse
import os

import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("collection", type=str)
    parser.add_argument(
        "--variant", "-v", choices=["overlap", "jaccard"], default="overlap"
    )
    args = parser.parse_args()
    collection = args.collection
    variant = args.variant

    alignmentpath = (
        f"../alignments/{collection}"
        if variant == "overlap"
        else f"../alignments-jaccard/{collection}"
    )
    if not os.path.exists(alignmentpath):
        raise Exception(f"Path {alignmentpath} does not exist!")

    alignment_files = sorted(
        [
            f
            for f in os.listdir(alignmentpath)
            if f.endswith(".csv") and "structure-alignments-overlap" not in f
        ]
    )
    dfs = [pd.read_csv(f"{alignmentpath}/{f}") for f in alignment_files]
    df = pd.concat(dfs, ignore_index=True)
    df.to_csv(f"{alignmentpath}/structure-alignments-{variant}.csv")
