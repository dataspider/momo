import argparse
import os

import networkx as nx


def get_reference_graph(file):
    G = nx.read_gpickle(file)
    G.remove_edges_from(
        [
            (u, v, k)
            for u, v, k, d in G.edges(data="edge_type", keys=True)
            if d != "reference"
        ]
    )
    return G


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("source_folder", type=str)
    parser.add_argument("target_folder", type=str)
    args = parser.parse_args()
    source_folder = args.source_folder
    target_folder = args.target_folder
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)

    source_files = sorted(
        [f for f in os.listdir(source_folder) if f.endswith("gpickle.gz")]
    )
    for file in source_files:
        prefix = "us" if "-12-31" not in file else "de"
        G = get_reference_graph(f"{source_folder}/{file}")
        nx.write_edgelist(
            G,
            f"{target_folder}/{prefix}_{file.replace('.gpickle.gz','').replace('-12-31','')}.txt",
            data=False,
        )
