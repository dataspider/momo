import argparse
import os

import networkx as nx
import pandas as pd


def get_tissue_graphs(filepath):
    df = pd.read_csv(filepath, sep="\t")
    df.columns = pd.Index(["protein1", "protein2", "tissue"], dtype="object")

    tissues = set(df.tissue.values)

    Gs = []
    for tissue in tissues:
        tdf = df.query("tissue == @tissue")
        edges = list(zip(tdf.protein1, tdf.protein2))
        G = nx.Graph()
        G.add_edges_from(edges)
        Gs.append((tissue, G))
    return Gs


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path_to_combined_edgelist", type=str)
    parser.add_argument("output_path", type=str)
    args = parser.parse_args()
    output_path = args.output_path
    filepath = args.path_to_combined_edgelist

    Gs = get_tissue_graphs(filepath)

    if not os.path.exists(output_path):
        os.makedirs(output_path)
    for tissue, G in Gs:
        nx.write_edgelist(G, f"{output_path}/tissue-{tissue}.txt", data=False)
