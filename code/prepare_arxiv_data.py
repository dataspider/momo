import argparse
import json
import os
import re
from datetime import datetime
from multiprocessing import Pool

import networkx as nx
import pandas as pd


def get_metadata(file):
    with open(file, "r") as f:
        for line in f:
            yield line


def get_parsed_papers(file):
    metadata = get_metadata(file)
    parsed_papers = []
    for idx, paper in enumerate(metadata):
        if idx % 1000 == 0:
            print(idx, end="\r")
        parsed = json.loads(paper)
        parsed_papers.append(parsed)
    return pd.DataFrame(parsed_papers)


def create_category_files(source_folder):
    df = get_parsed_papers(f"{source_folder}/arxiv-metadata-oai-snapshot.json")
    df.categories = df.categories.map(lambda x: x.split(" "))
    df["version_1"] = df.versions.map(
        lambda x: min(
            [
                datetime.strptime(y["created"], "%a, %d %b %Y %X %Z").strftime(
                    "%Y-%m-%d"
                )
                for y in x
            ]
        )
    )
    first_version_before = "2020-11-01"
    df = df.query("version_1 < @first_version_before")
    edf = df.explode("categories")
    edf = edf.reset_index(drop=True)
    grouped = edf.groupby("categories")
    category_folder = f"{source_folder}/category-files"
    if not os.path.exists(category_folder):
        os.makedirs(category_folder)
    for category in sorted(edf.categories.unique()):
        print("Category", category, end="\r")
        grouped.get_group(category).to_pickle(
            f"{category_folder}/arxiv_{category}_{first_version_before}.gpickle.gz"
        )


def write_graphs_for_category_file(
    category_file,
    source_path,
    author_graph_path,
    bipartite_graph_path,
):
    category = category_file.split("_")[1]
    print(category)
    cdf = pd.read_pickle(f"{source_path}/{category_file}")
    ecdf = cdf.explode("authors_parsed")
    ecdf["author_string"] = ecdf.authors_parsed.map(
        lambda x: re.sub("\W", "", "_".join(x[:2]).replace(" ", "_").lower())
    )
    author_paper_G = nx.Graph()
    author_paper_G.add_nodes_from(list(ecdf.author_string.unique()), node_type="author")
    author_paper_G.add_nodes_from(list(ecdf.id.unique()), node_type="paper")
    author_paper_G.add_edges_from(list(zip(ecdf.id, ecdf.author_string)))
    author_G = nx.bipartite.projection.weighted_projected_graph(
        author_paper_G,
        nodes=[n for n, d in author_paper_G.nodes(data="node_type") if d == "author"],
    )
    nx.write_edgelist(
        author_paper_G, f"{bipartite_graph_path}/{category}.txt", data=False
    )
    nx.write_edgelist(author_G, f"{author_graph_path}/{category}.txt", data=False)


def make_temporal_graphs(category_file, source_path, year_graph_path):
    category = category_file.split("_")[1]
    for y in range(2011, 2021):
        print(category, y, end="\r")
        year = f"{y}-11-01"
        cdf = pd.read_pickle(f"{source_path}/{category_file}")
        cdf = cdf.query("version_1 < @year")
        ecdf = cdf.explode("authors_parsed")
        ecdf["author_string"] = ecdf.authors_parsed.map(
            lambda x: re.sub("\W", "", "_".join(x[:2]).replace(" ", "_").lower())
        )
        author_paper_G = nx.Graph()
        author_paper_G.add_nodes_from(
            list(ecdf.author_string.unique()), node_type="author"
        )
        author_paper_G.add_nodes_from(list(ecdf.id.unique()), node_type="paper")
        author_paper_G.add_edges_from(list(zip(ecdf.id, ecdf.author_string)))
        author_G = nx.bipartite.projection.weighted_projected_graph(
            author_paper_G,
            nodes=[
                n for n, d in author_paper_G.nodes(data="node_type") if d == "author"
            ],
        )
        nx.write_edgelist(
            author_G, f"{year_graph_path}/{category}_{year}.txt", data=False
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("source_folder", type=str)
    parser.add_argument("target_folder", type=str)
    args = parser.parse_args()
    source_folder = args.source_folder
    target_folder = args.target_folder
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)
    create_category_files(source_folder)

    source_path = f"{source_folder}/category-files"
    author_graph_path = f"{source_folder}/arxiv-author"
    bipartite_graph_path = f"{source_folder}/bipartite-graphs"
    category_files = sorted([f for f in os.listdir(source_path) if f.endswith("gz")])

    for p in [author_graph_path, bipartite_graph_path]:
        if not os.path.exists(p):
            os.makedirs(p)
    with Pool(4) as p:
        p.starmap(
            write_graphs_for_category_file,
            [
                (c, source_path, author_graph_path, bipartite_graph_path)
                for c in category_files
            ],
        )

    cs_files = [
        "arxiv_cs.LG_2020-11-01.gpickle.gz",
        "arxiv_cs.SI_2020-11-01.gpickle.gz",
    ]
    for category_file in cs_files:
        make_temporal_graphs(category_file, source_path, target_folder)
