import argparse
import os
from itertools import combinations, product
from multiprocessing import Pool

import networkx as nx
import pandas as pd

from custom_utils import get_jaccard, get_structures_added, load_json
from statics import STRUCTURE_TYPES


def get_product_graph(g1, g2, graphicspath, resultspath):
    o1 = f"{graphicspath}/structure_overlap_matrix-{g1}.csv"
    o2 = f"{graphicspath}/structure_overlap_matrix-{g2}.csv"
    j1 = f"{resultspath}/{g1}.json"
    j2 = f"{resultspath}/{g2}.json"
    json1 = load_json(j1)
    json2 = load_json(j2)
    s1 = {
        idx: (s["n_nodes_total"], s["structure_type"], s)
        for idx, s in enumerate(get_structures_added(json1))
    }
    s2 = {
        idx: (s["n_nodes_total"], s["structure_type"], s)
        for idx, s in enumerate(get_structures_added(json2))
    }
    df1 = pd.read_csv(o1, index_col=0).fillna(0)
    df1.columns = df1.columns.astype(int)
    df2 = pd.read_csv(o2, index_col=0).fillna(0)
    df2.columns = df2.columns.astype(int)

    G1 = nx.Graph(df1)
    G2 = nx.Graph(df2)

    nodes = [(x, y) for x, y in product(s1.keys(), s2.keys()) if s1[x][1] == s2[y][1]]
    edges = [
        (
            tup[0],
            tup[1],
            G1.edges[tup[0][0], tup[1][0]]["weight"]
            * G2.edges[tup[0][-1], tup[1][-1]]["weight"],
        )
        for tup in filter(
            lambda tup: tup[0][0] != tup[1][0]
            and tup[0][-1] != tup[1][-1]
            and G1.edges.get([tup[0][0], tup[1][0]])
            and G2.edges.get([tup[0][-1], tup[1][-1]]),
            combinations(nodes, 2),
        )
    ]
    G12 = nx.Graph()
    G12.add_nodes_from(nodes)
    G12.add_weighted_edges_from(edges)
    return s1, s2, G1, G2, G12


def get_overlap_alignment_edges(G12):
    selected_edges = []
    while True:
        if nx.number_of_edges(G12) == 0:
            break
        new_edge = max(
            G12.edges(data="weight"),
            key=lambda tup: (
                tup[-1],
                -abs(tup[0][0] - tup[0][1]),
                -abs(tup[1][0] - tup[1][1]),
            ),
        )
        selected_edges.append(new_edge)
        nodes_to_remove = list(
            filter(
                lambda x: (x[0] == new_edge[0][0] and x[1] != new_edge[0][1])
                or (x[1] == new_edge[0][1] and x[0] != new_edge[0][0])
                or (x[0] == new_edge[1][0] and x[1] != new_edge[1][1])
                or (x[1] == new_edge[1][1] and x[0] != new_edge[1][0]),
                G12.nodes(),
            )
        )
        G12.remove_edge(*new_edge[:-1])
        G12.remove_nodes_from(nodes_to_remove)

        print(
            "Constructed G12 with",
            "n =",
            nx.number_of_nodes(G12),
            "and m =",
            nx.number_of_edges(G12),
            end="\r",
        )
    return selected_edges


def get_common_model_graph(G12):
    selected_edges = get_overlap_alignment_edges(G12)
    common_graph = nx.Graph()
    common_graph.add_weighted_edges_from(selected_edges)
    return common_graph


def complete_with_greedy(common_graph, G12):
    leftover_candidates = set(G12.nodes()) - set(common_graph.nodes())
    while leftover_candidates:
        x, y = min(leftover_candidates)
        common_graph.add_node((x, y))
        leftover_candidates -= {
            (u, v) for (u, v) in leftover_candidates if u == x or v == y
        }
    return common_graph


def get_overlap(g1, g2, structures1, structures2, structure_type):
    s1 = list(filter(lambda x: x["structure_type"] == structure_type, structures1))
    s2 = list(filter(lambda x: x["structure_type"] == structure_type, structures2))
    overlap = pd.DataFrame(
        index=list(map(lambda x: x["position"], s1)),
        columns=list(map(lambda x: x["position"], s2)),
    )
    for idx1, structure1 in enumerate(s1):
        for idx2, structure2 in enumerate(s2):
            overlap.iat[idx1, idx2] = get_jaccard(
                structure1, structure2, structure_type
            )
    return overlap


def match_structures_by_overlap(overlap):
    matched = set()
    pairs = []
    found_match = True
    while found_match:
        found_match = False
        remaining = (
            overlap.query("index not in @matched").T.query("index not in @matched").T
        )
        if min(remaining.shape) > 0:
            max_u = remaining.max().idxmax()
            other = (
                overlap.query("index not in @matched")
                .T.query("index not in @matched")
                .T[[max_u]]
                .T
            )
            if min(other.shape) > 0:
                max_v = other.max().idxmax()
                pairs.append((max_v, max_u, overlap.at[max_v, max_u]))
                found_match = True
                matched = matched | {max_u, max_v}
    return pairs


def get_matches_with_node_alignment(g1, g2, structures1, structures2):
    clique_overlap = get_overlap(g1, g2, structures1, structures2, "clique")
    clique_pairs = match_structures_by_overlap(clique_overlap)
    assert len(clique_pairs) == min(clique_overlap.shape), "Clique alignment problem"
    star_overlap = get_overlap(g1, g2, structures1, structures2, "star")
    star_pairs = match_structures_by_overlap(star_overlap)
    assert len(star_pairs) == min(star_overlap.shape), "Star alignment problem"
    biclique_overlap = get_overlap(g1, g2, structures1, structures2, "biclique")
    biclique_pairs = match_structures_by_overlap(biclique_overlap)
    assert len(biclique_pairs) == min(
        biclique_overlap.shape
    ), "Biclique alignment problem"
    starclique_overlap = get_overlap(g1, g2, structures1, structures2, "starclique")
    starclique_pairs = match_structures_by_overlap(starclique_overlap)
    assert len(starclique_pairs) == min(
        starclique_overlap.shape
    ), "Starclique alignment problem"
    matched = sorted(
        clique_pairs + star_pairs + biclique_pairs + starclique_pairs,
        key=lambda x: x[-1],
        reverse=True,
    )
    unmatched1 = (
        sorted(set(clique_overlap.index) - set(map(lambda x: x[0], matched)))
        + sorted(set(star_overlap.index) - set(map(lambda x: x[0], matched)))
        + sorted(set(biclique_overlap.index) - set(map(lambda x: x[0], matched)))
        + sorted(set(starclique_overlap.index) - set(map(lambda x: x[0], matched)))
    )
    unmatched2 = (
        sorted(set(clique_overlap.columns) - set(map(lambda x: x[1], matched)))
        + sorted(set(star_overlap.columns) - set(map(lambda x: x[1], matched)))
        + sorted(set(biclique_overlap.columns) - set(map(lambda x: x[1], matched)))
        + sorted(set(starclique_overlap.columns) - set(map(lambda x: x[1], matched)))
    )
    return matched, unmatched1, unmatched2


def get_structures_for_type(structures, structure_type):
    selected = list(filter(lambda x: x["structure_type"] == structure_type, structures))
    return list(map(lambda x: x["position"], selected))


def get_matches_greedy(structures1, structures2):
    matched = []
    unmatched1 = []
    unmatched2 = []
    for structure_type in STRUCTURE_TYPES:
        s1 = get_structures_for_type(structures1, structure_type)
        s2 = get_structures_for_type(structures2, structure_type)
        matched.extend(list(zip(s1, s2)))
        if len(s1) > len(s2):
            unmatched1.extend(s1[len(s2) :])
        elif len(s1) < len(s2):
            unmatched2.extend(s2[len(s1) :])
    return matched, unmatched1, unmatched2


def get_graph_names(path):
    return sorted(
        [
            x.split(".json")[0]
            for x in [f for f in os.listdir(path) if f.endswith("json")]
        ]
    )


def save_jaccard_alignments(g1, idx, columns, graph_names, mapped_path, alignmentpath):
    print(
        f"Starting alignments for graph {idx + 1}/{len(graph_names)} ({g1})...",
        end="\r",
    )
    large_df = pd.DataFrame(columns=columns)
    for g2 in graph_names[idx:]:
        alignment = f"{g1}_-_{g2}"
        json1 = load_json(f"{mapped_path}/{g1}.json")
        json2 = load_json(f"{mapped_path}/{g2}.json")
        structures1 = get_structures_added(json1)
        structures2 = get_structures_added(json2)

        matched, unmatched1, unmatched2 = get_matches_with_node_alignment(
            g1, g2, structures1, structures2
        )
        df = pd.concat(
            [
                pd.DataFrame(
                    matched, columns=["g1_structure", "g2_structure", "jaccard"]
                ),
                pd.DataFrame([unmatched1], index=["g1_structure"]).T,
                pd.DataFrame([unmatched2], index=["g2_structure"]).T,
            ],
            sort=True,
        )[["g1_structure", "g2_structure", "jaccard"]]
        df["alignment"] = alignment
        large_df = large_df.append(df, ignore_index=True, sort=True)
    large_df[columns].sort_values(
        ["alignment", "jaccard", "g1_structure", "g2_structure"],
        ascending=[True, False, True, True],
    ).to_csv(f"{alignmentpath}/{g1}-jaccard.csv", index=False)


def save_overlap_alignments(
    g1, idx, columns, graph_names, mapped_path, alignmentpath, graphicspath
):
    print(
        f"Starting alignments for graph {idx + 1}/{len(graph_names)} ({g1})...",
        end="\r",
    )
    large_df = pd.DataFrame(columns=columns)
    for g2 in graph_names[idx:]:
        alignment = f"{g1}_-_{g2}"
        print(f"Aligning {g1} and {g2}...", end="\r")
        if graphicspath is not None:
            s1, s2, G1, G2, G12 = get_product_graph(g1, g2, graphicspath, mapped_path)
            common_graph = get_common_model_graph(G12)
            common_graph = complete_with_greedy(common_graph, G12)
            unmatched_s1 = [
                (s1[x][-1]["position"], None, alignment)
                for x in sorted(s1.keys())
                if x not in map(lambda tup: tup[0], common_graph.nodes())
            ]
            unmatched_s2 = [
                (None, s2[x][-1]["position"], alignment)
                for x in sorted(s2.keys())
                if x not in map(lambda tup: tup[1], common_graph.nodes())
            ]
            common_graph = nx.relabel_nodes(
                common_graph,
                {
                    n: (s1[n[0]][-1]["position"], s2[n[1]][-1]["position"])
                    for n in common_graph.nodes()
                },
            )
            matched_nodes = [(*x, alignment) for x in list(common_graph.nodes())]
            ndf = pd.DataFrame(
                list(matched_nodes + unmatched_s1 + unmatched_s2),
                columns=["g1_structure", "g2_structure", "alignment"],
            )
            common_graph = nx.relabel_nodes(
                common_graph, {n: f"{n[0]}_-_{n[1]}" for n in common_graph.nodes()}
            )
            nx.write_graphml(common_graph, f"{alignmentpath}/{alignment}.graphml")
        else:  # just greedy - currently assumes synthetic models, ie, no number of edges and loop-freeness constraints at the beginning of macro_structures
            structures1 = load_json(f"{mapped_path}/{g1}.json")["macro_structures"]
            structures2 = load_json(f"{mapped_path}/{g2}.json")["macro_structures"]
            matched, unmatched1, unmatched2 = get_matches_greedy(
                structures1, structures2
            )
            matched_nodes = [(s1, s2, alignment) for (s1, s2) in matched]
            unmatched_s1 = [(s1, None, alignment) for s1 in unmatched1]
            unmatched_s2 = [(None, s2, alignment) for s2 in unmatched2]
            ndf = pd.DataFrame(
                list(matched_nodes + unmatched_s1 + unmatched_s2),
                columns=["g1_structure", "g2_structure", "alignment"],
            )
        if len(ndf) > 0:
            large_df = large_df.append(ndf, ignore_index=True, sort=True)
    large_df[columns].sort_values(
        ["alignment", "g1_structure", "g2_structure"], ascending=True
    ).to_csv(f"{alignmentpath}/{g1}-overlap.csv", index=False)


if __name__ == "__main__":
    # Parser setup
    parser = argparse.ArgumentParser()
    parser.add_argument("collection")
    parser.add_argument("--node_alignment", "-a", action="store_true")
    parser.add_argument("--parallelize", "-p", type=int, default=0)
    parser.add_argument("--no_overlap", "-o", action="store_true")
    args = parser.parse_args()
    collection = args.collection
    node_alignment_present = args.node_alignment
    parallelize = args.parallelize
    no_overlap = args.no_overlap

    # Path and file name setup
    mapped_path = f"../results_mapped/{collection}"
    if not os.path.exists(mapped_path):
        raise Exception(f"Path {mapped_path} does not exist!")
    graph_names = get_graph_names(mapped_path)

    if node_alignment_present:
        columns = ["g1_structure", "g2_structure", "jaccard", "alignment"]
        alignmentpath = f"../alignments-jaccard/{collection}"
    else:
        columns = ["g1_structure", "g2_structure", "alignment"]
        if no_overlap:
            graphicspath = None
        else:
            graphicspath = f"../graphics/{collection}"
            if not os.path.exists(graphicspath):
                raise Exception(f"Path {graphicspath} does not exist!")
        alignmentpath = f"../alignments/{collection}"
    if not os.path.exists(alignmentpath):
        os.makedirs(alignmentpath)

    # Actual computation
    if parallelize and node_alignment_present:
        with Pool(parallelize) as p:
            p.starmap(
                save_jaccard_alignments,
                [
                    (g1, idx, columns, graph_names, mapped_path, alignmentpath)
                    for idx, g1 in enumerate(graph_names)
                ],
            )
    elif parallelize and not node_alignment_present:
        with Pool(parallelize) as p:
            p.starmap(
                save_overlap_alignments,
                [
                    (
                        g1,
                        idx,
                        columns,
                        graph_names,
                        mapped_path,
                        alignmentpath,
                        graphicspath,
                    )
                    for idx, g1 in enumerate(graph_names)
                ],
            )
    elif node_alignment_present:
        for idx, g1 in enumerate(graph_names):
            save_jaccard_alignments(
                g1, idx, columns, graph_names, mapped_path, alignmentpath
            )
    else:
        for idx, g1 in enumerate(graph_names):
            save_overlap_alignments(
                g1, idx, columns, graph_names, mapped_path, alignmentpath, graphicspath
            )
