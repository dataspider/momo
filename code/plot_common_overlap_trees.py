import os
from copy import deepcopy

import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from custom_utils import (
    add_tree_layout,
    get_node_color,
    get_node_marker,
    load_json,
    get_structures_added,
)


def get_shared_label(n):
    x, y = n.split("_-_")
    pos1 = int(x.split("_")[-1])
    pos2 = int(y.split("_")[-1])
    return f"{pos1} | {pos2}"


def get_msts_and_common_trees(g1, g2, alignmentpath):
    common_graph = nx.read_graphml(f"{alignmentpath}/{g1}_-_{g2}.graphml")
    common_mst = nx.maximum_spanning_tree(common_graph)
    return_values = [[None, None], [None, None]]
    for i in [0, 1]:
        g = g1 if i == 0 else g2
        mstedges_g = list(
            map(
                lambda x: (x[0].split("_-_")[i], x[1].split("_-_")[i]),
                common_mst.edges(),
            )
        )
        mstedges_g_mapped = list(
            map(lambda x: (int(x[0][-3:]) - 1, int(x[1][-3:]) - 1), mstedges_g)
        )
        g_df = pd.read_csv(
            f"../graphics/law/structure_overlap_matrix-{g}.csv", index_col=0
        ).fillna(0)
        g_graph = nx.Graph(g_df.values)
        g_mst = forest_to_tree(nx.maximum_spanning_tree(g_graph), len(g_graph))
        g_common_tree = forest_to_tree(
            nx.Graph(nx.edge_subgraph(g_graph, mstedges_g_mapped)), len(g_graph)
        )
        return_values[i][0] = deepcopy(g_mst)
        return_values[i][1] = deepcopy(g_common_tree)
    return return_values


def forest_to_tree(maxst, artificial_root):
    ccs = list(nx.connected_components(maxst))
    for c in ccs:
        component_subgraph = maxst.subgraph(c)
        component_root = max(nx.degree(component_subgraph), key=lambda tup: tup[-1])[
            0
        ]  # node with max unweighted degree
        maxst.add_edge(artificial_root, component_root, weight=np.finfo(float).eps)
    tree = nx.traversal.bfs_tree(maxst, artificial_root)
    for e in tree.edges():
        tree.edges[e]["weight"] = maxst.edges[e]["weight"]
    tree = add_tree_layout(tree, artificial_root, 10, 10)
    return tree


def get_total_edge_weight(g):
    return sum(w for u, v, w in g.edges(data="weight"))


def get_common_tree_weight_fraction(g_common_tree, g_mst):
    return get_total_edge_weight(g_common_tree) / get_total_edge_weight(g_mst)


def make_weight_fraction_df(alignmentpath):
    comparisons = sorted(x for x in os.listdir(alignmentpath) if x.endswith("graphml"))
    pairs = [tuple(x.replace(".graphml", "").split("_-_")) for x in comparisons]
    unique_ids = sorted(set(list(zip(*pairs))[0]))
    df = pd.DataFrame(index=unique_ids, columns=unique_ids, data=0.0)
    msts = {}
    trees = {}
    for g1, g2 in pairs:
        (g1_mst, g1_common_tree), (g2_mst, g2_common_tree) = get_msts_and_common_trees(
            g1, g2, alignmentpath
        )
        w1 = get_common_tree_weight_fraction(g1_common_tree, g1_mst)
        w2 = get_common_tree_weight_fraction(g2_common_tree, g2_mst)
        df.at[g1, g2] = (w1 + w2) / 2
        msts[g1] = g1_mst
        msts[g2] = g2_mst
        trees[f"{g1}_-_{g2}"] = (g1_common_tree, g2_common_tree)
    df += np.triu(df, k=1).T
    return df, msts, trees


collection = "law"
resultpath = f"../results_mapped/{collection}"
alignmentpath = f"../alignments/{collection}"
df, msts, trees = make_weight_fraction_df(alignmentpath)

pairs = [
    t.split("_-_")
    for t in trees.keys()
    if (t.split("_-_")[0].endswith("2018") or t.split("_-_")[0].endswith("2019"))
    and (t.split("_-_")[1].endswith("2018") or t.split("_-_")[1].endswith("2019"))
]
for g1, g2 in pairs:
    common_graph = nx.read_graphml(f"{alignmentpath}/{g1}_-_{g2}.graphml")
    json1 = load_json(f"{resultpath}/{g1}.json")
    json2 = load_json(f"{resultpath}/{g2}.json")
    n1 = json1["n"]
    n2 = json2["n"]
    s1 = {x["position"]: x for x in get_structures_added(json1)}
    s2 = {x["position"]: x for x in get_structures_added(json2)}
    G = common_graph
    maxst = nx.tree.maximum_spanning_tree(G)
    artificial_root = G.number_of_nodes()
    ccs = list(nx.connected_components(maxst))
    for c in ccs:
        component_subgraph = maxst.subgraph(c)
        component_root = max(nx.degree(component_subgraph), key=lambda tup: tup[-1])[
            0
        ]  # node with max unweighted degree
        maxst.add_edge(artificial_root, component_root, weight=np.finfo(float).eps)
    tree = nx.traversal.bfs_tree(maxst, artificial_root)
    for e in tree.edges():
        tree.edges[e]["weight"] = maxst.edges[e]["weight"]
    tree = add_tree_layout(tree, artificial_root, 10, 10)
    _, ax = plt.subplots(1, 1, figsize=(12, 12))
    for node in tree.nodes():
        x = tree.nodes[node]["x"]
        y = tree.nodes[node]["y"]
        if type(node) != int:
            structure1, structure2 = node.split("_-_")
            size1, size2 = (
                s1[structure1]["n_nodes_total"],
                s2[structure2]["n_nodes_total"],
            )
            idx1, idx2 = int(structure1[-3:]), int(structure2[-3:])
            structure_type = structure1.split("_")[-2]
            marker = get_node_marker(structure_type)
        else:
            structure_type = None
            marker = "X"
        color = get_node_color(structure_type) if structure_type else "k"
        for succ in tree.successors(node):
            ax.plot(
                [x, tree.nodes[succ]["x"]],
                [y, tree.nodes[succ]["y"]],
                "-k",
                linewidth=max(tree.edges[node, succ]["weight"] * 10, 1),
                zorder=1,
                alpha=1,
            )
        ax.scatter(
            x,
            y,
            color=color,
            s=(size1 + size2) / 2 * 6 if structure_type else 300,
            marker=marker,
            zorder=2,
            alpha=1,
        )
    #             if structure_type:
    #                 ax.annotate(f'{idx1}|{idx2}' if idx1 != idx2 else str(idx1),(x,y),ha='center',va='center',fontsize=10)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(
        f"../graphics/common-models/shared_{g1}_-_{g2}.pdf",
        transparent=True,
        bbox_inches="tight",
    )
    plt.close()
