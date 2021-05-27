import random

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from custom_utils import get_node_color, load_json


def get_node_colors(x, clique_nodes, star_nodes, biclique_nodes, starclique_nodes):
    if x in clique_nodes and x not in star_nodes + biclique_nodes + starclique_nodes:
        return get_node_color("clique")
    elif x in starclique_nodes and x not in clique_nodes + star_nodes + biclique_nodes:
        return get_node_color("starclique")
    elif x in star_nodes and x not in clique_nodes + starclique_nodes + biclique_nodes:
        return get_node_color("star")
    elif x in biclique_nodes and x not in clique_nodes + starclique_nodes + star_nodes:
        return get_node_color("biclique")
    elif x not in clique_nodes + starclique_nodes + star_nodes + biclique_nodes:
        return "k"
    else:
        return "silver"


path = "../data/planted"
for fn in ["dummy_1-1-1-1", "dummy_1-1-0-1"]:

    figsize = (12, 9) if fn == "dummy_1-1-1-1" else 0.8 * np.array((12, 9))
    G = nx.read_edgelist(f"{path}/{fn}.txt", nodetype=int)
    res = load_json(f"../results_mapped/planted/{fn}.json")

    clique_nodes = [
        int(item)
        for s in filter(
            lambda s: s["structure_type"] == "clique", res["macro_structures"]
        )
        for item in s["nodes"]
    ]
    star_nodes = [
        int(item)
        for s in filter(
            lambda s: s["structure_type"] == "star", res["macro_structures"]
        )
        for item in [s["hub"]] + s["spokes"]
    ]
    biclique_nodes = [
        int(item)
        for s in filter(
            lambda s: s["structure_type"] == "biclique", res["macro_structures"]
        )
        for item in s["left_nodes"] + s["right_nodes"]
    ]
    starclique_nodes = [
        int(item)
        for s in filter(
            lambda s: s["structure_type"] == "starclique", res["macro_structures"]
        )
        for item in s["left_nodes"] + s["right_nodes"]
    ]

    random.seed(1234)
    np.random.seed(1234)
    pos = nx.kamada_kawai_layout(G, scale=3)
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    nx.draw_networkx(
        G,
        pos=pos,
        nodelist=G.nodes(),
        node_color=list(
            map(
                lambda x: get_node_colors(
                    x, clique_nodes, star_nodes, biclique_nodes, starclique_nodes
                ),
                G.nodes(),
            )
        ),
        with_labels=False,
    )
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(f"../graphics/{fn}.pdf", transparent=True)
