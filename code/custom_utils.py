import functools
import json
import re
from collections import Counter

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns

from statics import STRUCTURE_TYPES

sns.set_style("whitegrid")
plt.rcParams["figure.figsize"] = (18, 12)
plt.rcParams["font.size"] = 12
np.random.seed(1234)


def get_jaccard(structure1, structure2, structure_type):
    if structure_type == "clique":
        nodes1 = set(structure1["nodes"])
        nodes2 = set(structure2["nodes"])
        overlap = nodes1.intersection(nodes2)
        union = nodes1.union(nodes2)
        return len(overlap) / len(union)
    if structure_type in ["biclique", "starclique"]:
        left1, left2 = set(structure1["left_nodes"]), set(structure2["left_nodes"])
        right1, right2 = set(structure1["right_nodes"]), set(structure2["right_nodes"])
        left_overlap = left1.intersection(left2)
        left_union = left1.union(left2)
        right_overlap = right1.intersection(right2)
        right_union = right1.union(right2)
        return (
            len(left_overlap) / len(left_union) + len(right_overlap) / len(right_union)
        ) / 2
    if structure_type == "star":
        hub1, hub2 = {structure1["hub"]}, {structure2["hub"]}
        spokes1, spokes2 = set(structure1["spokes"]), set(structure2["spokes"])
        hub_overlap = hub1.intersection(hub2)
        hub_union = hub1.union(hub2)
        spoke_overlap = spokes1.intersection(spokes2)
        spoke_union = spokes1.union(spokes2)
        return (
            len(hub_overlap) / len(hub_union) + len(spoke_overlap) / len(spoke_union)
        ) / 2
    raise Exception(f"Unknown structure type: {structure_type}!")


def get_dataset_color(dataset):
    if dataset.startswith("ors") or dataset.startswith("asb"):
        return "dodgerblue"
    elif dataset.startswith("orp") or dataset.startswith("asp"):
        return "lightskyblue"
    elif dataset.startswith("usl") or dataset.startswith("lus"):
        return "r"
    elif dataset.startswith("del") or dataset.startswith("lde"):
        return "darkorange"
    elif dataset.startswith("clg"):
        return "purple"
    elif dataset.startswith("csi"):
        return "magenta"
    elif "bio$_{\mathcal{A}}" in dataset:
        return "green"
    elif dataset.startswith("bio\n") or dataset.startswith("bio"):
        return "g"
    elif dataset.startswith("bag") or dataset.startswith("rba"):
        return "gray"
    elif dataset.startswith("erg") or dataset.startswith("rer"):
        return "darkgray"
    else:
        raise Exception(dataset)


def load_json(file):
    """
    load a json file as a dictionary
    """
    with open(file) as f:
        model_json = json.load(f)
    return model_json


def load_log(file):
    """
    load a log file as a list of log file lines
    """
    with open(file) as f:
        model_log = f.read().split("\n")
    return model_log


def create_df(model_json):
    """
    convert the model json computed by julia into a pd.DataFrame
    """
    tuples = list(
        zip(
            model_json["macro_structures"],
            model_json["macro_structure_description_lengths"],
            model_json["description_lengths_over_time"],
        )
    )
    df = pd.DataFrame(
        tuples, columns=["structure", "structure_cost", "description_length"]
    )
    df["n_edges_total"] = [
        x.get("n_edges_total", model_json["m"]) for x in df.structure
    ]
    df["n_nodes_total"] = [
        x.get("n_nodes_total", model_json["n"]) for x in df.structure
    ]
    df["structure_type"] = [x.get("structure_type") for x in df.structure]
    df["structure_shape"] = [
        get_node_marker(x) if x in STRUCTURE_TYPES else "X" for x in df.structure_type
    ]
    df["structure_color"] = [
        get_node_color(x) if x in STRUCTURE_TYPES else "k" for x in df.structure_type
    ]
    return df


def create_progression_plot(df, save_path=None):
    """
    position of structure in the sequence on x, description length after adding structure on y, color signaling structure type, size signalling number of edges
    """
    scattertuples = list(
        zip(
            df.index - 1,
            df.description_length / df.description_length.max(),
            df.n_edges_total,
            df.structure_color,
            df.structure_shape,
        )
    )
    for t in reversed(scattertuples[1:]):
        plt.scatter(t[0], t[1], s=t[2] if t[3] != "k" else 10, c=t[3], marker="o")
    plt.xticks(range(0, len(scattertuples[1:]) + 1, 2))
    plt.xlim(-1, len(scattertuples[1:]) + 1)
    plt.xlabel("Selected structure")
    plt.ylabel("Total description length after structure selected")
    plt.title(save_path)
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path)
        plt.close()


def create_size_plot(model_json, x_granularity, y_granularity, save_path=None):
    """
    number of nodes on x, number of edges on y, color signaling structure type
    """
    structure_types, n_nodes, n_edges = list(
        zip(
            *(
                [
                    (s["structure_type"], s.get("n_nodes_total", 0), s["n_edges_total"])
                    for s in model_json["macro_structures"]
                ]
            )
        )
    )
    plt.scatter(
        n_nodes[2:],
        n_edges[2:],
        c=list(map(get_node_color, structure_types[2:])),
    )
    plt.xlabel("Number of Nodes")
    plt.xticks(range(0, max(n_nodes[2:]) + x_granularity, x_granularity))
    plt.yticks(range(0, max(n_edges[2:]) + y_granularity, y_granularity))
    plt.ylim(0, max(n_edges[2:]) + y_granularity)
    plt.ylabel("Number of Edges")
    plt.title(save_path)
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path)
        plt.close()


def get_structures_added(model_json):
    """
    return list of dicts, with each dict a structure added in the model building process (i.e., generic structures are excluded)
    """
    return model_json["macro_structures"][2:]


def get_node_sets(structures_added):
    """
    return a list of lists, with each inner list holding the nodes of a structure
    """
    return [_get_nodes(structure) for structure in structures_added]


def _get_nodes(structure):
    """
    helper for get_node_sets
    """
    if structure["structure_type"] in ["biclique", "starclique"]:
        return structure["left_nodes"] + structure["right_nodes"]
    elif structure["structure_type"] == "clique":
        return structure["nodes"]
    elif structure["structure_type"] == "star":
        return [structure["hub"]] + structure["spokes"]
    else:
        raise Exception(f"Unknown structure type {structure['structure_type']}!")


def get_structure_dfs(structures_added, node_sets):
    """
    return two pd.DataFrame objects encoding the node overlap between structures: abs_df (# nodes in the overlap), rel_df (jaccard similarity)
    """
    abs_df = pd.DataFrame(
        index=range(len(structures_added)),
        columns=range(len(structures_added)),
        data=np.nan,
    )
    rel_df = pd.DataFrame(
        index=range(len(structures_added)),
        columns=range(len(structures_added)),
        data=np.nan,
    )
    for idx in range(0, len(node_sets) - 1):
        for idx2 in range(idx + 1, len(node_sets)):
            abs_df.at[idx, idx2] = len(
                set(node_sets[idx]).intersection(set(node_sets[idx2]))
            )
            abs_df.at[idx2, idx] = abs_df.at[idx, idx2]
            rel_df.at[idx, idx2] = len(
                set(node_sets[idx]).intersection(set(node_sets[idx2]))
            ) / len(set(node_sets[idx]).union(set(node_sets[idx2])))
            rel_df.at[idx2, idx] = rel_df.at[idx, idx2]
    return abs_df, rel_df


def _get_n_nodes_covered(node_sets):
    """
    helper for get_fraction_nodes_covered
    """
    return len(set(functools.reduce(lambda x, y: x + y, node_sets, [])))


def get_fraction_nodes_covered(node_sets, model_json):
    return _get_n_nodes_covered(node_sets) / model_json["n"]


def plot_overlap_heatmap(df, save_path=None):
    """
    structures added to model on x and y, similarity as per df as color, default colormap, robust=False
    """
    sns.heatmap(df, square=True)
    if save_path is not None:
        plt.savefig(save_path)
        plt.close()


def create_rooted_bfs_tree(df, layout=False):
    G = nx.Graph(df.fillna(0))
    maxst = nx.tree.maximum_spanning_tree(G)
    artificial_root = G.number_of_nodes()
    ccs = list(nx.connected_components(G))
    for c in ccs:
        component_subgraph = maxst.subgraph(c)
        component_root = max(nx.degree(component_subgraph), key=lambda tup: tup[-1])[
            0
        ]  # node with max unweighted degree
        maxst.add_edge(artificial_root, component_root, weight=np.finfo(float).eps)
    tree = nx.traversal.bfs_tree(maxst, artificial_root)
    for e in tree.edges():
        tree.edges[e]["weight"] = maxst.edges[e]["weight"]
    if layout:
        pos = nx.layout.kamada_kawai_layout(maxst, weight=None)
        return tree, pos
    else:
        return tree


def add_tree_layout(G, root, node_sep, level_sep):
    for node in G.nodes():
        G.nodes[node]["y"] = -level_sep * nx.dijkstra_path_length(
            G, root, node, weight=None
        )
    base = 0
    for node in nx.dfs_postorder_nodes(G, root):
        succ = sorted(list(G.successors(node)), reverse=True)
        if len(succ) < 1:
            G.nodes[node]["x"] = base + node_sep
            base += node_sep
        else:
            xmin = min([G.nodes[node]["x"] for node in succ])
            xmax = max([G.nodes[node]["x"] for node in succ])
            G.nodes[node]["x"] = xmin + (xmax - xmin) / 2
    for node in G.nodes:
        G.nodes[node]["x"] = -G.nodes[node]["x"]
    return G


def add_color(G, df):
    for node in G.nodes():
        G.nodes[node]["color"] = (
            df.at[node + 2, "structure_color"] if node != len(df) - 2 else "k"
        )
    return G


def plot_tree(G, df, save_path=None):
    G = add_color(G, df)
    _, ax = plt.subplots(1, 1, figsize=(12, 12))
    for node in G.nodes():
        x = G.nodes[node]["x"]
        y = G.nodes[node]["y"]
        color = G.nodes[node]["color"]
        for succ in G.successors(node):
            ax.plot(
                [x, G.nodes[succ]["x"]],
                [y, G.nodes[succ]["y"]],
                "-k",
                linewidth=max(G.edges[node, succ]["weight"] * 10, 1),
                zorder=1,
                alpha=1,
            )
        ax.scatter(
            x,
            y,
            color=color,
            s=df.at[node + 2, "n_nodes_total"] * 6 if node != len(df) - 2 else 300,
            marker=df.at[node + 2, "structure_shape"] if node != len(df) - 2 else "X",
            zorder=2,
            alpha=1,
        )
        # if node != len(df) - 2:
        #     ax.annotate(node + 1, (x, y), fontsize=10, ha="center", va="center")
    plt.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.axis("off")
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, transparent=True, bbox_inches="tight")
        plt.close()


def plot_structure_tree(tree, layout, df, save_path=None):
    """
    plot structure tree in basic kamada kawai layout;
    structure identifiers in order of structure addition and color corresponding to structure type (artificial root node black)
    """
    nx.draw_networkx_edges(tree, pos=layout)
    for node, (x, y) in layout.items():
        plt.scatter(
            x,
            y,
            color=df.at[node + 2, "structure_color"] if node != len(df) - 2 else "k",
            s=df.at[node + 2, "n_nodes_total"] * 6 if node != len(df) - 2 else 100,
            marker=df.at[node + 2, "structure_shape"] if node != len(df) - 2 else "X",
            zorder=2,
            alpha=0.8,
        )
    labels = {idx: idx + 1 for idx in tree.nodes()}
    nx.draw_networkx_labels(tree, pos=layout, labels=labels)
    plt.axis("off")
    if save_path is not None:
        plt.savefig(save_path)
        plt.close()


def write_plots_for_model_json(
    json_path,
    save_base,
    x_granularity_size,
    y_granularity_size,
):
    """
    end-to-end plot generation for json file at given json_path
    """
    print(f"Starting {json_path}...")
    model_json = load_json(json_path)
    save_base = save_base.split("_size")[0]
    df = create_df(model_json)
    df.to_csv(re.sub("figure", "structures", save_base) + ".csv", index=False)
    structures_added = get_structures_added(model_json)
    node_sets = get_node_sets(structures_added)
    try:
        abs_df, rel_df = get_structure_dfs(structures_added, node_sets)
        rel_df.to_csv(re.sub("figure", "structure_overlap_matrix", save_base) + ".csv")
        tree, layout = create_rooted_bfs_tree(rel_df, layout=True)
        plot_tree(
            add_tree_layout(tree, tree.number_of_nodes() - 1, 10, 10),
            df,
            re.sub("figure", "tree-hierarchical", save_base) + ".pdf",
        )
        plot_structure_tree(
            tree, layout, df, re.sub("figure", "tree-kamada", save_base) + ".pdf"
        )
        G = create_overlap_quotient_graph(structures_added, abs_df, model_json["n"])
        plot_overlap_quotient_graph(
            G,
            df,
            model_json["n"],
            re.sub("figure", "overlap-quotient", save_base) + ".pdf",
        )
        G = create_structure_quotient_graph(node_sets, save_base)
        plot_structure_quotient_graph(
            G,
            node_sets,
            structures_added,
            save_path=re.sub("figure", "structure-quotient", save_base) + ".pdf",
        )
    except:
        print(
            f"Error for overlap dataframes or graph plots: {json_path} - moving on..."
        )
    try:
        create_progression_plot(
            df,
            re.sub("figure", "progress", save_base) + ".pdf",
        )
    except:
        print(f"Error for progression plot: {json_path} - moving on...")
    try:
        create_size_plot(
            model_json,
            x_granularity_size,
            y_granularity_size,
            re.sub("figure", "sizes", save_base) + ".pdf",
        )
    except:
        print(f"Error for size plot: {json_path} - moving on...")


def get_edgelist_separator(edgelist_path):
    with open(edgelist_path) as f:
        for line in f:
            if not line.startswith("#"):
                if "\t" in line:
                    return "\t"
                elif "," in line:
                    return ","
                elif " " in line:
                    return " "
                else:
                    raise


def create_structure_quotient_graph(nodes, save_base):
    nodemap_path = (
        re.sub("figure-", "", re.sub("graphics/", "results/", save_base))
        + "-nodemap.csv"
    )
    nodemap = pd.read_csv(nodemap_path)
    edgelist_path = (
        re.sub("figure-", "", re.sub("graphics/", "data/", save_base)) + ".txt"
    )
    edges = pd.read_csv(
        edgelist_path,
        sep=get_edgelist_separator(edgelist_path),
        comment="#",
        header=None,
        usecols=[0, 1],
    ).rename({0: "u", 1: "v"}, axis=1)
    new_edges = edges.merge(nodemap, left_on="u", right_on="original_id").merge(
        nodemap, left_on="v", right_on="original_id", suffixes=("_u", "_v")
    )[["julia_id_u", "julia_id_v"]]
    assert len(edges) == len(new_edges)
    nodes_to_structures = get_nodes_to_structures(nodes)
    G = nx.MultiGraph()
    G.add_nodes_from(range(1, len(nodes) + 1))
    for u, v in zip(new_edges.julia_id_u, new_edges.julia_id_v):
        u_structures = nodes_to_structures.get(u, [])
        v_structures = nodes_to_structures.get(v, [])
        if (
            u_structures
            and v_structures
            and not set(u_structures).intersection(v_structures)
        ):
            for us in u_structures:
                for vs in v_structures:
                    G.add_edge(us, vs)
    wG = nx.Graph()
    wG.add_nodes_from(G.nodes())
    wG.add_weighted_edges_from([(*k, v) for k, v in dict(Counter(G.edges())).items()])
    return wG


def get_nodes_to_structures(nodes):
    nodes_to_structures = {}
    for idx, nodeset in enumerate(nodes, start=1):
        for node in nodeset:
            nodes_to_structures[node] = nodes_to_structures.get(node, []) + [idx]
    return nodes_to_structures


def get_node_color(node_type):
    if node_type == "star":
        return "orange"
    elif node_type == "clique":
        return "dodgerblue"
    elif node_type == "biclique":
        return "#BE271A"  # red3
    elif node_type == "starclique":
        return "orchid"
    else:
        raise


def get_node_marker(node_type):
    if node_type == "star":
        return "^"
    elif node_type == "clique":
        return "o"
    elif node_type == "biclique":
        return "s"
    elif node_type == "starclique":
        return "d"
    else:
        raise


def plot_structure_quotient_graph(wG, nodes, structures, save_path=None):
    pos = nx.layout.fruchterman_reingold_layout(wG, k=2.5, seed=0)
    _ = plt.figure(figsize=(12, 12))
    nx.draw_networkx_edges(
        wG,
        pos=pos,
        edgelist=wG.edges(),
        width=[w / 100 for u, v, w in wG.edges(data="weight")],
    )
    for node in wG.nodes():
        plt.scatter(
            *pos[node],
            s=len(nodes[node - 1]) * 5,
            c=get_node_color(structures[node - 1]["structure_type"]),
            marker=get_node_marker(structures[node - 1]["structure_type"]),
        )
    nx.draw_networkx_labels(wG, pos, zorder=100)
    plt.axis("off")
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path)
        plt.close()


def create_overlap_quotient_graph(structures_added, abs_df, n_total):
    G = nx.Graph()
    for idx, structure in enumerate(structures_added, start=1):
        G.add_node(
            idx, **{**structure, "n_relative": structure["n_nodes_total"] / n_total}
        )
    for i in range(len(abs_df)):
        for j in range(i + 1, len(abs_df)):
            edge_weight = abs_df.at[i, j] / n_total
            if edge_weight > 0:
                G.add_edge(i + 1, j + 1, weight=edge_weight)
    return G


def plot_overlap_quotient_graph(G, df, n_total, save_path=None):
    np.random.seed(1234)
    pos = nx.layout.fruchterman_reingold_layout(G)
    _, ax = plt.subplots(1, 1, figsize=(12, 12))
    for x, y, w in G.edges(data="weight"):
        if w * n_total > 1:
            ax.plot(
                [pos[x][0], pos[y][0]],
                [pos[x][1], pos[y][1]],
                "-k",
                linewidth=w * n_total / 100,
                zorder=-10,
                alpha=0.5,
            )
    for node in G.nodes(data=True):
        ax.scatter(
            *pos[node[0]],
            s=5.0 * node[1]["n_nodes_total"],
            c=df.at[node[0] + 1, "structure_color"],
            marker=df.at[node[0] + 1, "structure_shape"],
            zorder=1,
        )
    nx.draw_networkx_labels(G, pos, zorder=100)
    plt.axis("off")
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, transparent=True)
        plt.close()
