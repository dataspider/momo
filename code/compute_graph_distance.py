import argparse
import math
import os

import pandas as pd

from custom_utils import get_structures_added, load_json
from statics import MAPPED_RESULT_PATH


def compute_nmd(individual_model_cost1, individual_model_cost2, common_model_cost):
    if  max(individual_model_cost1, individual_model_cost2) > 0:
        nmd_when_encoding_together = (common_model_cost - min(individual_model_cost1, individual_model_cost2)) / max(individual_model_cost1, individual_model_cost2)
        return min(nmd_when_encoding_together, 1.)
    else:
        return 0.

def get_structures_indexed(json, mapped_result_path):
    if "synthetic" in mapped_result_path:
        structures_added = json['macro_structures']
    else:
        structures_added = get_structures_added(json)
    structures = {s['position']:s for s in structures_added}
    return structures

def compute_graph_distance_from_alignments(alignment_file, mapped_result_path):
    json_paths = {x.split(".json")[0]: f"{mapped_result_path}/{x}" for x in [f for f in os.listdir(mapped_result_path) if f.endswith('json')]}
    alignments = pd.read_csv(alignment_file)
    alignments['g1'] = alignments.alignment.map(lambda x: x.split("_-_")[0])
    alignments['g2'] = alignments.alignment.map(lambda x: x.split("_-_")[1])
    nmds = pd.DataFrame(columns=["g1","g2","nmd"])
    unique_alignments = alignments.alignment.unique()
    n_unique_alignments = len(unique_alignments)
    for idx,alignment in enumerate(unique_alignments):
        print(f"Starting alignment {idx+1:05} of {n_unique_alignments:05}", end="\r")
        reversed_alignment = False
        g1, g2 = alignment.split("_-_")
        json1 = load_json(json_paths[g1])
        json2 = load_json(json_paths[g2])
        n1, n2 = json1['n'], json2['n']
        if n1 < n2:
            n1, n2 = n2, n1
            g1, g2 = g2, g1
            json1, json2 = json2, json1
            reversed_alignment = True

        structures1 = get_structures_indexed(json1, mapped_result_path)
        structures2 = get_structures_indexed(json2, mapped_result_path)
        for s in structures1.values():
            s['individual_cost'] = get_unpaired_structures_cost([s], n1)
        for s in structures2.values():
            s['individual_cost'] = get_unpaired_structures_cost([s], n2)

        df = alignments.query("alignment == @alignment")
        matched_structures = df.query("not @pd.isna(g1_structure) and not @pd.isna(g2_structure)")
        if not reversed_alignment:
            paired_structures = list(map(
                lambda x: (structures1[x[0]], structures2[x[1]]),
                list(zip(matched_structures.g1_structure, matched_structures.g2_structure))))
            unpaired_structures1 = list(map(
                lambda x: structures1[x],
                list(df.query("@pd.isna(g2_structure)").g1_structure)))
            unpaired_structures2 = list(map(
                lambda x: structures2[x],
                list(df.query("@pd.isna(g1_structure)").g2_structure)))
        else:
            paired_structures = list(map(
                lambda x: (structures1[x[0]], structures2[x[1]]),
                list(zip(matched_structures.g2_structure, matched_structures.g1_structure))))
            unpaired_structures1 = list(map(
                lambda x: structures1[x],
                list(df.query("@pd.isna(g1_structure)").g2_structure)))
            unpaired_structures2 = list(map(
                lambda x: structures2[x],
                list(df.query("@pd.isna(g2_structure)").g1_structure)))

        common_paired_structures_cost = get_paired_structures_cost(paired_structures, n1, n2)
        unpaired_structures_cost1 = get_unpaired_structures_cost(unpaired_structures1, n1)
        unpaired_structures_cost2 = get_unpaired_structures_cost(unpaired_structures2, n2)

        if len(paired_structures) > 0:
            paired_structures1, paired_structures2 = [
                list(x) for x in list(zip(*paired_structures))
            ]
        else:
            paired_structures1, paired_structures2 = [], []
        paired_structures_cost1 = get_unpaired_structures_cost(paired_structures1, n1)
        paired_structures_cost2 = get_unpaired_structures_cost(paired_structures2, n2)

        individual_model_cost1 = paired_structures_cost1 + unpaired_structures_cost1
        individual_model_cost2 = paired_structures_cost2 + unpaired_structures_cost2
        common_model_cost = common_paired_structures_cost + unpaired_structures_cost1 + unpaired_structures_cost2

        if not reversed_alignment:
            nmds = nmds.append(pd.DataFrame(
                [[g1, g2, compute_nmd(individual_model_cost1, individual_model_cost2, common_model_cost)]],columns=["g1","g2","nmd"]))
        else:
            nmds = nmds.append(pd.DataFrame(
                [[g2, g1, compute_nmd(individual_model_cost1, individual_model_cost2, common_model_cost)]],
                columns=["g1", "g2", "nmd"]))
    return nmds


def get_unpaired_structures_cost(unpaired_structures, n):
    cost = 0
    for s1 in unpaired_structures:
        s1type = s1["structure_type"]
        if s1type == "clique":
            structure_cost = get_clique_cost(n, n, s1, s1)
        elif s1type == "star":
            structure_cost = get_star_cost(n, n, s1, s1)
        elif s1type in ["biclique", "starclique"]:
            structure_cost = get_biclique_starclique_cost(n, n, s1, s1, s1type)
        else:
            raise NotImplementedError(f"Scenario not implemented: s1 of type {s1type}")
        if structure_cost < 0:
            raise ValueError(s1)
        cost += structure_cost
    return cost


def get_paired_structures_cost(paired_structures, n1, n2):
    cost = 0
    n_delta_directions = 0
    for s1, s2 in paired_structures:
        s1type = s1["structure_type"]
        s2type = s2["structure_type"]
        if s1type == s2type == "clique":
            structure_cost = get_clique_cost(n1, n2, s1, s2)
            n_delta_directions += 1
        elif s1type == s2type == "star":
            structure_cost = get_star_cost(n1, n2, s1, s2)
            n_delta_directions += 1
        elif s1type == s2type in ["biclique", "starclique"]:
            structure_cost = get_biclique_starclique_cost(n1, n2, s1, s2, s1type)
            n_delta_directions += 3
        else:
            raise NotImplementedError(f"Scenario not implemented: s1 of type {s1type}, s2 of type {s2type}")
        if structure_cost < 0:
            raise ValueError(s1, s2)
        cost += structure_cost
    cost += math.log2(max(n_delta_directions,1))
    return cost


def get_star_cost(n1, n2, s1, s2):
    edge_fraction1, edge_fraction2, node_fraction1, node_fraction2 = get_fractions(
        n1, n2, s1, s2, suffix="in_spokes"
    )
    return get_cost_for_node_edge_set(
        node_fraction1,
        node_fraction2,
        edge_fraction1,
        edge_fraction2,
        n1,
        s1,
        dense=False,
    )


def get_clique_cost(n1, n2, s1, s2):
    edge_fraction1, edge_fraction2, node_fraction1, node_fraction2 = get_fractions(
        n1, n2, s1, s2, suffix="total"
    )
    return get_cost_for_node_edge_set(
        node_fraction1,
        node_fraction2,
        edge_fraction1,
        edge_fraction2,
        n1,
        s1,
        dense=True,
    )


def get_biclique_starclique_cost(n1, n2, s1, s2, structure_type):
    edge_fraction1, edge_fraction2, node_fraction1, node_fraction2 = get_fractions(
        n1, n2, s1, s2, suffix="in_left"
    )
    cost_left = get_cost_for_node_edge_set(
        node_fraction1,
        node_fraction2,
        edge_fraction1,
        edge_fraction2,
        n1,
        s1,
        dense=True if structure_type == "starclique" else False,
    )
    edge_fraction1, edge_fraction2, node_fraction1, node_fraction2 = get_fractions(
        n1, n2, s1, s2, suffix="in_right"
    )
    cost_right = get_cost_for_node_edge_set(
        node_fraction1,
        node_fraction2,
        edge_fraction1,
        edge_fraction2,
        n1,
        s1,
        dense=False,
    )

    across_fraction1 = s1["n_edges_across"] * get_maximum_number_of_edges(
        s1["n_nodes_in_left"], s1["n_nodes_in_right"]
    )
    across_fraction2 = s2["n_edges_across"] * get_maximum_number_of_edges(
        s2["n_nodes_in_left"], s2["n_nodes_in_right"]
    )
    across_communicated = round(get_average([across_fraction1, across_fraction2]) * n1)
    across_delta = abs(s1["n_edges_across"] - across_communicated)
    cost_across_edges = math.log2(math.log2(get_maximum_number_of_edges(s2["n_nodes_in_left"], s2["n_nodes_in_right"]))) # bits for number of edges
    cost_across_delta = math.log2(max(get_maximum_number_of_edges(s2["n_nodes_in_left"], s2["n_nodes_in_right"]) - across_communicated,1)) + math.log2(
            math.log2(
                max(get_maximum_number_of_edges(
                    s2["n_nodes_in_left"], s2["n_nodes_in_right"]
                ) - across_communicated, across_communicated
                    )
            )
        ) + math.log2(across_delta)
    cost_across = cost_across_edges + cost_across_delta
    return cost_left + cost_right + cost_across


def get_fractions(n1, n2, s1, s2, suffix):
    node_fraction1 = get_fraction(s1[f"n_nodes_{suffix}"], n1)
    node_fraction2 = get_fraction(s2[f"n_nodes_{suffix}"], n2)
    edge_fraction1 = get_fraction(
        s1[f"n_edges_{suffix}"], get_maximum_number_of_edges(s1[f"n_nodes_{suffix}"])
    )
    edge_fraction2 = get_fraction(
        s2[f"n_edges_{suffix}"], get_maximum_number_of_edges(s2[f"n_nodes_{suffix}"])
    )
    return edge_fraction1, edge_fraction2, node_fraction1, node_fraction2


def get_cost_for_node_edge_set(
    node_fraction1, node_fraction2, edge_fraction1, edge_fraction2, n1, s1, dense=False
):
    node_fraction = get_average([node_fraction1, node_fraction2])
    edge_fraction = get_average([edge_fraction1, edge_fraction2])
    node_fraction_communicated = round(node_fraction * n1)
    edge_fraction_communicated = round(
        edge_fraction * get_maximum_number_of_edges(node_fraction_communicated)
    )
    node_fraction_communicated_delta = abs(
        s1["n_nodes_total"] - node_fraction_communicated
    )
    edge_fraction_communicated_delta = abs(
        s1["n_edges_total"] - edge_fraction_communicated
    )
    cost_node_fraction = math.log2(math.log2(n1)) + math.log2(
        max(node_fraction_communicated, 1)
    )
    if not dense:
        cost_edge_fraction = math.log2(
            max(math.log2(max(get_maximum_number_of_edges(node_fraction_communicated),1)),1)
        ) + math.log2(max(edge_fraction_communicated, 1))
    else:
        cost_edge_fraction = math.log2(
            math.log2(get_maximum_number_of_edges(node_fraction_communicated))
        ) + math.log2(
            max(
                get_maximum_number_of_edges(node_fraction_communicated)
                - edge_fraction_communicated,
                1,
            )
        )
    cost_node_delta = math.log2(math.log2(max(n1-node_fraction_communicated,node_fraction_communicated))) + math.log2(
        max(node_fraction_communicated_delta, 1)
    )
    cost_edge_delta = math.log2(
        math.log2(max(get_maximum_number_of_edges(node_fraction_communicated)-node_fraction_communicated,node_fraction_communicated))
    ) + math.log2(max(edge_fraction_communicated_delta, 1))
    return cost_node_fraction + cost_edge_fraction + cost_node_delta + cost_edge_delta


def get_fraction(concrete_value, max_value):
    return concrete_value / max_value


def get_maximum_number_of_edges(n_nodes, n_other=None):
    if n_other is None:
        return (n_nodes * (n_nodes - 1)) / 2
    else:
        return n_nodes * n_other


def get_average(values: list):
    return sum(values) / len(values)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("collection", type=str)
    parser.add_argument("--variant", "-v", type=str, default='overlap',choices=["overlap","jaccard"])
    args = parser.parse_args()
    dataset = args.collection
    alignment = args.variant
    if alignment == 'jaccard':
        alignment_file = f"../alignments-jaccard/{dataset}/structure-alignments-{alignment}.csv"
    else:
        alignment_file = f"../alignments/{dataset}/structure-alignments-{alignment}.csv"
    if not os.path.exists(alignment_file):
        raise Exception(f"Necessary alignment file does not exist: {alignment_file}")
    mapped_result_path = f"{MAPPED_RESULT_PATH}/{dataset}"
    nmds = compute_graph_distance_from_alignments(alignment_file, mapped_result_path)
    nmds.to_csv(f"{mapped_result_path}/nmd-{alignment}.csv", index=False)
