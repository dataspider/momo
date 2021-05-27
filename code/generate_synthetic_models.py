import itertools
import json
import math
import os

import numpy as np


def get_max_n_edges(n_nodes, n_other_nodes=None):
    if n_other_nodes is not None:
        return n_nodes * n_other_nodes
    else:
        return int(n_nodes * (n_nodes - 1) / 2)


def get_n_edges(n_nodes, density, n_other_nodes=None):
    return int(math.floor(density * get_max_n_edges(n_nodes, n_other_nodes)))


def get_n_nodes(n, fraction):
    return int(math.floor(fraction * n))


def get_clique(n_nodes, n_edges):
    return {
        "structure_type": "clique",
        "n_nodes_total": n_nodes,
        "n_edges_total": n_edges,
    }


def get_star(n_nodes_in_spokes, n_edges_in_spokes):
    return {
        "structure_type": "star",
        "n_nodes_total": n_nodes_in_spokes + 1,
        "n_edges_total": n_edges_in_spokes + n_nodes_in_spokes,
        "n_nodes_in_spokes": n_nodes_in_spokes,
        "n_edges_in_spokes": n_edges_in_spokes,
    }


def get_generic_bipartite(
    structure_type,
    n_nodes_in_left,
    n_nodes_in_right,
    n_edges_in_left,
    n_edges_in_right,
    n_edges_across,
):
    assert structure_type in ["biclique", "starclique"]
    return {
        "structure_type": structure_type,
        "n_nodes_total": n_nodes_in_left + n_nodes_in_right,
        "n_edges_total": n_edges_in_left + n_edges_in_right + n_edges_across,
        "n_nodes_in_right": n_nodes_in_right,
        "n_nodes_in_left": n_nodes_in_left,
        "n_edges_in_right": n_edges_in_right,
        "n_edges_in_left": n_edges_in_left,
        "n_edges_across": n_edges_across,
    }


def get_model(n: int, m: int, macro_structures: list, parameters: dict):
    return {
        "n": n,
        "m": m,
        "macro_structures": macro_structures,
        "parameters": parameters,
    }


def get_total_edges(structures):
    return sum(map(lambda x: x["n_edges_total"], structures))


class Model:
    def __init__(self, **entries):
        self.__dict__.update(entries)

    def __repr__(self):
        return f"<Model for graph with n={self.n}, m={self.m}, edge noise={self.parameters['noise_edge_fraction']}, |cl|={self.parameters['n_cliques']}>"

    def save(self, path):
        filename = (
            f"synthetic_cl{self.parameters.get('n_cliques', 0):03}_st{self.parameters.get('n_stars', 0):03}_bc{self.parameters.get('n_bicliques', 0):03}_sc{self.parameters.get('n_starcliques', 0):03}"
            + f"_n{self.n:015}_m{self.m:015}_noise{str(self.parameters['noise_edge_fraction']).replace('.', '-')}.json"
        )
        if not os.path.exists(path):
            os.makedirs(path)
        with open(f"{path}/{filename}", "w") as f:
            json.dump(self.__dict__, f, indent=4)


def generate_model_from_structures(
    n, structures, noise_edge_fraction=0.05, generating_function=None
):
    m = get_total_edges(structures) + math.floor(
        noise_edge_fraction * get_total_edges(structures)
    )
    n_cliques, n_stars, n_bicliques, n_starcliques = [
        len(list(filter(lambda x: x["structure_type"] == y, structures)))
        for y in ["clique", "star", "biclique", "starclique"]
    ]
    M = get_model(
        n,
        m,
        structures,
        parameters=dict(
            noise_edge_fraction=noise_edge_fraction,
            n_cliques=n_cliques,
            n_stars=n_stars,
            n_bicliques=n_bicliques,
            n_starcliques=n_starcliques,
            generating_function=generating_function.__name__
            if generating_function is not None
            else generate_model_from_structures.__name__,
        ),
    )
    return M


def create_clique_structures(n, n_cliques):
    fractions = np.linspace(0.001, 0.01, n_cliques)
    densities = np.linspace(0.9, 0.5, n_cliques)
    cliques = [
        get_clique(
            get_n_nodes(n, fraction), get_n_edges(get_n_nodes(n, fraction), density)
        )
        for fraction, density in zip(fractions, densities)
    ]
    return cliques


def create_star_structures(n, n_stars):
    fractions = np.linspace(0.001, 0.01, n_stars)
    densities = np.linspace(0.01, 0.05, n_stars)
    stars = [
        get_star(
            get_n_nodes(n, fraction), get_n_edges(get_n_nodes(n, fraction), density)
        )
        for fraction, density in zip(fractions, densities)
    ]
    return stars


def create_generic_bipartite_structures(n, n_bicliques=0, n_starcliques=0):
    bicliques, starcliques = [], []
    if n_starcliques > 0:
        l_fractions = np.linspace(0.001, 0.01, n_starcliques)
        r_fractions = np.linspace(0.01, 0.05, n_starcliques)
        l_ds = np.linspace(0.9, 0.5, n_starcliques)
        r_ds = np.linspace(0.001, 0.05, n_starcliques)
        a_ds = np.linspace(0.9, 0.5, n_starcliques)
        starcliques = [
            get_generic_bipartite(
                "starclique",
                get_n_nodes(n, l_fraction),
                get_n_nodes(n, r_fraction),
                get_n_edges(get_n_nodes(n, l_fraction), l_d),
                get_n_edges(get_n_nodes(n, r_fraction), r_d),
                get_n_edges(
                    get_n_nodes(n, l_fraction), a_d, get_n_nodes(n, r_fraction)
                ),
            )
            for l_fraction, r_fraction, l_d, r_d, a_d in zip(
                l_fractions, r_fractions, l_ds, r_ds, a_ds
            )
        ]
    if n_bicliques > 0:
        l_fractions = np.linspace(0.0001, 0.001, n_bicliques)
        r_fractions = np.linspace(0.001, 0.01, n_bicliques)
        l_ds = np.linspace(0.001, 0.05, n_bicliques)
        r_ds = np.linspace(0.001, 0.05, n_bicliques)
        a_ds = np.linspace(0.9, 0.5, n_bicliques)
        bicliques = [
            get_generic_bipartite(
                "biclique",
                get_n_nodes(n, l_fraction),
                get_n_nodes(n, r_fraction),
                get_n_edges(get_n_nodes(n, l_fraction), l_d),
                get_n_edges(get_n_nodes(n, r_fraction), r_d),
                get_n_edges(
                    get_n_nodes(n, l_fraction), a_d, get_n_nodes(n, r_fraction)
                ),
            )
            for l_fraction, r_fraction, l_d, r_d, a_d in zip(
                l_fractions, r_fractions, l_ds, r_ds, a_ds
            )
        ]
    return bicliques, starcliques


def nonempty_powerset(iterable):
    "nonempty_powerset([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(
        itertools.combinations(s, r) for r in range(1, len(s) + 1)
    )


n_structures = 100
for n in np.linspace(100000, 1000000, 10):
    n = int(n)
    cliques = create_clique_structures(n, n_structures)
    stars = create_star_structures(n, n_structures)
    bicliques, _ = create_generic_bipartite_structures(
        n, n_bicliques=n_structures, n_starcliques=0
    )
    _, starcliques = create_generic_bipartite_structures(
        n, n_bicliques=0, n_starcliques=n_structures
    )
    for combo in nonempty_powerset([cliques, stars, bicliques, starcliques]):
        structures = []
        cutoff = int(math.floor(100 / len(combo)))
        for elem in combo:
            structures.extend(elem[:cutoff])
        M = generate_model_from_structures(
            n, structures, noise_edge_fraction=0.05, generating_function=None
        )
        Model(**M).save("../results/synthetic")
