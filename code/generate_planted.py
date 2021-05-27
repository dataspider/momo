import itertools
import os

import networkx as nx

from custom_utils import get_node_color


def create_biclique(l, r):
    G = nx.Graph()
    G.add_nodes_from(list(range(l + r)))
    edges = list(itertools.product(list(range(l)), list(range(l, l + r))))
    G.add_edges_from(edges)
    return G


def create_starclique(l, r):
    G = nx.Graph()
    G.add_nodes_from(list(range(l + r)))
    edges = list(itertools.product(list(range(l)), list(range(l, l + r))))
    G.add_edges_from(edges)
    G.add_edges_from(
        [
            (u, v)
            for u, v in list(itertools.product(list(range(l)), list(range(l))))
            if u != v
        ]
    )
    return G


path = "../data/planted"
if not os.path.exists(path):
    os.makedirs(path)

# larger graph
star = nx.star_graph(19)
clique = nx.complete_graph(20)
biclique = create_biclique(10, 10)
starclique = create_starclique(10, 10)

G = nx.disjoint_union_all([star, clique, biclique, starclique])
G = nx.contracted_nodes(G, len(star) - 1, len(star))
G = nx.contracted_nodes(G, len(star) + len(clique) - 2, len(star) + len(clique))
G = nx.contracted_nodes(G, 0, len(G) - 1)
G.add_edges_from([(6, 7), (70, 71), (29, len(G)), (70, len(G)), (41, 100), (11, 101)])

nx.write_edgelist(G, f"{path}/dummy_1-1-1-1.txt", data=False)

# smaller graph
star = nx.star_graph(14)
clique = nx.complete_graph(15)
starclique = create_starclique(5, 10)

G = nx.disjoint_union_all([star, clique, starclique])
G = nx.contracted_nodes(G, 3, 43)
G.add_edges_from([(44, 26), (18, 100), (19, 13)])

nx.write_edgelist(G, f"{path}/dummy_1-1-0-1.txt", data=False)
