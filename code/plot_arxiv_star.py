import textwrap

import matplotlib.pyplot as plt
import networkx as nx

from custom_utils import get_structures_added, load_json


def get_name(x):
    return " ".join(map(lambda x: x.capitalize(), x.split("_")[1:] + [x.split("_")[0]]))


years = list(range(2011, 2021))

csis = [
    get_structures_added(
        load_json(f"../results_mapped/arxiv-temporal/cs.SI_{x}-11-01.json")
    )
    for x in years
]

G = nx.read_edgelist("../data/arxiv-temporal/cs.SI_2020-11-01.txt")

hub = "faloutsos_christos"
fc_star = list(map(lambda x: [c for c in x if c.get("hub", "").startswith(hub)], csis))[
    -1
][0]
fc_nodes = [fc_star["hub"]] + fc_star["spokes"]
Gsub = G.subgraph(fc_nodes)
# damn middle names
nodes_to_contract = [
    ["almeida_jussara_m", "almeida_jussara"],
    ["rodrigues_jose", "rodrigues_jose_f"],
    ["traina_agma", "traina_agma_j_m"],
]
for nodes in nodes_to_contract:
    Gsub = nx.minors.contracted_nodes(Gsub, *nodes)
fig, ax = plt.subplots(figsize=(12, 12))
pos = nx.fruchterman_reingold_layout(Gsub, seed=12, scale=0.8)
nx.draw_networkx_nodes(Gsub, pos=pos, node_color="silver")
nx.draw_networkx_edges(Gsub, pos=pos)
for x in Gsub.nodes():
    if x != hub:
        nx.draw_networkx_labels(
            Gsub,
            pos={
                k: (
                    [
                        v[0] - 0.02 if v[0] > 0 else v[0] + 0.02,
                        v[1] - 0.02
                        if k.startswith("dighe")
                        else (v[1] + 0.02 if k.startswith("kova") else v[1]),
                    ]
                )
                for k, v in pos.items()
            },
            labels={
                x: textwrap.fill(get_name(x), 10)
                if x != "eliassirad_tina"
                else textwrap.fill("Tina Eliassi-Rad", 15)
            },
            font_size=22,
            font_weight="bold",
            font_color="k",
        )
    else:
        nx.draw_networkx_labels(
            Gsub,
            pos=pos,
            labels={
                x: textwrap.fill(get_name(x), 15) for x in Gsub.nodes() if x == hub
            },
            font_size=34,
            font_color="k",
            font_weight="bold",
        )
plt.axis("off")
plt.tight_layout()
plt.savefig("../graphics/other/star-cs-si-2020.pdf")
