import os
import textwrap
from copy import deepcopy
from itertools import product

import matplotlib.pyplot as plt
import pandas as pd

import squarify
from custom_utils import get_node_color, get_structures_added, load_json


def get_name(x):
    return " ".join(map(lambda x: x.capitalize(), x.split("_")[1:] + [x.split("_")[0]]))


def get_single_label(x, stype, gstructures):
    if stype != "star":
        return ""
    else:
        return textwrap.fill(get_name(gstructures[x]["hub"]), 15)


def get_common_label(x, y, stype, g1structures, g2structures):
    if stype != "star":
        return ""
    else:
        return (
            textwrap.fill(get_name(g1structures[x]["hub"]), 15)
            + "\n"
            + "→"
            + "\n"
            + textwrap.fill(get_name(g2structures[y]["hub"]), 15)
        )


def get_hatch(structure_type):
    if structure_type == "clique":
        return "..."
    elif structure_type == "star":
        return "xxx"
    elif structure_type == "biclique":
        return "///"
    elif structure_type == "starclique":
        return "\\\\\\"


def plot_treemaps(g1, g2, alignment, jsons):
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["font.size"] = 20
    alignment_1 = f"{g1}_-_{g2}"
    alignment_2 = f"{g2}_-_{g1}"
    if alignment_1 not in alignment.alignment.unique():
        g1, g2 = g2, g1
    sdf = deepcopy(alignment.query("g1 == @g1 and g2 == @g2"))

    g1json = jsons[g1]
    g2json = jsons[g2]
    n1 = g1json["n"]
    n2 = g2json["n"]
    m1 = g1json["m"]
    m2 = g2json["m"]
    g1structures = {x["position"]: x for x in get_structures_added(g1json)}
    g2structures = {x["position"]: x for x in get_structures_added(g2json)}

    sdf["g1_structure_n"] = [
        g1structures[x]["n_nodes_total"] / n1 if not pd.isna(x) else 0
        for x in sdf.g1_structure
    ]
    sdf["g2_structure_n"] = [
        g2structures[x]["n_nodes_total"] / n2 if not pd.isna(x) else 0
        for x in sdf.g2_structure
    ]
    sdf["structure_fraction_n"] = (sdf.g1_structure_n + sdf.g2_structure_n) / 2
    sdf["structure_type"] = [
        x.split("_")[-2] if not pd.isna(x) else y.split("_")[-2]
        for x, y in zip(sdf.g1_structure, sdf.g2_structure)
    ]

    g1values = sdf.query("not @pd.isna(g1_structure)").g1_structure_n.values
    g2values = sdf.query("not @pd.isna(g2_structure)").g2_structure_n.values
    g12values = sdf.query(
        "not @pd.isna(g1_structure) and not @pd.isna(g2_structure)"
    ).structure_fraction_n.values
    g1labels = (
        []
    )  # list(map(lambda x: int(x.split("_")[-1]), sdf.query("not @pd.isna(g1_structure)").g1_structure))
    g2labels = (
        []
    )  # list(map(lambda x: int(x.split("_")[-1]), sdf.query("not @pd.isna(g2_structure)").g2_structure))
    g12labels = map(
        lambda x: round(x, 2),
        sdf.query(
            "not @pd.isna(g1_structure) and not @pd.isna(g2_structure)"
        ).jaccard.values,
    )
    g1densities = [
        0.4 if not pd.isna(x) else 1
        for x in sdf.query("not @pd.isna(g1_structure)").jaccard
    ]
    g2densities = [
        0.4 if not pd.isna(x) else 1
        for x in sdf.query("not @pd.isna(g2_structure)").jaccard
    ]
    g12densities = [1] * len(
        sdf.query("not @pd.isna(g1_structure) and not @pd.isna(g2_structure)")
    )  # [max(x,0.15) for x in sdf.query("not @pd.isna(g1_structure) and not @pd.isna(g2_structure)").jaccard]
    g1colors = list(
        map(
            get_node_color,
            sdf.query("not @pd.isna(g1_structure)").structure_type.values,
        )
    )
    g1colors = [c if d == 1 else "black" for c, d in zip(g1colors, g1densities)]
    g2colors = list(
        map(
            get_node_color,
            sdf.query("not @pd.isna(g2_structure)").structure_type.values,
        )
    )
    g2colors = [c if d == 1 else "black" for c, d in zip(g2colors, g2densities)]
    g12colors = list(
        map(
            get_node_color,
            sdf.query(
                "not @pd.isna(g1_structure) and not @pd.isna(g2_structure)"
            ).structure_type.values,
        )
    )
    g1hatches = list(
        map(get_hatch, sdf.query("not @pd.isna(g1_structure)").structure_type.values)
    )
    g2hatches = list(
        map(get_hatch, sdf.query("not @pd.isna(g2_structure)").structure_type.values)
    )
    g12hatches = list(
        map(
            get_hatch,
            sdf.query(
                "not @pd.isna(g1_structure) and not @pd.isna(g2_structure)"
            ).structure_type.values,
        )
    )
    g12norm = 100
    g1norm = 100  # n1/((n1 + n2)/2)*100
    g2norm = 100  # n2/((n1 + n2)/2)*100

    fig, ax = plt.subplots(figsize=(10, 10))  # (g1norm/100*12,g1norm/100*12))
    squarify.plot(
        sizes=g1values,
        pad=True,
        label=g1labels,
        color=g1colors,
        norm_x=g1norm,
        norm_y=g1norm,
        ax=ax,
    )
    for patch, d, h in zip(ax.patches, g1densities, g1hatches):
        patch.set_alpha(d)
    #         patch.set_hatch(h)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(
        f"../graphics/treemaps/{g1}_-_{g2}_M1.pdf",
        transparent=True,
        bbox_inches="tight",
    )
    plt.close()

    fig, ax = plt.subplots(figsize=(g2norm / 100 * 12, g2norm / 100 * 12))
    squarify.plot(
        sizes=g2values,
        pad=True,
        label=g2labels,
        color=g2colors,
        norm_x=g2norm,
        norm_y=g2norm,
        ax=ax,
    )
    for patch, d, h in zip(ax.patches, g2densities, g2hatches):
        patch.set_alpha(d)
    #         patch.set_hatch(h)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(
        f"../graphics/treemaps/{g1}_-_{g2}_M2.pdf",
        transparent=True,
        bbox_inches="tight",
    )
    plt.close()

    fig, ax = plt.subplots(figsize=(12, 12))
    squarify.plot(
        sizes=g12values,
        pad=True,
        label=g12labels,
        color=g12colors,
        norm_x=g12norm,
        norm_y=g12norm,
        ax=ax,
        text_kwargs={"size": "40", "weight": "bold", "fontstretch": "ultra-condensed"},
    )
    plt.axis("off")
    for patch, d, h in zip(ax.patches, g12densities, g12hatches):
        patch.set_alpha(d)
    #         patch.set_hatch(h)
    plt.tight_layout()
    plt.savefig(
        f"../graphics/treemaps/{g1}_-_{g2}_M12.pdf",
        transparent=True,
        bbox_inches="tight",
    )
    plt.close()


path = "../results_mapped/bio"
alignment = pd.read_csv(f"{path}/structure-alignments-jaccard.csv")
jsons = {
    p.split(".json")[0]: load_json(f"{path}/{p}")
    for p in os.listdir(path)
    if p.endswith("json")
}
alignment["g1"], alignment["g2"] = list(
    zip(*[x.split("_-_") for x in alignment.alignment])
)
first_row = [
    f"tissue-{x}"
    for x in ["colon", "small_intestine", "large_intestine", "pancreas", "liver"]
]
second_row = [
    f"tissue-{x}"
    for x in ["stomach", "vermiform_appendix", "cecum", "duodenum", "esophagus"]
]
digestive_tract = first_row + second_row

combinations = sorted(
    {tuple(sorted(x)) for x in product(digestive_tract, digestive_tract)}
)

if not os.path.exists("../graphics/treemaps"):
    os.makedirs("../graphics/treemaps")
for combination in combinations:
    plot_treemaps(*combination, alignment, jsons)
