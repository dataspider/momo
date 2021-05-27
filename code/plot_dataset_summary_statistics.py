import os

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns

from custom_utils import get_dataset_color
from statics import DATA_PATHS_FOR_ABBRS

sns.set_style("whitegrid")

df = pd.DataFrame(columns=["collection", "n", "m"])
for dataset, path in DATA_PATHS_FOR_ABBRS.items():
    print(dataset, end="\r")
    if dataset == "asb":
        files = sorted(
            [x for x in os.listdir(path) if "oregon1" in x and x.endswith("txt")]
        )
    elif dataset == "asp":
        files = sorted(
            [x for x in os.listdir(path) if "oregon2" in x and x.endswith("txt")]
        )
    elif dataset == "lus":
        files = sorted(
            [x for x in os.listdir(path) if "us_" in x and x.endswith("txt")]
        )
    elif dataset == "lde":
        files = sorted(
            [x for x in os.listdir(path) if "de_" in x and x.endswith("txt")]
        )
    elif dataset == "clg":
        files = sorted(
            [x for x in os.listdir(path) if "cs.LG" in x and x.endswith("txt")]
        )
    elif dataset == "csi":
        files = sorted(
            [x for x in os.listdir(path) if "cs.SI" in x and x.endswith("txt")]
        )
    elif dataset == "rba":
        files = sorted(
            [x for x in os.listdir(path) if "ba-" in x and x.endswith("txt")]
        )
    elif dataset == "rer":
        files = sorted(
            [x for x in os.listdir(path) if "er-" in x and x.endswith("txt")]
        )
    else:
        files = sorted([x for x in os.listdir(path) if x.endswith("txt")])
    for file in files:
        G = nx.read_edgelist(f"{path}/{file}")
        df.loc[file] = [dataset, G.number_of_nodes(), G.number_of_edges()]
df.n = df.n.astype(float)
df.m = df.m.astype(float)

descriptions = [
    (dataset, df.query("collection == @collection"))
    for collection in DATA_PATHS_FOR_ABBRS.keys()
]


fontsize = 28

for name, plot_func in [("boxen", sns.boxenplot)]:
    plt.figure(figsize=(8, 6))
    plot_func(
        data=df,
        x="collection",
        y="n",
        palette={x: get_dataset_color(x) for x in df.collection.unique()},
        hue_order=[],
    )
    plt.xlabel("Collection", fontsize=fontsize + 4)
    plt.ylabel("Number of Nodes", fontsize=fontsize + 4)
    plt.yticks(
        np.arange(0, 160001, 20000),
        [f"{x//1000}K" for x in np.arange(0, 160001, 20000)],
        fontsize=fontsize,
    )
    plt.ylim(-100, 160000)
    plt.xticks(plt.gca().get_xticks(), fontsize=fontsize)
    plt.tight_layout()
    plt.savefig(f"../graphics/other/figure-collection-statistics-{name}-n.pdf")
    plt.close()

for name, plot_func in [("boxen", sns.boxenplot)]:
    plt.figure(figsize=(8, 6))
    plot_func(
        data=df,
        x="collection",
        y="m",
        palette={x: get_dataset_color(x) for x in df.collection.unique()},
    )
    plt.xlabel("Collection", fontsize=fontsize + 4)
    plt.ylabel("Number of Edges", fontsize=fontsize + 4)
    plt.yticks(
        np.arange(0, 525001, 75000),
        [f"{x//1000}K" for x in np.arange(0, 525001, 75000)],
        fontsize=fontsize,
    )
    plt.xticks(plt.gca().get_xticks(), fontsize=fontsize)
    plt.ylim(-100, 525000)
    plt.tight_layout()
    plt.savefig(f"../graphics/other/figure-collection-statistics-{name}-m.pdf")
    plt.close()
