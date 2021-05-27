import os
import re
from datetime import datetime

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

from custom_utils import get_dataset_color

sns.set_style("whitegrid")


def parse_time(x):
    try:
        return datetime.strptime(x, "%Y-%m-%dT%H:%M:%S.%f")
    except:
        return datetime.strptime(x, "%Y-%m-%dT%H:%M:%S")


def get_stats(log, path):
    with open(f"{path}/{log}") as f:
        lines = f.readlines()
    ct = parse_time(lines[-1].split()[0]) - parse_time(lines[0].split()[0])
    graph_size = re.search(
        r"(\d+)\snodes.+\s(\d+)\sedges",
        [l for l in lines if "Created graph" in l][0].split(" ", 2)[-1],
    )
    n = int(graph_size.group(1))
    m = int(graph_size.group(2))
    n_successes = len([l for l in lines if "success:" in l])
    n_nonstars = len(
        [
            l
            for l in lines
            if "added: Clique " in l
            or "added: Biclique " in l
            or "added: Starclique " in l
        ]
    )
    n_failures = len([n for l in lines if "failure:" in l])
    return [n, m, n_successes, n_nonstars, n_failures, ct]


def get_collection(x):
    if x.startswith("tissue"):
        return "bio"
    if x.startswith("cs.LG"):
        return "clg"
    if x.startswith("cs.SI"):
        return "csi"
    if x.startswith("oregon1"):
        return "asb"
    if x.startswith("oregon2"):
        return "asp"
    if x.startswith("us_"):
        return "lus"
    if x.startswith("de_"):
        return "lde"
    if x.startswith("ba-"):
        return "rba"
    if x.startswith("er-"):
        return "rer"
    else:
        raise


def get_performance_df(paths, all_logs):
    df = pd.DataFrame(
        columns=["n", "m", "n_structures", "n_nonstars", "n_failures", "ct", "dataset"]
    )
    for path, logs in zip(paths, all_logs):
        for log in logs:
            df.loc[log] = get_stats(log, path) + [path.split("/")[-1]]
    df.n = df.n.astype(int)
    df.m = df.m.astype(int)
    df.n_structures = df.n_structures.astype(int)
    df.n_failures = df.n_failures.astype(int)
    df.ct = df.ct.astype("timedelta64[s]").map(lambda x: max(int(x), 1))
    df["collection"] = df.index.map(get_collection)
    df["n_tries"] = df.n_structures + df.n_failures
    df["n_plus_m"] = df.n + df.m
    return df


paths = [
    f"../results/{x}" for x in ["as-data", "law", "arxiv-temporal", "bio", "random"]
]
logs = [sorted([f for f in os.listdir(path) if f.endswith("log")]) for path in paths]

df = get_performance_df(paths, logs)
df = df.rename(dict(collection="Collection", n_structures="|S|"), axis=1)

fontsize = 36
patch = mpatches.Patch(color="w", label="")
cpatch = Line2D(
    range(1),
    range(1),
    color="white",
    marker="o",
    markersize=15,
    markerfacecolor="slategray",
)
f, ax = plt.subplots(figsize=(15, 9))
ax.set(xscale="log", yscale="log")
sns.scatterplot(
    data=df,
    x="m",
    y="ct",
    size="|S|",
    sizes=(20, 200),
    hue="Collection",
    palette={x: get_dataset_color(x) for x in df["Collection"].unique()},
)
x = np.linspace(2 ** 10.75, 2 ** 19.25, 1000)
ax.plot(x, x / 25, color="k")
plt.text(2 ** 16.5, 2 ** 14.25, r"$f(x) = \frac{1}{25}\cdot x$", fontsize=fontsize - 10)
handles, labels = ax.get_legend_handles_labels()
plt.xlim(2 ** 10.75, 2 ** 19.25)
plt.ylim(2 ** 2.75, 2 ** 15.25)
plt.yticks(
    [2 ** x for x in range(3, 16)],
    [f"$2^{'{'+str(x)+'}'}$" for x in range(3, 16)],
    fontsize=fontsize,
)
plt.xticks(
    [2 ** x for x in range(11, 20)],
    [f"$2^{'{'+str(x)+'}'}$" for x in range(11, 20)],
    fontsize=fontsize,
)
plt.xlabel("Number of edges", fontsize=fontsize + 2)
plt.ylabel("Computation time in seconds", fontsize=fontsize + 2)
plt.legend(
    handles + [patch, patch, patch, patch],
    labels + ["", "", "", "", ""],
    columnspacing=0,
    fontsize=fontsize,
    ncol=2,
    labelspacing=0.2,
    bbox_to_anchor=(1, 1),
    frameon=False,
    handletextpad=0.01,
    borderpad=0,
)
plt.tight_layout()
plt.savefig("../graphics/other/performance.pdf", transparent=True, bbox_inches="tight")
