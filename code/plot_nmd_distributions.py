import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from custom_utils import get_dataset_color

sns.set_style("whitegrid")
heuristic = "overlap"
descriptions = [
    (
        "asb",
        pd.read_csv(f"../results_mapped/as-data/nmd-{heuristic}.csv").query(
            'g1 != g2 and g1.str.startswith("oregon1") and g2.str.startswith("oregon1")'
        ),
    ),
    (
        "asp",
        pd.read_csv(f"../results_mapped/as-data/nmd-{heuristic}.csv").query(
            'g1 != g2 and g1.str.startswith("oregon2") and g2.str.startswith("oregon2")'
        ),
    ),
    ("bio", pd.read_csv(f"../results_mapped/bio/nmd-jaccard.csv").query("g1 != g2")),
    (
        "clg",
        pd.read_csv(f"../results_mapped/arxiv-temporal/nmd-{heuristic}.csv").query(
            'g1 != g2 and g1.str.startswith("cs.LG") and g2.str.startswith("cs.LG")'
        ),
    ),
    (
        "csi",
        pd.read_csv(f"../results_mapped/arxiv-temporal/nmd-{heuristic}.csv").query(
            'g1 != g2 and g1.str.startswith("cs.SI") and g2.str.startswith("cs.SI")'
        ),
    ),
    (
        "lde",
        pd.read_csv(f"../results_mapped/law/nmd-{heuristic}.csv").query(
            'g1 != g2 and g1.str.startswith("de_") and g2.str.startswith("de_")'
        ),
    ),
    (
        "lus",
        pd.read_csv(f"../results_mapped/law/nmd-{heuristic}.csv").query(
            'g1 != g2 and g1.str.startswith("us_") and g2.str.startswith("us_")'
        ),
    ),
]
dfs = pd.concat(
    [
        pd.DataFrame({"dataset": [dataset] * len(nmds), "nmd": nmds["nmd"]})
        for dataset, nmds in descriptions
    ]
)
number_of_observations = dfs.groupby("dataset").count().reset_index()
number_of_observations[
    "label"
] = number_of_observations.dataset  # + "\nn=" + number_of_observations.nmd.map(str)
number_of_observations.set_index("dataset", inplace=True)
dfs.dataset = dfs.dataset.map(number_of_observations.label)

palette = {k: get_dataset_color(k) for k in set(dfs.dataset.values)}
fontsize = 30
for name, plot_func in [("box", sns.boxplot)]:
    plt.figure(figsize=(9, 6))
    plot_func(data=dfs, x="dataset", y="nmd", palette=palette, linewidth=3)
    plt.xlabel("", fontsize=fontsize)
    plt.ylabel("NMD", fontsize=fontsize + 6)
    plt.yticks(np.arange(0, 1.01, 0.1), fontsize=fontsize)
    plt.xticks(plt.gca().get_xticks(), fontsize=fontsize + 6)
    plt.ylim(-0.05, 1.05)
    plt.tight_layout()
    plt.savefig(
        f"../graphics/other/{heuristic}-{name}.pdf",
        transparent=True,
        bbox_inches="tight",
    )
    plt.close()
