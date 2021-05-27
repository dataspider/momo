from collections import Counter

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from statics import STRUCTURE_TYPES

sns.set_style("whitegrid")

first_row = [
    f"tissue-{x}"
    for x in ["colon", "small_intestine", "large_intestine", "pancreas", "liver"]
]
second_row = [
    f"tissue-{x}"
    for x in ["stomach", "vermiform_appendix", "cecum", "duodenum", "esophagus"]
]
digestive_tract = first_row + second_row

large_df = pd.read_csv("../results_mapped/bio/structure-alignments-jaccard.csv")
large_df["g1"] = large_df.alignment.map(lambda x: x.split("_-_")[0])
large_df["g2"] = large_df.alignment.map(lambda x: x.split("_-_")[1])

single_tissue = "tissue-esophagus"
pairs = [
    sorted([g1, g2])
    for idx, g1 in enumerate(digestive_tract)
    for g2 in digestive_tract[idx:]
]

common_model_distributions = pd.DataFrame.from_records(
    [
        dict(
            Counter(
                large_df.query(
                    "g1 == @g1 and g2 == @g2 and not @pd.isna(jaccard)"
                ).g1_structure.map(lambda x: x.split("_")[-2])
            )
        )
        for g1, g2 in pairs
    ]
)
common_model_distributions["g1"], common_model_distributions["g2"] = zip(*pairs)
common_model_distributions.index = (
    common_model_distributions["g1"] + "_-_" + common_model_distributions["g2"]
)
common_model_distributions["total_matched"] = common_model_distributions[
    common_model_distributions.columns[:4]
].sum(axis=1)

common_model_distributions["priority"] = [
    0 if g1 == single_tissue and g2 == single_tissue else 1
    for g1, g2 in zip(common_model_distributions.g1, common_model_distributions.g2)
]

g1_model_leftovers = pd.DataFrame.from_records(
    [
        dict(
            Counter(
                large_df.query(
                    "g1 == @g1 and g2 == @g2 and @pd.isna(jaccard) and not @pd.isna(g1_structure)"
                ).g1_structure.map(lambda x: x.split("_")[-2])
            )
        )
        for g1, g2 in pairs
    ]
).fillna(0)
g1_model_leftovers["g1"], g1_model_leftovers["g2"] = zip(*pairs)
g1_model_leftovers.index = g1_model_leftovers["g1"] + "_-_" + g1_model_leftovers["g2"]
g1_model_leftovers["total_matched"] = common_model_distributions[
    common_model_distributions.columns[:4]
].sum(axis=1)
g1_model_leftovers["priority"] = [
    0 if g1 == single_tissue and g2 == single_tissue else 1
    for g1, g2 in zip(g1_model_leftovers.g1, g1_model_leftovers.g2)
]

g2_model_leftovers = pd.DataFrame.from_records(
    [
        dict(
            Counter(
                large_df.query(
                    "g1 == @g1 and g2 == @g2 and @pd.isna(jaccard) and not @pd.isna(g2_structure)"
                ).g2_structure.map(lambda x: x.split("_")[-2])
            )
        )
        for g1, g2 in pairs
    ]
).fillna(0)
g2_model_leftovers["g1"], g2_model_leftovers["g2"] = zip(*pairs)
g2_model_leftovers.index = g2_model_leftovers["g1"] + "_-_" + g2_model_leftovers["g2"]
g2_model_leftovers["total_matched"] = common_model_distributions[
    common_model_distributions.columns[:4]
].sum(axis=1)
g2_model_leftovers["priority"] = [
    0 if g1 == single_tissue and g2 == single_tissue else 1
    for g1, g2 in zip(g2_model_leftovers.g1, g2_model_leftovers.g2)
]

swaps = [x for x in g1_model_leftovers.index if x.endswith(single_tissue)]
for s in swaps:
    g1_model_leftovers.loc[s], g2_model_leftovers.loc[s] = (
        g2_model_leftovers.loc[s],
        g1_model_leftovers.loc[s],
    )  # ensuring that the single tissue is always G1


fontsize = 40
ax = (
    g1_model_leftovers.query("g1 == @single_tissue or g2 == @single_tissue")
    .sort_values(["total_matched", "priority"], ascending=[True, False])[
        STRUCTURE_TYPES
    ]
    .plot.barh(
        stacked=True,
        color=["dodgerblue", "orange", "#BE271A", "orchid"],
        legend=False,
        figsize=(7, 9),
    )
)
xupper = 60
ax.set_xlim(0, 60)
ax.set_xticks(range(0, xupper + 1, 10))
ax.set_xticklabels(range(0, xupper + 1, 10), fontsize=fontsize)
ax.set_yticklabels([])
plt.xlabel("Number of structures", fontsize=fontsize)
plt.tight_layout()
plt.savefig(
    f"../graphics/common-models/digestion_{single_tissue}_M1.pdf",
    transparent=True,
    bbox_inches="tight",
)
plt.close()

ax = (
    g2_model_leftovers.query("g1 == @single_tissue or g2 == @single_tissue")
    .sort_values(["total_matched", "priority"], ascending=[True, False])[
        STRUCTURE_TYPES
    ]
    .plot.barh(
        stacked=True,
        color=["dodgerblue", "orange", "darkred", "orchid"],
        legend=False,
        figsize=(7, 9),
    )
)
xupper = 60
ax.set_xlim(0, xupper)
ax.set_xticklabels(range(0, xupper + 1, 10), fontsize=fontsize)
ax.set_yticklabels([])
ax.legend(loc=(0.1425, 0.15), fontsize=fontsize, labelspacing=0.25)
plt.xlabel("Number of structures", fontsize=fontsize)
plt.tight_layout()
plt.savefig(
    f"../graphics/common-models/digestion_{single_tissue}_M2.pdf",
    transparent=True,
    bbox_inches="tight",
)
plt.close()

ax = common_model_distributions.query("g1 == @single_tissue or g2 == @single_tissue").sort_values(['total_matched','priority'], ascending=[True,False]
)[STRUCTURE_TYPES].plot.barh(stacked=True, color=['dodgerblue','orange', 'darkred','orchid'],
                                                     legend=False, figsize=(8.25,9)
                                                    )
xupper = 60
ax.set_xlim(0,xupper)
ax.set_xticks(range(0,xupper+1,10))
ax.set_xticklabels(range(0,xupper+1,10), fontsize=fontsize)
ax.set_yticklabels(common_model_distributions.query("g1 == @single_tissue or g2 == @single_tissue").sort_values(['total_matched','priority'], ascending=[True,False]
).index.map(lambda x:f"{x.split('_-_')[0].split('tissue-')[1][:2]}/{x.split('_-_')[-1].split('tissue-')[1].replace('vermiform_','').replace('large_intestine','il').replace('small_intestine','is')[:2]}"
            if x not in swaps else f"{x.split('_-_')[-1].split('tissue-')[1][:2]}/{x.split('_-_')[0].split('tissue-')[1][:2]}"
           ), fontsize=fontsize)
plt.xlabel("Number of structures", fontsize=fontsize)
plt.tight_layout()
plt.savefig(f"../graphics/common-models/digestion_{single_tissue}_M12.pdf", transparent=True, bbox_inches='tight')
plt.close()
