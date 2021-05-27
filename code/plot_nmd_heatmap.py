import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_style("whitegrid")


def plot_heatmap(df_to_plot, select, collection, save_path):
    fig, ax = plt.subplots(figsize=(20, 20))
    pivoted_df = df_to_plot.pivot(index="g1", columns="g2", values="nmd").fillna(0)
    pivoted_df += pivoted_df.T
    cbar = (
        True
        if select == "all" or collection not in ["law", "as-data", "arxiv-temporal"]
        else False
    )
    mask = np.triu(np.ones_like(pivoted_df, dtype=bool), k=1)
    fig = sns.heatmap(
        pivoted_df.mask(mask)
        if not cbar
        else pivoted_df[len(pivoted_df) // 2 :].T[: len(pivoted_df) // 2],
        square=True,
        cbar=cbar,
        cmap="viridis",
        vmin=0.0,
        cbar_kws={"shrink": 0.7},
        xticklabels="",
        yticklabels="",
    )
    if cbar:
        fig.collections[0].colorbar.ax.tick_params(labelsize=80)
        fig.collections[0].colorbar.ax.set_ylabel("", fontsize=80)
    ax.set_facecolor("white")
    plt.xlabel("")
    plt.ylabel("")
    plt.tight_layout()

    plt.savefig(
        f"{save_path}/nmd-{select}.pdf",
        transparent=True,
        bbox_inches="tight",
    )
    plt.close()


def get_synthetic_clustermap_label(x):
    return f'{x.replace("synthetic_", "").split("_n")[0].replace("_", "/")}/i={(int(int(x.split("_")[-3][1:]) / 100000)):02}'


def plot_square_heatmap(df_to_plot, save_path):
    sns.set_style()
    plt.rcParams["xtick.major.pad"] = 20
    plt.rcParams["ytick.minor.pad"] = 20
    plt.rcParams["ytick.major.pad"] = 20
    fontsize = 90
    fig, ax = plt.subplots(figsize=(20, 20))
    pivoted_df = df_to_plot.pivot(index="g1", columns="g2", values="nmd").fillna(0)
    pivoted_df += pivoted_df.T
    cbar = False
    fig = sns.heatmap(
        pivoted_df,
        square=True,
        cbar=cbar,
        cmap="viridis",
        vmin=0.0,
        cbar_kws={"shrink": 0.76, "fraction": 0.1},
        xticklabels=1,
        yticklabels=1,
        vmax=1.0,
    )
    xlabels = list(
        map(lambda x: x.get_text(), fig.collections[0].axes.get_xticklabels())
    )
    ylabels = list(
        map(lambda x: x.get_text(), fig.collections[0].axes.get_yticklabels())
    )
    fig.collections[0].axes.set_xticklabels(xlabels, fontsize=fontsize)
    fig.collections[0].axes.set_yticklabels(ylabels, fontsize=fontsize)
    plt.yticks(
        np.arange(len(fig.get_yticklabels())) + 0.5,
        fig.get_yticklabels(),
        rotation=90,
        va="center",
    )

    ax.set_facecolor("white")
    plt.xlabel("")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(f"{save_path}/nmd-square.pdf", transparent=True, bbox_inches="tight")
    plt.close()


def get_row_color(x):
    if "cl100" in x:
        return "blue"  # "#0047ab"
    elif "st100" in x:
        return "#FFED00"
    elif "bc100" in x:
        return "#FF0000"
    elif "sc100" in x:
        return "#da70d6"  # "#FF00AB"
    elif "cl050" in x and "st050" in x:
        return "#809A56"
    elif "cl050" in x and "bc050" in x:
        return "#802456"
    elif "cl050" in x and "sc050" in x:
        return "#8024AB"
    elif "st050" in x and "bc050" in x:
        return "#FF7700"
    elif "st050" in x and "sc050" in x:
        return "#FF7756"
    elif "bc050" in x and "sc050" in x:
        return "#FF0056"
    elif "cl033" in x and "st033" in x and "bc033" in x:
        return "#AA6739"
    elif "cl033" in x and "st033" in x and "sc033" in x:
        return "#AA6772"
    elif "cl033" in x and "bc033" in x and "sc033" in x:
        return "#AA1872"
    elif "st033" in x and "bc033" in x and "sc033" in x:
        return "#FF4F39"
    elif "cl025" in x and "st025" in x and "bc025" in x and "sc025" in x:
        return "#BF4D56"
    else:
        raise Exception(x)


def plot_colored_clustermap(df_to_plot, save_path):
    sns.set(font_scale=0.25)
    pivoted_df = df_to_plot.pivot(index="g1", columns="g2", values="nmd").fillna(0)
    pivoted_df += pivoted_df.T
    if collection == "synthetic":
        colcolors = df_to_plot.g2.map(get_row_color).values
        pivoted_df = df_to_plot.pivot(index="g1", columns="g2", values="nmd").fillna(0)
        pivoted_df += pivoted_df.T
        fig = sns.clustermap(
            pivoted_df,
            method="single",
            cmap="Reds",
            xticklabels=pivoted_df.index.map(get_synthetic_clustermap_label),
            yticklabels=pivoted_df.columns.map(get_synthetic_clustermap_label),
            row_colors=colcolors,
            col_colors=colcolors,
            cbar_kws=dict(fraction=0.1),
        )
        fig.ax_cbar.set_ylabel("NMD", fontsize=16)
        fig.ax_cbar.tick_params(labelsize=16)
        fig.ax_heatmap.set_xlabel("$G_1$", fontsize=16)
        fig.ax_heatmap.set_ylabel("$G_2$", fontsize=16)
        for a in fig.ax_row_dendrogram.collections:
            a.set_linewidth(2)
        for a in fig.ax_col_dendrogram.collections:
            a.set_linewidth(2)
        fig.savefig(f"{save_path}/nmd-clustermap-single.pdf")
        plt.close()


def get_selectors(collection):
    if collection == "law":
        return ["all", "de", "us"]
    if collection == "arxiv-temporal":
        return ["all", "cs.LG", "cs.SI"]
    if collection == "as-data":
        return ["all", "oregon1", "oregon2"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "collection",
        type=str,
        choices=[
            "as-data",
            "arxiv-temporal",
            "bio",
            "law",
            "portrait-divergence",
            "synthetic",
        ],
    )
    parser.add_argument(
        "--variant", "-v", choices=["overlap", "jaccard"], default="overlap"
    )
    args = parser.parse_args()
    collection = args.collection
    variant = args.variant
    save_path = f"../graphics/heatmaps/{collection}"
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    df = pd.read_csv(f"../results_mapped/{collection}/nmd-{variant}.csv")

    if collection == "portrait-divergence":
        df.g1 = [x.replace("dev-dev_", "") for x in df.g1]
        df.g2 = [x.replace("dev-dev_", "") for x in df.g2]
        plot_square_heatmap(df, save_path)
    elif collection == "synthetic":
        df.g1 = [x[:-5].replace("synthetic-", "") for x in df.g1]
        df.g2 = [x[:-5].replace("synthetic-", "") for x in df.g2]
        plot_colored_clustermap(df, save_path)
    else:
        selectors = get_selectors(collection)
        for select in selectors:
            df_to_plot = (
                df.query(
                    f"g1.str.startswith('{select}') and g2.str.startswith('{select}')"
                )
                if select != "all"
                else df
            )
            df_to_plot = df_to_plot.copy()
            df_to_plot.g1 = [
                x.split("_size")[0].replace("tissue-", "")
                + x.split("max-100.0")[-1][:-5]
                for x in df_to_plot.g1
            ]
            df_to_plot.g2 = [
                x.split("_size")[0].replace("tissue-", "")
                + x.split("max-100.0")[-1][:-5]
                for x in df_to_plot.g2
            ]
            plot_heatmap(df_to_plot, select, collection, save_path)
