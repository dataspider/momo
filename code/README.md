# code

## Input data preparation

We provide all preprocessed data in the `data` folder.
To produce this data from the raw data (or generate it from scratch in the case of random graphs), for *clg* and *csi* (*arxiv-temporal*), *bio*, *lde* and *lus* (*law*), as well as *rba* and *rer*, we proceed as follows.

To prepare the *arxiv-temporal* data,
- download the November 2020 snapshot of arXiv from Kaggle and store it in the arxiv-raw subfolder of the raw-data directory as `arxiv-metadata-oai-snapshot.json` (we do not provide this file directly because it is publicly available and rather large), then
- run `python prepare_arxiv_data.py raw-data/arxiv-raw data/arxiv-temporal`
    (this will also create category files, bipartite graphs, and author projections of the bipartite graphs for all arXiv categories as of 2020-11-01).

To prepare the *bio* data:

`python prepare_bio_data.py raw-data/bio-raw/PPT-Ohmnet_tissues-combined.edgelist data/bio`

To prepare the *law* data:

`python prepare_law_data.py raw-data/law-raw data/law`

To create the *random* graphs:

`python generate_random_graphs.py`

Finally, to generate the toy graphs from Figure 1 and the synthetic models for Figure 12:

`python generate_planted.py` and `python generate_synthetic_models.py`

## Beppo

### Computing individual models in Julia
For all experiments except the comparison of NMD and NPD, to compute models for all graphs in `folder`:

With parallelization:

`julia -p {number of parallel processes} run_models.jl ../data/{folder} --stop_after_n_structures 100. --stop_threshold 300.`

Without parallelization:

`julia run_models.jl ../data/{folder} -s true --stop_after_n_structures 100. --stop_threshold 300.`

For the comparison of NMD and NPD:

`julia run_models.jl ../data/portrait-divergence -s true --stop_after_n_structures 100. --stop_threshold 300. --size_threshold 3.`

This creates a `folder` directory in the `results` directory (and also create that directory, if it is not already present) and write three files per graph to this directory:
- A JSON file with the individual model
- A log file documenting how the model is computed
- A CSV file mapping Julia node IDs to the original node IDs

To overwrite results that are already present (e.g., if you have unpacked the results folder in this repository), add an `-o true`.

### Mapping node IDs in the Julia models back to original IDs

To map results stored in `results/{folder}` to their original IDs, run:

`python map_results.py {folder}`

This creates a `folder` directory in the `results_mapped` directory (and also create that directory, if it is not already present) and write one JSON file with the mapped individual model to that directory per individual model in the `results/{folder}` directory.

## Gigi

### Precomputing Jaccard similarities between structures in individual models

This step is needed if the Gigi alignment is to be computed without node alignment and we want to consider overlap between structures *within* graphs.

To create Jaccard similarity matrices for all graphs within `results/{folder}`:

`python plot_models.py {folder}`

This also creates all plots for each graph in the `graphics` folder (e.g., node overlap trees).

### Computing alignments

To compute Gigi alignments based on the Jaccard overlap of structures *within* the graphs in `results_mapped/{folder}`:

`python compute_structure_alignment.py {folder}`

If the node IDs of the graphs are (partially) aligned, such that nodes with identical IDs should be treated as identical and nodes with different IDs should be treated as different, add an `-a` flag to compute Gigi alignments based on the Jaccard overlap of structures *between* the graphs.

To parallelize over x processes, add `-p x`.
We do this for the `bio` collection.

To guarantee that structures do not overlap on nodes (i.e., forgo the product graph construction and skip directly to the greedy matching part of MaximalGreedy if no node alignment is present), add an `-o` flag. 
We do this for the `synthetic` collection.

To allow parallelization, the alignments are stored in several files (one per graph).
The storage folder is `alignments` if the `-a` flag is not set, and `alignments_jaccard` otherwise. 

### Merging alignments

To merge alignments stored in `alignments/{folder}`:

`python merge_alignments.py {folder}`

To merge alignments stored in `alignments_jaccard/{folder}`:

`python merge_alignments.py {folder} -v jaccard`

This stores a file `structure-alignments-{overlap|jaccard}.csv` with all alignments in the source path.

## NMD computation

To compute NMDs based on merged alignments in `alignments/{folder}`:

`python compute_graph_distance {folder}`

To compute NMDs based on merged alignments in `alignments_jaccard/{folder}`:

`python compute_graph_distance {folder} -v jaccard`

This stores a file `nmd-{overlap|jaccard}.csv` in `results_mapped/{folder}`.


## Figures

The given directives assume that all necessary input results from Beppo, Gigi, or NMD computation are present (e.g., as provided in the other folders of this repository).

Most of the plotting code was originally developed in Jupyter notebook and then ported to Python scripts to facilitate reproduction.

| Figure            | Description    | Command|       
| -------------     |-------------| ---|
|1|Toy graphs|`python plot_planted.py`|
|2|Structure types|- (created manually)|
|3|Beppo performance|`python plot_performance.py`|
|4|*csi* *2020* example star|`python plot_arxiv_star.py`|
|5|Individual node overlap trees|`python plot_models.py bio` |
|6|Individual and common model treemaps|`python plot_treemaps.py`|
|7|Common model and transformation stacked bar charts|`python plot_multilateral_comparison.py`|
|8|NMD distribution per collection|`python plot_nmd_distributions.py`|
|9|NMD heatmaps for *clg* and *csi*|`python plot_nmd_heatmap.py arxiv-temporal`|
|10|Distribution of *n* and *m* per collection|`python plot_dataset_summary_statistics.py`|
|11|Common node overlap trees for *lde* 2018 and 2019|`python plot_common_overlap_trees.py`|
|12|NMD clustermap for *synthetic*|`python plot_nmd_heatmap.py synthetic`|
|13|NMD vs. NPD comparison|`plot_nmd_heatmap.py portrait-divergence` (for right part); `python portrait-data.py` in the `network-portrait-divergence` folder (for left part)|