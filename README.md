# 21KDD-Momo

## Folders

| Folder            | Description           
| -------------     |-------------| 
| alignments                | Gigi alignments computed *without* node alignments       | 
| alignments-jaccard        | Gigi alignments computed *with* node alignments       | 
| code                      | Code for Beppo, Gigi, NMD computation, plotting, and table generation       | 
| data                      | Preprocessed input data as used in all experiments       | 
| graphics                  | Figures from paper and further visualizations       | 
| network-portrait-divergence   | NPD code and NPD heatmap      | 
| raw-data                  | Raw input data as provided by the stated sources (if different from the preprocessed data) |
| results                   | Beppo results with Julia IDs      | 
| results_mapped            | Beppo results with original IDs   | 

All folders have a separate README.md that describes their contents or usage.


## Setup

The code runs on Debian Debian GNU/Linux 10 (buster) and macOS Catalina; other operating systems were not tested.

We need the following packages when running with Julia version 1.0.3 (this also works with Julia 1.3.0 and slightly older versions of some packages):

```
       "ArgParse" => v"1.1.0"
            "CSV" => v"0.7.4"
     "DataFrames" => v"0.21.4"
 "DataStructures" => v"0.17.19"
        "GraphIO" => v"0.5.0"
      "GraphPlot" => v"0.4.2"
           "JSON" => v"0.21.0"
    "LightGraphs" => v"1.3.3"
   "LineSearches" => v"7.0.1"
          "Optim" => v"0.22.0"
      "StatsBase" => v"0.33.0"
           "YAML" => v"0.4.0"
```

This is the package setup with which the results showcased in the paper were generated.
A better practice for distributing Julia code, which we only learned about later, is documented [here](https://pkgdocs.julialang.org/v1.6/environments/).
For convenience, we also included requirements-minimal Manifest.toml and Project.toml files that should work, e.g., with Julia 1.6, in this folder.

> Julia allows multiple versions of itself to coexist, i.e., we can download Julia 1.6 and set up an alias like `julia16` in the startup script of our preferred shell (e.g., `~/.zshrc`), even if we are normally working with a different version.

To run the code using the more professional setup (assuming Julia 1.6), do the following:
1. Activate the `julia16` REPL
2. Type `]` to go to the package REPL
3. Type `activate .` (this should cause the string between the round brackets to the left of your prompt to change)
4. Type `instantiate` (this will install all packages in the versions required by `Manifest.toml` and allow you to run Julia code)

To run the Python code: 

1. Set up a virtual environment with Python 3.7 and activate it (recommended, not required)
2. Run `pip install -r requirements.txt`

Detailed instructions on how to reproduce the results are given in the README.md of the `code` folder.