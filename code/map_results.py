import argparse
import json
import os

import pandas as pd

from custom_utils import get_structures_added, load_json


def set_original_nodes(structure, nodemap):
    stype = structure["structure_type"]
    if stype == "clique":
        structure["nodes"] = sorted(
            map(lambda x: str(nodemap.at[x, "original_id"]), structure["nodes"])
        )
    elif stype == "star":
        structure["hub"] = str(nodemap.at[structure["hub"], "original_id"])
        structure["spokes"] = sorted(
            map(lambda x: str(nodemap.at[x, "original_id"]), structure["spokes"])
        )
    elif stype in ["biclique", "starclique"]:
        structure["left_nodes"] = sorted(
            map(lambda x: str(nodemap.at[x, "original_id"]), structure["left_nodes"])
        )
        structure["right_nodes"] = sorted(
            map(lambda x: str(nodemap.at[x, "original_id"]), structure["right_nodes"])
        )


def map_structures_to_original(g1, path, jsons, nodemaps=None):
    json = load_json(f"{path}/{jsons[g1]}")
    all_structures = (
        get_structures_added(json)
        if nodemaps is not None
        else json[
            "macro_structures"
        ]  # for synthetic data, we don't have the base constraints in the models
    )
    for idx, structure in enumerate(all_structures):
        pos = idx + 1
        structure["position"] = f"{g1}_{structure['structure_type']}_{pos:03}"
    if nodemaps is not None:
        nodemap = pd.read_csv(f"{path}/{nodemaps[g1]}").set_index("julia_id")
        for structure in all_structures:
            set_original_nodes(structure, nodemap)
    return json


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("collections", type=str, nargs="+")
    args = parser.parse_args()
    collections = args.collections

    for dataset in collections:
        path = f"../results/{dataset}"
        result_path = f"../results_mapped/{dataset}"
        if not os.path.exists(result_path):
            os.makedirs(result_path)

        jsons = {
            x.split("_size")[0]: x
            for x in [f for f in os.listdir(path) if f.endswith("json")]
        }

        if dataset != "synthetic":
            nodemaps = {
                x.split("-nodemap")[0]: x
                for x in [f for f in os.listdir(path) if "nodemap" in f]
            }
            assert nodemaps.keys() == jsons.keys()
        else:
            nodemaps = None
        for g1 in sorted(jsons.keys()):
            mapped_model = map_structures_to_original(g1, path, jsons, nodemaps)
            with open(
                f"{result_path}/{g1}.json"
                if dataset != "synthetic"
                else f"{result_path}/{g1}",
                "w",
            ) as f:
                json.dump(mapped_model, f, indent=4, ensure_ascii=False)
