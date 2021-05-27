import argparse
import os

from custom_utils import write_plots_for_model_json

parser = argparse.ArgumentParser()
parser.add_argument("collection", type=str, help="the name of the collection")
parser.add_argument(
    "--x_granularity_size",
    type=int,
    default=100,
    help="x-axis granularity for size plot",
)
parser.add_argument(
    "--y_granularity_size",
    type=int,
    default=1000,
    help="y-axis granularity for size plot",
)

args = parser.parse_args()
collection = args.collection
path = f"../results/{collection}"
save_directory = f"../graphics/{collection}"

files = sorted(
    [f"../results/{collection}/{f}" for f in os.listdir(path) if f.endswith("json")]
)
base_names = [os.path.split(os.path.splitext(f)[0])[1] for f in files]

if not os.path.exists(save_directory):
    os.makedirs(save_directory)

for file, base_name in zip(files, base_names):
    write_plots_for_model_json(
        file,
        f"{save_directory}/figure-{base_name}",
        args.x_granularity_size,
        args.y_granularity_size,
    )
