import argparse, toml, numpy as np
import structure_clustering
from structure_clustering import element_to_atomic_number

# parsing command-line arguments
parser = argparse.ArgumentParser(
    prog="structure_clustering",
    description="Cluster molecular structures into groups.",
)
parser.add_argument(
    "xyz_file",
    type=str,
    help="path of the multi-xyz-file containing the structures",
)
parser.add_argument(
    "--disconnected",
    help="if you want to include disconnected graphs",
    action="store_true",
)
parser.add_argument(
    "--config",
    type=str,
    help="path of the config TOML file",
    required=False,
)
parser.add_argument(
    "--output",
    type=str,
    help="path of the resulting output file, defaults to <xyz_file>.sc.dat",
    required=False,
)

args = parser.parse_args()

# creating sc machine
sc_machine = structure_clustering.Machine()

# always allow O-H distance up to 2.3 Ang
# if needed otherwise, change via config.TOML
sc_machine.addPairDistance(1, 8, 2.3)

# parsing config file, if provided
if args.config:

    print(f"Loading configuration from {args.config}")
    config_toml = toml.load(args.config)

    # processing connected graphs option
    if "options" in config_toml:
        args.disconnected = not config_toml["options"].get(
            "only_connected_graphs", args.disconnected
        )

    # storing covalent radii within the sc machine, if provided
    if "covalent" in config_toml:
        for element, radius in config_toml["covalent"].items():
            atomic_number = element_to_atomic_number(element)
            print(f"Using covalent radius of {radius} for {element}")
            sc_machine.setCovalentRadius(atomic_number, radius)

    # storing pair distances within the sc machine, if provided
    if "pair" in config_toml:
        for pair, distance in config_toml["pair"].items():
            elements = pair.split("-")
            atomic_numbers = [element_to_atomic_number(elem) for elem in elements]
            print(f"Using pair distance of {distance} for {elements[0]}-{elements[1]}")
            sc_machine.addPairDistance(*atomic_numbers, distance)

if args.disconnected:
    print("Clustering includes disconnected graphs")
else:
    print("Clustering does not include disconnected graphs")
sc_machine.setOnlyConnectedGraphs(not args.disconnected)

# loading structures from multi-xyz file
structures = structure_clustering.import_multi_xyz(args.xyz_file)
print(f"\nUsing {len(structures)} structures from {args.xyz_file}")

# clustering
result = sc_machine.cluster(structures)

# stats
num_structures = len(structures)
num_resulting_structures = sum(map(len, result.clusters)) - len(result.singles)
num_sorted_out = num_structures - num_resulting_structures
percentage_sorted_out = num_sorted_out / num_structures * 100
cluster_sizes_list = list(map(len, result.clusters))
num_connections_list = [struct.numConnections for struct in result.structures]
remaining_structures_idxs = result.singles + [idx for c in result.clusters for idx in c]
num_connections_list_remaining = [
    num_connections_list[i]
    for i in remaining_structures_idxs
    if i < len(num_connections_list)
]


# stats formatter
def format(list):
    list = np.asarray(list)

    avg = np.mean(list)
    median = np.median(list)
    q1 = np.percentile(list, 25)
    q3 = np.percentile(list, 75)

    return f"Avg={avg:.1f} Med={median:.1f} Q1={q1:.1f} Q3={q3:.1f}"


# printing stats
print("Clustering finished", result)
print(f"  {len(result.clusters)} clusters (total {sum(cluster_sizes_list)} structures)")
print(f"  {len(result.singles)} unique single structures")
print(
    f"  {num_sorted_out} ({percentage_sorted_out:.2f}%) structures sorted out ({num_resulting_structures} remaining)"
)
print(f"  cluster size: {format(cluster_sizes_list)}")
print(f"  connections/structure: {format(num_connections_list)} (all {num_structures})")
print(
    f"  connections/structure: {format(num_connections_list_remaining)} (remaining {num_resulting_structures})"
)

# writing output file
if not args.output:
    args.output = args.xyz_file + ".sc.dat"
print(f"Writing output file to {args.output} ...")
result.export(args.output)
msg = "Open https://photophys.github.io/cluster-vis/ to visualize your results"
print(f"\n🚀 \033[1m{msg}\033[0m")
