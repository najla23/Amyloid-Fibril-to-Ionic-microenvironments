#!/usr/bin/env python3


"""
Author: Najla Hosseini
Email: najla.hosseini@gmail.com| najla.hosseini@ki.se
Affiliation: Department of Oncology-Pathology, Karolinska Institutet,
SE-171 77 Stockholm, Sweden

Description:
This code implements a computational framework combining molecular dynamics,
the Jarzynski estimator,
to analyze amyloid protofibril stability.

Citation:
If you use this code in your research, please cite my article.

"""



import matplotlib.pyplot as plt
import numpy as np
import glob, argparse, csv, hashlib, json, os
from pathlib import Path
import platform, re, shutil,sys, tempfile,warnings



DIRECTORIES = [
    "NNQQ_Water0.0001",
    "NNQQ_WaterNa0.0001",
    "NNQQ_WaterK0.0001",
]

LABEL_DICT = {
    "NNQQ_Water0.0001": "NNQQ+Water",
    "NNQQ_WaterNa0.0001": "NNQQ+Water_NaCl",
    "NNQQ_WaterK0.0001": "NNQQ+Water_KCl",
}

COLOR_DICT = {
    "NNQQ_Water0.0001": "orange",
    "NNQQ_WaterNa0.0001": "magenta",
    "NNQQ_WaterK0.0001": "purple",
}

R_GAS = 0.008314462618 
TEMPERATURE = 298.15   

KT = R_GAS 
BETA = 1.0 / KT
MAX_ROWS = 790000
N_GRID_POINTS = 3000
SMOOTHING_WINDOW = 200
N_BOOTSTRAP = 2000
RANDOM_SEED = 12345
TIME_DECIMALS = 6
MIN_TRAJECTORIES = 2
OUTPUT_PDF = "all_freeenergy_backward_bootstrap.pdf"
OUTPUT_SUMMARY = "jarzynski_bootstrap_summary.txt"
ERROR_DISPLAY = "sem"



def natural_sort_key(path):
    basename = os.path.basename(path)

    return [
        int(token) if token.isdigit() else token.lower()
        for token in re.split(r"(\d+)", basename)
    ]


def extract_file_index(path):
    basename = os.path.basename(path)
    matches = re.findall(r"(\d+)", basename)

    if not matches:
        return None

    return int(matches[-1])


def moving_average(data, window_size):
    data = np.asarray(data, dtype=float)

    if data.ndim != 1:
        raise ValueError("moving_average expects a one-dimensional array.")

    if window_size <= 1:
        return data.copy()

    if window_size > len(data):
        raise ValueError(
            f"Smoothing window ({window_size}) is larger than "
            f"the array length ({len(data)})."
        )

    kernel = np.ones(window_size, dtype=float) / float(window_size)

    return np.convolve(data, kernel, mode="valid")






def smooth_profiles(profiles, window_size):
    profiles = np.asarray(profiles, dtype=float)

    if profiles.ndim != 2:
        raise ValueError("smooth_profiles expects a two-dimensional array.")

    return np.asarray(
        [
            moving_average(profile, window_size)
            for profile in profiles
        ],
        dtype=float
    )



def read_xvg(filename, max_rows=None):
    rows = []

    with open(filename, "r", encoding="utf-8", errors="replace") as handle:

        for line in handle:

            stripped = line.strip()

            if not stripped:
                continue

            if stripped.startswith("#") or stripped.startswith("@"):
                continue

            fields = stripped.split()

            if len(fields) < 2:
                continue

            try:
                time_value = float(fields[0])
                data_value = float(fields[1])

            except ValueError:
                continue

            if not np.isfinite(time_value) or not np.isfinite(data_value):
                continue

            rows.append((time_value, data_value))

            if max_rows is not None and len(rows) >= max_rows:
                break

    if not rows:
        return np.empty((0, 2), dtype=float)

    return np.asarray(rows, dtype=float)


def find_trajectory_pairs(directory):
    force_files = sorted(
        glob.glob(os.path.join(directory, "slow_pullf*.xvg")),
        key=natural_sort_key,
    )

    distance_files = sorted(
        glob.glob(os.path.join(directory, "slow_pullx*.xvg")),
        key=natural_sort_key,
    )

    force_map = {}
    distance_map = {}

    for filename in force_files:
        index = extract_file_index(filename)

        if index is None:
            print(
                f"WARNING: could not determine trajectory index for "
                f"{filename}; skipping it."
            )
            continue

        if index in force_map:
            raise RuntimeError(
                f"Duplicate force-file index {index} in {directory}:\n"
                f"  {force_map[index]}\n"
                f"  {filename}"
            )

        force_map[index] = filename

    for filename in distance_files:
        index = extract_file_index(filename)

        if index is None:
            print(
                f"WARNING: could not determine trajectory index for "
                f"{filename}; skipping it."
            )
            continue

        if index in distance_map:
            raise RuntimeError(
                f"Duplicate distance-file index {index} in {directory}:\n"
                f"  {distance_map[index]}\n"
                f"  {filename}"
            )

        distance_map[index] = filename

    force_indices = set(force_map)
    distance_indices = set(distance_map)

    missing_distance = sorted(force_indices - distance_indices)
    missing_force = sorted(distance_indices - force_indices)

    if missing_distance:
        print(
            f"WARNING: force files without matching distance files "
            f"in {directory}: {missing_distance}"
        )

    if missing_force:
        print(
            f"WARNING: distance files without matching force files "
            f"in {directory}: {missing_force}"
        )

    common_indices = sorted(force_indices & distance_indices)

    pairs = [
        (
            index,
            force_map[index],
            distance_map[index],
        )
        for index in common_indices
    ]

    return pairs


def match_coordinate_force_by_time(
    distance_data,
    force_data,
    decimals=6,
):
    if distance_data.ndim != 2 or distance_data.shape[1] < 2:
        raise ValueError("Invalid distance-data array.")

    if force_data.ndim != 2 or force_data.shape[1] < 2:
        raise ValueError("Invalid force-data array.")

    rounded_distance_times = np.round(
        distance_data[:, 0],
        decimals=decimals,
    )

    rounded_force_times = np.round(
        force_data[:, 0],
        decimals=decimals,
    )

    common_times, distance_indices, force_indices = np.intersect1d(
        rounded_distance_times,
        rounded_force_times,
        assume_unique=False,
        return_indices=True,
    )

    if len(common_times) < 2:
        raise ValueError(
            "Fewer than two common time points were found between "
            "the coordinate and force files."
        )

    distance = distance_data[distance_indices, 1]
    force = force_data[force_indices, 1]

    finite_mask = (
        np.isfinite(common_times)
        & np.isfinite(distance)
        & np.isfinite(force)
    )

    common_times = common_times[finite_mask]
    distance = distance[finite_mask]
    force = force[finite_mask]

    if len(common_times) < 2:
        raise ValueError(
            "Fewer than two finite matched points remained after filtering."
        )

    return common_times, distance, force



def cumulative_trapezoidal_integral(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.ndim != 1 or y.ndim != 1:
        raise ValueError("x and y must be one-dimensional arrays.")

    if len(x) != len(y):
        raise ValueError(
            f"x and y have unequal lengths: {len(x)} and {len(y)}."
        )

    if len(x) < 2:
        raise ValueError("At least two points are needed for integration.")

    result = np.zeros(len(x), dtype=float)

    dx = np.diff(x)
    average_y = 0.5 * (y[:-1] + y[1:])

    result[1:] = np.cumsum(dx * average_y)

    return result


def remove_duplicate_coordinates(distance, force):
    distance = np.asarray(distance, dtype=float)
    force = np.asarray(force, dtype=float)

    unique_distance, inverse_indices = np.unique(
        distance,
        return_inverse=True,
    )

    force_sum = np.zeros(len(unique_distance), dtype=float)
    force_count = np.zeros(len(unique_distance), dtype=float)

    np.add.at(force_sum, inverse_indices, force)
    np.add.at(force_count, inverse_indices, 1.0)

    average_force = force_sum / force_count

    return unique_distance, average_force


def calculate_work_profile(force_file, distance_file):
    force_data = read_xvg(
        force_file,
        max_rows=MAX_ROWS,
    )

    distance_data = read_xvg(
        distance_file,
        max_rows=MAX_ROWS,
    )

    if len(force_data) == 0:
        raise ValueError(
            f"No numerical data were found in {force_file}"
        )

    if len(distance_data) == 0:
        raise ValueError(
            f"No numerical data were found in {distance_file}"
        )

    common_times, distance, force = match_coordinate_force_by_time(
        distance_data,
        force_data,
        decimals=TIME_DECIMALS,
    )

    coordinate_order = np.argsort(distance)

    sorted_distance = distance[coordinate_order]
    sorted_force = force[coordinate_order]

    unique_distance, average_force = remove_duplicate_coordinates(
        sorted_distance,
        sorted_force,
    )

    if len(unique_distance) < 2:
        raise ValueError(
            "Fewer than two unique reaction-coordinate values were found."
        )

    work = cumulative_trapezoidal_integral(
        unique_distance,
        average_force,
    )

    work -= work[0]

    return {
        "time": common_times,
        "distance": unique_distance,
        "force": average_force,
        "work": work,
        "n_force_points": len(force_data),
        "n_distance_points": len(distance_data),
        "n_matched_points": len(common_times),
    }



def interpolate_work_profiles(trajectory_results, n_grid_points):
    common_minimum = max(
        result["distance"].min()
        for result in trajectory_results
    )

    common_maximum = min(
        result["distance"].max()
        for result in trajectory_results
    )

    if common_minimum >= common_maximum:
        ranges = "\n".join(
            (
                f"  trajectory {result['trajectory_index']}: "
                f"{result['distance'].min():.6f} to "
                f"{result['distance'].max():.6f} nm"
            )
            for result in trajectory_results
        )

        raise RuntimeError(
            "The trajectories do not share a common reaction-coordinate "
            f"interval.\n{ranges}"
        )

    common_distance = np.linspace(
        common_minimum,
        common_maximum,
        n_grid_points,
    )

    interpolated_work = []

    for result in trajectory_results:

        work_on_common_grid = np.interp(
            common_distance,
            result["distance"],
            result["work"],
        )

        work_on_common_grid -= work_on_common_grid[0]

        interpolated_work.append(work_on_common_grid)

    interpolated_work = np.asarray(
        interpolated_work,
        dtype=float,
    )

    return (
        common_distance,
        interpolated_work,
        common_minimum,
        common_maximum,
    )



def jarzynski_free_energy(work_profiles, beta):


    work_profiles = np.asarray(work_profiles, dtype=float)

    if work_profiles.ndim != 2:
        raise ValueError(
            "work_profiles must have shape "
            "(number_of_trajectories, number_of_coordinates)."
        )

    exponent = -beta * work_profiles

    maximum_exponent = np.max(
        exponent,
        axis=0,
    )

    shifted_exponentials = np.exp(
        exponent - maximum_exponent
    )

    mean_shifted_exponential = np.mean(
        shifted_exponentials,
        axis=0,
    )

    tiny = np.finfo(float).tiny

    mean_shifted_exponential = np.maximum(
        mean_shifted_exponential,
        tiny,
    )

    log_mean_exponential = (
        maximum_exponent
        + np.log(mean_shifted_exponential)
    )

    free_energy = -log_mean_exponential / beta

    return free_energy



def bootstrap_jarzynski(
    work_profiles,
    beta,
    n_bootstrap,
    random_seed,
):
    work_profiles = np.asarray(work_profiles, dtype=float)

    if work_profiles.ndim != 2:
        raise ValueError("work_profiles must be a two-dimensional array.")

    n_trajectories = work_profiles.shape[0]

    if n_trajectories < 2:
        raise ValueError(
            "At least two trajectories are required for bootstrapping."
        )

    random_generator = np.random.default_rng(random_seed)

    bootstrap_profiles = np.empty(
        (
            n_bootstrap,
            work_profiles.shape[1],
        ),
        dtype=float,
    )

    for bootstrap_index in range(n_bootstrap):

        sampled_indices = random_generator.integers(
            low=0,
            high=n_trajectories,
            size=n_trajectories,
        )

        sampled_work = work_profiles[sampled_indices, :]

        bootstrap_profiles[bootstrap_index, :] = (
            jarzynski_free_energy(
                sampled_work,
                beta,
            )
        )

    return bootstrap_profiles



def reverse_and_shift_profile(distance, free_energy):
    reversed_distance = distance[::-1].copy()
    reversed_free_energy = free_energy[::-1].copy()

    reversed_free_energy -= reversed_free_energy[0]

    return reversed_distance, reversed_free_energy


def reverse_and_shift_bootstrap_profiles(bootstrap_profiles):
    reversed_profiles = bootstrap_profiles[:, ::-1].copy()

    reversed_profiles -= reversed_profiles[:, [0]]

    return reversed_profiles



def process_system(directory):
    print("\n" + "=" * 80)
    print(f"PROCESSING SYSTEM: {directory}")
    print("=" * 80)

    if not os.path.isdir(directory):
        raise FileNotFoundError(
            f"Directory does not exist: {directory}"
        )

    trajectory_pairs = find_trajectory_pairs(directory)

    if not trajectory_pairs:
        raise RuntimeError(
            f"No matching slow_pullf*.xvg and slow_pullx*.xvg "
            f"pairs were found in {directory}."
        )

    trajectory_results = []

    for trajectory_index, force_file, distance_file in trajectory_pairs:

        print(
            f"\nTrajectory {trajectory_index}\n"
            f"  Force file:    {force_file}\n"
            f"  Distance file: {distance_file}"
        )

        try:
            result = calculate_work_profile(
                force_file,
                distance_file,
            )

        except (
            FileNotFoundError,
            OSError,
            ValueError,
            FloatingPointError,
        ) as error:

            print(
                f"  WARNING: trajectory {trajectory_index} "
                f"was skipped.\n"
                f"  Reason: {error}"
            )
            continue

        result["trajectory_index"] = trajectory_index
        result["force_file"] = force_file
        result["distance_file"] = distance_file

        trajectory_results.append(result)


    if len(trajectory_results) < MIN_TRAJECTORIES:
        raise RuntimeError(
            f"Only {len(trajectory_results)} valid trajectories were "
            f"found for {directory}. At least {MIN_TRAJECTORIES} are "
            f"required."
        )

    (
        common_distance,
        work_profiles,
        common_minimum,
        common_maximum,
    ) = interpolate_work_profiles(
        trajectory_results,
        N_GRID_POINTS,
    )

    print(
        f"\nValid trajectories used: {len(trajectory_results)}"
    )

    print(
        f"Common coordinate range: "
        f"{common_minimum:.6f} to "
        f"{common_maximum:.6f} nm"
    )


    free_energy = jarzynski_free_energy(
        work_profiles,
        BETA,
    )


    bootstrap_profiles = bootstrap_jarzynski(
        work_profiles,
        BETA,
        N_BOOTSTRAP,
        RANDOM_SEED,
    )


    reversed_distance, reversed_free_energy = (
        reverse_and_shift_profile(
            common_distance,
            free_energy,
        )
    )

    reversed_bootstrap = (
        reverse_and_shift_bootstrap_profiles(
            bootstrap_profiles
        )
    )


    if SMOOTHING_WINDOW > 1:

        smoothed_distance = moving_average(
            reversed_distance,
            SMOOTHING_WINDOW,
        )

        smoothed_free_energy = moving_average(
            reversed_free_energy,
            SMOOTHING_WINDOW,
        )

        smoothed_bootstrap = smooth_profiles(
            reversed_bootstrap,
            SMOOTHING_WINDOW,
        )

    else:
        smoothed_distance = reversed_distance.copy()
        smoothed_free_energy = reversed_free_energy.copy()
        smoothed_bootstrap = reversed_bootstrap.copy()


    bootstrap_mean = np.mean(
        smoothed_bootstrap,
        axis=0,
    )

    bootstrap_sem = np.std(
        smoothed_bootstrap,
        axis=0,
        ddof=1,
    )

    maximum_bootstrap_difference = np.max(
        np.abs(
            smoothed_free_energy
            - bootstrap_mean
        )
    )

    print(
        f"Maximum difference between the full-sample profile and "
        f"the bootstrap mean: "
        f"{maximum_bootstrap_difference:.6f} kJ/mol"
    )

    final_free_energy = smoothed_free_energy[-1]
    final_bootstrap_sem = bootstrap_sem[-1]

    print(
        f"Final relative free energy: "
        f"{final_free_energy:.6f} kJ/mol"
    )

    print(
        f"Bootstrap standard error at final point: "
        f"{final_bootstrap_sem:.6f} kJ/mol"
    )


    return {
        "directory": directory,
        "trajectory_results": trajectory_results,
        "common_distance": common_distance,
        "work_profiles": work_profiles,
        "n_trajectories": len(trajectory_results),
        "distance": smoothed_distance,
        "free_energy": smoothed_free_energy,
        "bootstrap_sem": bootstrap_sem,
        "final_free_energy": final_free_energy,
        "final_bootstrap_sem": final_bootstrap_sem,
        "common_minimum": common_minimum,
        "common_maximum": common_maximum,
    }



def main():
    if ERROR_DISPLAY not in {"sem", "ci"}:
        raise ValueError(
            'ERROR_DISPLAY must be either "sem" or "ci".'
        )

    if N_BOOTSTRAP < 2:
        raise ValueError(
            "N_BOOTSTRAP must be at least 2."
        )

    if N_GRID_POINTS < 2:
        raise ValueError(
            "N_GRID_POINTS must be at least 2."
        )

    if SMOOTHING_WINDOW < 1:
        raise ValueError(
            "SMOOTHING_WINDOW must be at least 1."
        )

    if SMOOTHING_WINDOW > N_GRID_POINTS:
        raise ValueError(
            "SMOOTHING_WINDOW cannot be greater than "
            "N_GRID_POINTS."
        )

    system_results = {}

    plt.figure(figsize=(8, 6))

    for directory in DIRECTORIES:

        try:
            result = process_system(directory)

        except (
            FileNotFoundError,
            RuntimeError,
            ValueError,
            OSError,
        ) as error:

            warnings.warn(
                f"System {directory} was not plotted: {error}"
            )
            continue

        system_results[directory] = result

        distance = result["distance"]
        free_energy = result["free_energy"]

        if ERROR_DISPLAY == "sem":

            lower_bound = (
                free_energy
                - result["bootstrap_sem"]
            )

            upper_bound = (
                free_energy
                + result["bootstrap_sem"]
            )


        plt.fill_between(
            distance,
            lower_bound,
            upper_bound,
            color=COLOR_DICT[directory],
            alpha=0.25,
            linewidth=0,
        )

        plt.plot(
            distance,
            free_energy,
            linewidth=4,
            color=COLOR_DICT[directory],
            label=LABEL_DICT[directory],
        )

    if not system_results:
        raise RuntimeError(
            "No valid systems were processed; no figure was created."
        )


    plt.xlabel(
        "Distance (nm)",
        fontsize=32,
    )

    plt.ylabel(
        "PMF (kJ/mol)",
        fontsize=32,
    )

    plt.tick_params(
        axis="both",
        labelsize=32,
    )

    plt.legend(
        frameon=False,
        fontsize=20,
    )

    plt.text(
        4.5,
        -50,
        "Solution",
        fontsize=16,
        horizontalalignment="center",
        verticalalignment="center",
    )

    plt.text(
        1.09,
        -40.0,
        "Fibril",
        fontsize=18,
        horizontalalignment="center",
        verticalalignment="center",
    )

    plt.annotate(
        "",
        xy=(3.0, 0.0),
        xytext=(4.5, 0.0),
        arrowprops={
            "arrowstyle": "->",
            "linewidth": 3,
            "color": "black",
        },
    )

    plt.tight_layout()

    plt.savefig(
        OUTPUT_PDF,
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()


    with open(
        OUTPUT_SUMMARY,
        "w",
        encoding="utf-8",
    ) as summary_file:

        summary_file.write(
            "Jarzynski free-energy analysis with "
            "trajectory-level bootstrap uncertainty\n"
        )

        summary_file.write(
            "=" * 78 + "\n\n"
        )

        summary_file.write(
            f"Temperature: {TEMPERATURE:.2f} K\n"
        )

        summary_file.write(
            f"kT: {KT:.8f} kJ/mol\n"
        )

        summary_file.write(
            f"Bootstrap samples: {N_BOOTSTRAP}\n"
        )

        summary_file.write(
            f"Common grid points: {N_GRID_POINTS}\n"
        )

        summary_file.write(
            f"Smoothing window: {SMOOTHING_WINDOW}\n"
        )

        summary_file.write(
            f"Displayed uncertainty: {ERROR_DISPLAY}\n\n"
        )

        for directory in DIRECTORIES:

            if directory not in system_results:
                continue

            result = system_results[directory]

            summary_file.write(
                f"System: {LABEL_DICT[directory]}\n"
            )

            summary_file.write(
                f"Directory: {directory}\n"
            )

            summary_file.write(
                f"Valid trajectories: "
                f"{result['n_trajectories']}\n"
            )

            summary_file.write(
                f"Common coordinate range: "
                f"{result['common_minimum']:.8f} to "
                f"{result['common_maximum']:.8f} nm\n"
            )

            summary_file.write(
                f"Final relative free energy: "
                f"{result['final_free_energy']:.8f} kJ/mol\n"
            )

            summary_file.write(
                f"Bootstrap standard error: "
                f"{result['final_bootstrap_sem']:.8f} kJ/mol\n"
            )


            summary_file.write("\n")

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"Figure written to:  {OUTPUT_PDF}")
    print(f"Summary written to: {OUTPUT_SUMMARY}")



REPRODUCIBILITY_ATOL = 1.0e-10  # absolute comparison
SETTING_NAMES = (
    "DIRECTORIES", "LABEL_DICT", "COLOR_DICT", "R_GAS", "TEMPERATURE",
    "KT", "BETA", "MAX_ROWS", "N_GRID_POINTS", "SMOOTHING_WINDOW",
    "N_BOOTSTRAP", "RANDOM_SEED", "TIME_DECIMALS", "MIN_TRAJECTORIES",
    "ERROR_DISPLAY",
)
PROFILE_KEYS = (
    "distance", "free_energy", "bootstrap_sem",
)
SUMMARY_KEYS = (
    "n_trajectories", "final_free_energy", "final_bootstrap_sem",
)


def _write_json(path, value):
    with Path(path).open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def _read_json(path):
    with Path(path).open(encoding="utf-8") as handle:
        return json.load(handle)


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_matrix(path, values, headers, fmt="%.17g"):
    with Path(path).open("w", encoding="utf-8") as handle:
        handle.write("\n".join(headers) + "\n")
        np.savetxt(handle, values, fmt=fmt)


def _read_matrix(path):
    with Path(path).open(encoding="utf-8") as handle:
        lines = (line for line in handle
                 if line.strip() and not line.lstrip().startswith(("#", "@")))
        return np.loadtxt(lines, ndmin=2)


def export_processed_system(result, destination):
    trajectories = result["trajectory_results"]
    if len(trajectories) != 15:
        raise ValueError(
            f"Export requires 15 valid trajectory pairs for indices 0--14; "
            f"found {len(trajectories)}. No trajectories are dropped or duplicated."
        )
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    mappings, file_hashes = [], {}

    for new_index, trajectory in enumerate(trajectories, start=0):
        coordinate = trajectory["distance"]
        force = trajectory["force"]
        sample_index = np.arange(len(coordinate), dtype=np.int64)
        xname = f"slow_pullx{new_index}.xvg"
        fname = f"slow_pullf{new_index}.xvg"
        common_headers = [
            f"# Original trajectory index: {trajectory['trajectory_index']}",
            '@ xaxis label "Sample index"',
            "@ TYPE xy",
        ]
        _write_matrix(
            destination / xname, np.column_stack((sample_index, coordinate)),
            common_headers + ['@ yaxis label "Position (nm)"'],
            fmt=["%d", "%.17g"],
        )
        _write_matrix(
            destination / fname, np.column_stack((sample_index, force)),
            common_headers + [r'@ yaxis label "Force (kJ mol\S-1\N nm\S-1\N)"'],
            fmt=["%d", "%.17g"],
        )
        file_hashes[xname] = _sha256(destination / xname)
        file_hashes[fname] = _sha256(destination / fname)
        mappings.append({
            "new_index": new_index,
            "original_index": trajectory["trajectory_index"],
            "original_pullx": str(Path(trajectory["distance_file"]).resolve()),
            "original_pullf": str(Path(trajectory["force_file"]).resolve()),
            "exported_pullx": xname, "exported_pullf": fname,
            "original_read_x_rows": trajectory["n_distance_points"],
            "original_read_f_rows": trajectory["n_force_points"],
            "matched_time_points": trajectory["n_matched_points"],
            "exported_unique_positions": len(coordinate),
            "original_first_matched_time": float(trajectory["time"][0]),
            "original_last_matched_time": float(trajectory["time"][-1]),
        })

    with (destination / "mapping.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(mappings[0]))
        writer.writeheader()
        writer.writerows(mappings)

    _write_matrix(
        destination / "reference_profile.xvg",
        np.column_stack([result[key] for key in PROFILE_KEYS]),
        ["# Columns: distance_nm free_energy_kJmol bootstrap_SE",
         '@ xaxis label "Distance (nm)"', '@ yaxis label "Free energy (kJ/mol)"'],
    )
    _write_matrix(
        destination / "reference_work_profiles.xvg",
        np.column_stack((result["common_distance"], result["work_profiles"].T)),
        ["# Columns: distance_nm followed by W_0 ... W_N in kJ/mol.",
         '@ xaxis label "Distance (nm)"', '@ yaxis label "Integrated force (kJ/mol)"'],
    )
    _write_json(destination / "reference_summary.json",
                {key: result[key] for key in SUMMARY_KEYS})
    _write_json(destination / "input_hashes.json", file_hashes)


def verify_processed_system(result, directory):
    directory = Path(directory)
    expected_hashes = _read_json(directory / "input_hashes.json")
    present = {p.name for pattern in ("slow_pullx*.xvg", "slow_pullf*.xvg")
               for p in directory.glob(pattern)}
    if present != set(expected_hashes):
        raise RuntimeError(f"Input file set changed in {directory}")
    for name, expected_hash in expected_hashes.items():
        if _sha256(directory / name) != expected_hash:
            raise RuntimeError(f"Exported input was changed: {directory / name}")

    profile = _read_matrix(directory / "reference_profile.xvg")
    work = _read_matrix(directory / "reference_work_profiles.xvg")
    summary = _read_json(directory / "reference_summary.json")
    report = {"passed": True, "absolute_tolerance": REPRODUCIBILITY_ATOL,
              "maximum_absolute_differences": {}, "numerically_exact": {}}
    comparisons = [(key, result[key], profile[:, i])
                   for i, key in enumerate(PROFILE_KEYS)]
    comparisons += [
        ("common_distance", result["common_distance"], work[:, 0]),
        ("work_profiles", result["work_profiles"], work[:, 1:].T),
    ]
    comparisons += [(key, result[key], summary[key]) for key in SUMMARY_KEYS]
    for key, actual, expected in comparisons:
        actual, expected = np.asarray(actual), np.asarray(expected)
        if actual.shape != expected.shape:
            raise RuntimeError(f"Replay shape mismatch for {key}: {actual.shape} != {expected.shape}")
        difference = float(np.max(np.abs(actual - expected)))
        if not np.isfinite(difference) or difference > REPRODUCIBILITY_ATOL:
            raise RuntimeError(f"Replay mismatch for {key}: maximum difference {difference}")
        report["maximum_absolute_differences"][key] = difference
        report["numerically_exact"][key] = bool(np.array_equal(actual, expected))
    return report


def _release_intermediate_arrays(result):
    for key in ("trajectory_results", "common_distance", "work_profiles"):
        result.pop(key, None)
    return result


def run_export_cli(argv=None):
    parser = argparse.ArgumentParser(
        description="Run the original profile calculation and export reproducible XVG inputs.",
        epilog="Requires NumPy and Matplotlib. Exported first columns are sample indices, not MD times.",
    )
    parser.add_argument("--input-root", type=Path, default=Path("."),
                        help="Parent of the system directories in DIRECTORIES")
    parser.add_argument("--export-root", type=Path, default=Path("NNQQ_Solution_0.0001"),
                        help="New directory for all three systems; must not already exist")
    parser.add_argument("--replay", type=Path,
                        help="Reanalyse an exported bundle using its saved settings")
    parser.add_argument("--show", action="store_true", help="Also open the plot window")
    args = parser.parse_args(argv)

    original_process = process_system
    original_show = plt.show
    original_pdf, original_summary = OUTPUT_PDF, OUTPUT_SUMMARY
    saved_globals = {name: globals()[name] for name in SETTING_NAMES}
    if not args.show:
        plt.show = lambda: None

    try:
        if args.replay is not None:
            root = args.replay.resolve()
            metadata = _read_json(root / "analysis_settings.json")
            if metadata.get("format_version") != 1:
                raise ValueError("Unsupported processed-input format")
            for name in SETTING_NAMES:
                globals()[name] = metadata["settings"][name]
            checked = []
            reports = {}

            def process_replay(directory):
                folder = root / metadata["system_folders"][directory]
                result = original_process(str(folder))
                reports[directory] = verify_processed_system(result, folder)
                checked.append(directory)
                print(f"REPLAY VERIFIED: {directory}")
                return _release_intermediate_arrays(result)

            globals()["process_system"] = process_replay
            globals()["OUTPUT_PDF"] = str(root / ("replayed_" + Path(original_pdf).name))
            globals()["OUTPUT_SUMMARY"] = str(root / ("replayed_" + Path(original_summary).name))
            main()
            if checked != DIRECTORIES:
                raise RuntimeError("Replay did not verify every configured system")
            _write_json(root / "replay_verification.json", reports)
            print(f"All {len(checked)} systems reproduce the stored profiles and uncertainty bands.")
            return

        source_root = args.input_root.resolve()
        destination = args.export_root.resolve()
        if destination.exists():
            raise FileExistsError(f"Export directory already exists: {destination}")
        names = {directory: Path(directory).name for directory in DIRECTORIES}
        if len(set(names.values())) != len(DIRECTORIES) or any(n in ("", ".", "..") for n in names.values()):
            raise ValueError("Each system must have a unique, nonempty folder name")
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".pmf_export_", dir=destination.parent) as temp:
            staging = Path(temp) / "bundle"
            staging.mkdir()
            completed = []

            def process_and_export(directory):
                source = source_root / directory
                result = original_process(str(source))
                exported = staging / names[directory]
                export_processed_system(result, exported)
                # Run the same original calculation again, on the written XVGs.
                replayed = original_process(str(exported))
                report = verify_processed_system(replayed, exported)
                _write_json(exported / "verification.json", report)
                completed.append(directory)
                print(f"EXPORT VERIFIED: {directory}")
                return _release_intermediate_arrays(result)

            globals()["process_system"] = process_and_export
            globals()["OUTPUT_PDF"] = str(staging / Path(original_pdf).name)
            globals()["OUTPUT_SUMMARY"] = str(staging / Path(original_summary).name)
            main()
            if completed != DIRECTORIES:
                raise RuntimeError("Not every configured system passed; no export bundle was published")
            _write_json(staging / "analysis_settings.json", {
                "format_version": 1,
                "settings": saved_globals,
                "system_folders": names,
                "python_version": platform.python_version(),
                "numpy_version": np.__version__,
                "analysis_script_sha256": _sha256(__file__),
                "input_kind": "sorted unique positions and duplicate-averaged forces",
                "first_column": "sample index, not simulation time",
                "extra_position_cutoff": None,
            })
            shutil.copyfile(__file__, staging / "analysis_used.py")
            (staging / "README.txt").write_text(
                "REPRODUCIBLE PROCESSED XVG INPUTS\n\n"
                "Each system contains paired slow_pullxN.xvg / slow_pullfN.xvg files.\n"
                "Temperature and all estimator/smoothing/bootstrap settings are restored on replay.\n\n"
                "From inside this directory run:\n"
                "    python3 analysis_used.py --replay .\n\n",
                encoding="utf-8",
            )
            if destination.exists():
                raise FileExistsError(f"Export directory appeared during processing: {destination}")
            staging.rename(destination)
        print(f"Verified processed XVG directories saved to: {destination}")
    finally:
        globals()["process_system"] = original_process
        globals()["OUTPUT_PDF"], globals()["OUTPUT_SUMMARY"] = original_pdf, original_summary
        for name, value in saved_globals.items():
            globals()[name] = value
        plt.show = original_show
        plt.close("all")


if __name__ == "__main__":
    try:
        run_export_cli()
    except (OSError, ValueError, RuntimeError, KeyError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
