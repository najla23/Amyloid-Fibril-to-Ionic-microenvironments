#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import os

"""
Author: Najla Hosseini
Email: najla.hosseini@gmail.com| najla.hosseini@ki.se
Affiliation: Department of Oncology-Pathology, Karolinska Institutet,
SE-171 77 Stockholm, Sweden

Description:
This code implements a computational framework combining molecular dynamics,
the Jarzynski estimator, and symmetry-adapted perturbation theory (SAPT)
to analyze amyloid protofibril stability.

Citation:
If you use this code in your research, please cite my article.

"""


directories = [
    "NNQQ_Water",
    "NNQQ_WaterNaCl",
    "NNQQ_WaterKCl",
]

label_dict = {
    "NNQQ_Water": "NNQQ+Water",
    "NNQQ_WaterNaCl": "NNQQ+Water_NaCl",
    "NNQQ_WaterKCl": "NNQQ+Water_KCl",
}

color_dict = {
    "NNQQ_Water": "orange",
    "NNQQ_WaterNaCl": "magenta",
    "NNQQ_WaterKCl": "purple",
}

panel_labels = ["A", "B", "C"]

force_pattern = "slow_pullf*.xvg"
distance_pattern = "slow_pullx*.xvg"

def moving_average(data, window_size):
    return np.convolve(data, np.ones(window_size) / window_size, mode="valid")


def get_index(fname):
    m = re.search(r"(\d+)", fname)
    return int(m.group(1)) if m else -1


fig, axes = plt.subplots(3, 1, figsize=(6, 12), sharex=True)
current_dir = os.getcwd()

for idx, dir_name in enumerate(directories):

    print(f"\nProcessing forces in: {dir_name}")
    os.chdir(dir_name)

    force_files = sorted(glob.glob(force_pattern), key=get_index)
    distance_files = sorted(glob.glob(distance_pattern), key=get_index)

    all_forces = []
    all_distances = []

    for force_file, distance_file in zip(force_files, distance_files):

        force = np.loadtxt(force_file, skiprows=17, usecols=1, max_rows=79000)
        dist = np.loadtxt(distance_file, skiprows=17, usecols=1, max_rows=79000)

        order = np.argsort(dist)
        force = force[order]
        dist = dist[order]

        force_sm = moving_average(force, 3000)
        dist_sm = moving_average(dist, 3000)

        all_forces.append(force_sm)
        all_distances.append(dist_sm)

        axes[idx].plot(
            dist_sm,
            force_sm,
            color=color_dict[dir_name],
            alpha=0.25,
            lw=3
        )

    min_len = min(map(len, all_forces))
    forces_trunc = np.array([f[:min_len] for f in all_forces])
    dists_trunc = np.array([d[:min_len] for d in all_distances])

    mean_force = np.mean(forces_trunc, axis=0)
    mean_dist = np.mean(dists_trunc, axis=0)

    axes[idx].plot(
        mean_dist,
        mean_force,
        color=color_dict[dir_name],
        lw=4,
        label=label_dict[dir_name]
    )

    axes[idx].set_title(label_dict[dir_name], fontsize=28)
    axes[idx].tick_params(axis="both", labelsize=32)

    axes[idx].text(
        0.02,
        0.93,
        panel_labels[idx],
        transform=axes[idx].transAxes,
        fontsize=34,
        fontweight="bold",
        va="top",
        ha="left"
    )

    os.chdir(current_dir)

axes[-1].set_xlabel("Distance (nm)", fontsize=32)
axes[1].set_ylabel("Force (kJ/mol/nm)", fontsize=32)

plt.tight_layout()
plt.savefig("all_force.pdf", dpi=300, bbox_inches="tight")
plt.show()

