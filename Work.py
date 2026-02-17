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


def trapezoidal_int(x, f):
    cumulative_int = np.zeros_like(x)
    for i in range(1, len(x)):
        dx = x[i] - x[i - 1]
        cumulative_int[i] = (
            cumulative_int[i - 1]
            + dx * (f[i - 1] + f[i]) / 2.0
        )
    return cumulative_int

fig_work, axes_work = plt.subplots(3, 1, figsize=(6, 12), sharex=True)

current_dir = os.getcwd()

for idx, dir_name in enumerate(directories):

    print(f"\nProcessing {dir_name}")
    os.chdir(dir_name)

    force_files = sorted(glob.glob(force_pattern), key=get_index)
    distance_files = sorted(glob.glob(distance_pattern), key=get_index)

    k = 0.008314
    T = 298.15

    w_all = []
    x_all = []

    for force_file, distance_file in zip(force_files, distance_files):

        force = np.loadtxt(force_file, skiprows=17, usecols=1, max_rows=79000)
        dist = np.loadtxt(distance_file, skiprows=17, usecols=1, max_rows=79000)

        order = np.argsort(dist)
        dist = dist[order]
        force = force[order]

        w = trapezoidal_int(dist, force)

        w_all.append(w)
        x_all.append(dist)

        axes_work[idx].plot(
            dist,
            w,
            color=color_dict[dir_name],
            alpha=0.7,
            lw=4
        )

    x = np.mean(np.array(x_all), axis=0)
    w_mean = np.mean(np.array(w_all), axis=0)

    axes_work[idx].set_title(label_dict[dir_name], fontsize=28)
    axes_work[idx].tick_params(axis="both", labelsize=32)


    axes_work[idx].text(
    4.5, 360,
    "Water",
    fontsize=16,
    ha="center",
    va="center"
)

    axes_work[idx].text(
    1.2, 135.0,
    "Fibril",
    fontsize=18,
    ha="center",
    va="center"
)

    axes_work[idx].annotate(
    "",
    xy=(2.1, 40.0),      
    xytext=(1.1, 40.0),  
    arrowprops=dict(
        arrowstyle="->",
        lw=3,
        color="black"
    )
)

    axes_work[idx].text(
        0.02,
        0.93,
        panel_labels[idx],
        transform=axes_work[idx].transAxes,
        fontsize=34,
        fontweight="bold",
        va="top",
        ha="left"
    )


    os.chdir(current_dir)


axes_work[-1].set_xlabel("Distance (nm)", fontsize=32)
axes_work[1].set_ylabel("Work (kJ/mol)", fontsize=32)

fig_work.tight_layout()
fig_work.savefig("all_work.pdf", dpi=300, bbox_inches="tight")

plt.show()
