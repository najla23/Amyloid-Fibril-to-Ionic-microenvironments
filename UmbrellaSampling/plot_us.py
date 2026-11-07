#!/usr/bin/env python3
"""
Author: Najla Hosseini
Email: najla.hosseini@gmail.com| najla.hosseini@ki.se
Affiliation: Department of Oncology-Pathology, Karolinska Institutet,
SE-171 77 Stockholm, Sweden

Panel A: umbrella histograms. 
Panel B: US PMF for NNQQ1-water system"""

import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np


HISTOGRAM_FILE = "histo23.xvg"
PROFILE_FILE = "profile23.xvg"
OUTPUT_BASENAME = "US_histograms_PMF"


def read_xvg(path, usecols=None):
    try:
        data = np.loadtxt(
            path, comments=("#", "@"), usecols=usecols, ndmin=2,
        )
    except (OSError, ValueError) as error:
        raise ValueError(f"Cannot read {path}: {error}") from error
    if data.shape[0] < 2 or data.shape[1] < 2:
        raise ValueError(f"{path}: use the complete XVG file with at least two data rows.")
    if not np.isfinite(data).all():
        raise ValueError(f"{path}: NaN or infinite values were found.")
    return data


def make_figure(histograms, profile):
    distance = histograms[:, 0]
    counts = histograms[:, 1:]
    if np.any(np.diff(distance) <= 0):
        raise ValueError("Histogram distances must increase from row to row.")
    if np.any(counts < 0):
        raise ValueError("Histogram counts cannot be negative.")

    us_distance = profile[:, 0]
    reference_index = np.argmax(us_distance)
    us_pmf = profile[:, 1] - profile[reference_index, 1]

    fig, axes = plt.subplots(2, 1, figsize=(10, 12), sharex=True)
    palette = LinearSegmentedColormap.from_list(
        "umbrella_purples", ["#4B0082", "#9370DB", "#DA70D6", "#C71585"]
    )
    colors = palette(np.linspace(0, 1, counts.shape[1]))

    for i in range(counts.shape[1]):
        axes[0].plot(
            distance, counts[:, i], color=colors[i], linewidth=1.3,
            alpha=0.9, label=f"Series {i + 1}",
        )
    axes[0].set_ylabel("Count", fontsize=32)
    axes[0].set_ylim(bottom=0)
    if counts.shape[1] <= 12:
        axes[0].legend(frameon=False, fontsize=10, ncol=2)

    axes[1].plot(
        us_distance, us_pmf, color="black", linestyle="-",
        linewidth=4, label="US",
    )
    axes[1].set_xlabel("Distance (nm)", fontsize=32)
    axes[1].set_ylabel("PMF (kJ/mol)", fontsize=32)
    axes[1].legend(frameon=False, fontsize=20)

    region_label_y = -50
    axes[1].text(
        4.5, region_label_y, "Solution", fontsize=16,
        horizontalalignment="center", verticalalignment="center",
    )
    axes[1].text(
        0.06, region_label_y, "Fibril", fontsize=18,
        transform=axes[1].get_yaxis_transform(),
        horizontalalignment="left", verticalalignment="center",
    )
    axes[1].annotate(
        "", xy=(3.0, 0.0), xytext=(4.5, 0.0),
        arrowprops={"arrowstyle": "->", "linewidth": 3, "color": "black"},
    )

    minimum_index = int(np.argmin(us_pmf))
    minimum_x = float(us_distance[minimum_index])
    minimum_pmf = float(us_pmf[minimum_index])
    x_left, x_right = axes[1].get_xlim()
    label_on_right = minimum_x <= 0.5 * (x_left + x_right)
    axes[1].scatter(
        [minimum_x], [minimum_pmf], color="black", s=40, zorder=5,
    )
    axes[1].annotate(
        f"Min = {minimum_pmf:.2f} kJ/mol",
        xy=(minimum_x, minimum_pmf),
        xytext=(20 if label_on_right else -20, 20),
        textcoords="offset points", fontsize=16,
        ha="left" if label_on_right else "right", va="bottom",
        arrowprops={"arrowstyle": "->", "linewidth": 1.2, "color": "black"},
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85},
    )

    for ax, panel in zip(axes, ("A", "B")):
        ax.text(
            0.0, 1.02, panel, transform=ax.transAxes,
            fontsize=32, fontweight="bold", ha="left", va="bottom",
        )
        ax.tick_params(axis="both", labelsize=32)
        ax.ticklabel_format(axis="both", style="plain", useOffset=False)
        ax.grid(False, which="both", axis="both")

    fig.tight_layout()
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "directory", nargs="?", type=Path,
        help="Folder containing both XVG files; defaults to us/ or the current folder.",
    )
    args = parser.parse_args()
    directory = args.directory
    if directory is None:
        directory = next(
            (d for d in (Path("us"), Path("."))
             if (d / HISTOGRAM_FILE).is_file() and (d / PROFILE_FILE).is_file()),
            None,
        )
        if directory is None:
            parser.error(
                "Both histo23.xvg and profile23.xvg must be in us/ or the current "
                "folder. Alternatively, supply their folder as an argument."
            )
    try:
        histograms = read_xvg(directory / HISTOGRAM_FILE)
        profile = read_xvg(directory / PROFILE_FILE, usecols=(0, 1))
        fig = make_figure(histograms, profile)
    except ValueError as error:
        parser.error(str(error))

    for extension in ("pdf", "png"):
        output = directory / f"{OUTPUT_BASENAME}.{extension}"
        fig.savefig(output, dpi=300, bbox_inches="tight")
        print(f"Saved: {output.resolve()}")
    plt.close(fig)


if __name__ == "__main__":
    main()

