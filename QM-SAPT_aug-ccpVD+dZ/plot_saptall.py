#!/usr/bin/env python3
"""Author: Najla Hosseini
Email: najla.hosseini@gmail.com| najla.hosseini@ki.se
Affiliation: Karolinska Institutet,
SE-171 77 Stockholm, Sweden

Description:Read the two SAPT CSV files and plot them."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator







ION_FILE = Path("sapt_summary.csv")
PEPTIDE_FILE = Path("PEP_SAPT0.csv")
OUTPUT_DIR = Path("sapt_figure")

ION_FRAMES = [35,36,37,38,39]
PEPTIDE_IDS = ("NNQQ1", "NNQQ2")
MATCH_ION_Y_LIMITS = True
DPI = 600

COMPONENTS = {
    "Electrostatics": "Electrostatics (kJ/mol)",
    "Exchange": "Exchange (kJ/mol)",
    "Induction": "Induction (kJ/mol)",
    "Dispersion": "Dispersion (kJ/mol)",
}

ENERGY_COLS = list(COMPONENTS.values()) + [
    "Total SAPT (kJ/mol)"
]


ALPHA = 1
COLORS = {
    "Electrostatics": "#A580C1",  # Pastel purple
    "Exchange":      "#C9B8ED",  # Lavender
    "Induction":     "#EDB8EB",  # Pastel orchid
    "Dispersion":    "#E38AC2",  # Soft pink
}

STYLE = {
    "font.size": 24,
    "text.usetex": False,
    "text.color": "black",
    "axes.labelcolor": "black",
    "axes.titlecolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "axes.titlesize": 28,
    "axes.labelsize": 24,
    "axes.linewidth": 0.9,
    "axes.grid": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}



def read_table(path):
    if not path.is_file():
        raise FileNotFoundError(
            f"Cannot find {path.name} in the current directory:\n"
            f"{Path.cwd()}"
        )

    df = pd.read_csv(
        path,
        encoding="utf-8-sig",
        dtype=str,
        keep_default_na=False,
    )

    df.columns = df.columns.str.strip()

    required = ["Ion", "Frame", "Directory"] + ENERGY_COLS
    missing = [
        column for column in required
        if column not in df.columns
    ]

    if missing:
        raise ValueError(
            f"{path.name}: missing columns: {missing}"
        )

    if df.empty:
        raise ValueError(f"{path.name}: no data rows.")

    for column in ("Ion", "Frame"):
        df[column] = df[column].str.strip()

        if df[column].eq("").any():
            raise ValueError(
                f"{path.name}: empty {column} labels."
            )

    for column in ENERGY_COLS:
        df[column] = pd.to_numeric(
            df[column], errors="coerce"
        )

    invalid = ~np.isfinite(
        df[ENERGY_COLS].to_numpy(dtype=float)
    )

    if invalid.any():
        details = [
            f"CSV row {r + 2}: {ENERGY_COLS[c]}"
            for r, c in np.argwhere(invalid)
        ]

        raise ValueError(
            f"{path.name}: invalid energies:\n"
            + "\n".join(details)
        )

    return df


def load_inputs():
    ions = read_table(ION_FILE)

    aliases = {
        "na": "Na",
        "na+": "Na",
        "k": "K",
        "k+": "K",
    }

    labels = (
        ions["Ion"]
        .str.replace("⁺", "+", regex=False)
        .str.lower()
    )

    ions["Ion"] = labels.map(aliases)

    if ions["Ion"].isna().any():
        raise ValueError(
            f"{ION_FILE.name}: expected Na or K "
            "in the Ion column."
        )

    frames = pd.to_numeric(
        ions["Frame"], errors="coerce"
    )

    if (
        not np.isfinite(frames.to_numpy(dtype=float)).all()
        or (frames % 1 != 0).any()
    ):
        raise ValueError(
            f"{ION_FILE.name}: ion Frame values "
            "must be integers."
        )

    ions["Frame"] = frames.astype(int)

    selected = (
        sorted(ions["Frame"].unique())
        if ION_FRAMES is None
        else list(ION_FRAMES)
    )

    if not selected or len(selected) != len(set(selected)):
        raise ValueError(
            "ION_FRAMES must be nonempty "
            "and contain no duplicates."
        )

    ions = ions[
        ions["Frame"].isin(selected)
    ].copy()

    if ions.duplicated(["Ion", "Frame"]).any():
        raise ValueError(
            "More than one ion calculation per Ion/Frame; "
            "select the intended rows first."
        )

    ion_groups = {}

    for ion in ("Na", "K"):
        group = (
            ions[ions["Ion"] == ion]
            .set_index("Frame")
        )

        absent = [
            frame for frame in selected
            if frame not in group.index
        ]

        if absent:
            raise ValueError(
                f"{ION_FILE.name}: missing "
                f"{ion} frames {absent}."
            )

        group = group.loc[selected].reset_index()

        group["Conformation"] = [
            f"C{i + 1}"
            for i in range(len(selected))
        ]

        ion_groups[ion] = group

    peptides = read_table(PEPTIDE_FILE)

    if (
        len(peptides) != 2
        or set(peptides["Ion"]) != set(PEPTIDE_IDS)
    ):
        raise ValueError(
            f"{PEPTIDE_FILE.name}: expected exactly "
            "one NNQQ1 row and one NNQQ2 row."
        )

    peptide_rows = {
        name: peptides.loc[
            peptides["Ion"] == name
        ].iloc[0]
        for name in PEPTIDE_IDS
    }

    return ion_groups, peptide_rows



def energy_limits(values, padding=0.13):
    values = np.asarray(values, dtype=float)

    low = min(0.0, values.min())
    high = max(0.0, values.max())
    margin = padding * max(high - low, 1.0)

    return low - margin, high + margin


def format_axis(ax, letter, title, ylim):
    ax.set_title(title, pad=17)
    ax.set_ylim(*ylim)

    ax.set_ylabel(
        r"SAPT energy (kJ mol$^{-1}$)",
        labelpad=9,
    )

    ax.yaxis.set_major_locator(
        MaxNLocator(nbins=6)
    )

    ax.ticklabel_format(
        axis="y",
        style="plain",
        useOffset=False,
    )

    ax.tick_params(
        axis="both",
        pad=6,
        length=4,
        width=0.9,
    )

    ax.axhline(
        0,
        color="black",
        linewidth=0.8,
        alpha=0.75,
        zorder=2,
    )

    ax.minorticks_off()

    ax.text(
        -0.13,
        1.045,
        letter,
        transform=ax.transAxes,
        fontsize=21,
        fontweight="bold",
        color="black",
        ha="left",
        va="bottom",
        clip_on=False,
    )


def draw_ion(ax, df, ion, letter, ylim):
    x = np.arange(len(df))
    width = 0.185

    for i, (name, column) in enumerate(COMPONENTS.items()):
        ax.bar(
            x + (i - 1.5) * width,
            df[column].to_numpy(),
            width=width,
            color=COLORS[name],
            alpha=ALPHA,
            edgecolor="none",
            zorder=3,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(df["Conformation"])
    ax.set_xlim(-0.6, len(df) - 0.4)
    ax.set_xlabel("Conformation", labelpad=9)

    format_axis(
        ax,
        letter,
        rf"NNQQ–{ion}$^+$",
        ylim,
    )


def draw_peptide(ax, row, letter, ylim):
    values = row[
        list(COMPONENTS.values())
    ].to_numpy(dtype=float)

    bars = ax.bar(
        np.arange(4),
        values,
        width=0.62,
        color=[
            COLORS[name] for name in COMPONENTS
        ],
        alpha=ALPHA,
        edgecolor="none",
        zorder=3,
    )

    ax.set_xticks(np.arange(4))
    ax.set_xticklabels(list(COMPONENTS))
    ax.set_xlim(-0.6, 3.6)
    ax.set_xlabel("SAPT component", labelpad=9)

    format_axis(
        ax,
        letter,
        f"{row['Ion']}-{row['Ion']}",
        #f"{row['Ion']} (frame {row['Frame']})",
        ylim,
    )

    for bar, value in zip(bars, values):
        ax.annotate(
            f"{value:.1f}",
            (
                bar.get_x() + bar.get_width() / 2,
                value,
            ),
            xytext=(0, 6 if value >= 0 else -6),
            textcoords="offset points",
            ha="center",
            va="bottom" if value >= 0 else "top",
            fontsize=11,
        )

    total = row["Total SAPT (kJ/mol)"]

    ax.text(
        0.98,
        0.96,
        rf"Total SAPT = {total:.2f} kJ mol$^{{-1}}$",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=11,
    )



def make_figure(ion_groups, peptide_rows):
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(13.2, 9.5),
    )

    fig.subplots_adjust(
        left=0.09,
        right=0.985,
        bottom=0.10,
        top=0.85,
        hspace=0.55,
        wspace=0.30,
    )

    ion_values = np.concatenate([
        group[
            list(COMPONENTS.values())
        ].to_numpy().ravel()
        for group in ion_groups.values()
    ])

    for column, ion in enumerate(("Na", "K")):
        values = (
            ion_values
            if MATCH_ION_Y_LIMITS
            else ion_groups[ion][
                list(COMPONENTS.values())
            ]
        )

        draw_ion(
            axes[0, column],
            ion_groups[ion],
            ion,
            "AB"[column],
            energy_limits(values),
        )

    peptide_values = np.concatenate([
        row[
            list(COMPONENTS.values())
        ].to_numpy(dtype=float)
        for row in peptide_rows.values()
    ])

    shared_limits = energy_limits(
        peptide_values,
        padding=0.18,
    )

    for column, name in enumerate(PEPTIDE_IDS):
        draw_peptide(
            axes[1, column],
            peptide_rows[name],
            "CD"[column],
            shared_limits,
        )

    handles = [
        Patch(
            facecolor=COLORS[name],
            alpha=ALPHA,
            edgecolor="none",
            label=name,
        )
        for name in COMPONENTS
    ]

    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.54, 0.985),
        ncol=4,
        frameon=False,
        fontsize=11.5,
        columnspacing=2.1,
        handlelength=1.5,
        handletextpad=0.65,
    )

    fig.align_ylabels(axes[:, 0])
    fig.align_ylabels(axes[:, 1])

    return fig



def main():
    ion_groups, peptide_rows = load_inputs()

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    with plt.rc_context():
        plt.rcdefaults()
        plt.rcParams.update(STYLE)

        fig = make_figure(
            ion_groups,
            peptide_rows,
        )

        try:
            for extension in ("pdf", "svg", "png"):
                path = (
                    OUTPUT_DIR
                    / f"SAPT_abcd.{extension}"
                )

                fig.savefig(
                    path,
                    dpi=DPI,
                    bbox_inches="tight",
                    pad_inches=0.15,
                )

                print(f"Saved: {path.resolve()}")

        finally:
            plt.close(fig)

    exports = []

    for panel, ion in (("A", "Na"), ("B", "K")):
        data = ion_groups[ion].copy()
        data.insert(0, "Panel", panel)
        data.insert(1, "Source file", ION_FILE.name)
        exports.append(data)

    for panel, name in zip(("C", "D"), PEPTIDE_IDS):
        data = peptide_rows[name].to_frame().T.copy()
        data.insert(0, "Panel", panel)
        data.insert(1, "Source file", PEPTIDE_FILE.name)
        data["Conformation"] = data["Frame"]
        exports.append(data)

    pd.concat(
        exports,
        ignore_index=True,
    ).to_csv(
        OUTPUT_DIR / "SAPT_abcd_plotted_values.csv",
        index=False,
    )

    print("a: Na; b: K; c: NNQQ1; d: NNQQ2.")

    print(
        "Ion frames: "
        + ", ".join(
            map(str, ion_groups["Na"]["Frame"])
        )
    )

    print(
        "a/b y limits: "
        + (
            "matched"
            if MATCH_ION_Y_LIMITS
            else "independent"
        )
    )

    print("c/d y limits: matched. No grid lines.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as exc:
        raise SystemExit(f"Error: {exc}")
