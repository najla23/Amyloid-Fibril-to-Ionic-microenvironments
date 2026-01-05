#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

############################################
# INPUT FILES
############################################
sapt_csv = "sapt_summary.csv"
pt_dist_csv = "extracted_nnqq_atomID_pt/Pt_COM_distances.csv"
na_dist_csv = "extracted_nnqq_atomID_na/Na_COM_distances.csv"

############################################
# LOAD DATA
############################################
sapt = pd.read_csv(sapt_csv)
pt_dist = pd.read_csv(pt_dist_csv)
na_dist = pd.read_csv(na_dist_csv)

sapt["Frame"] = sapt["Frame"].astype(int)
pt_dist["frame"] = pt_dist["frame"].astype(int)
na_dist["frame"] = na_dist["frame"].astype(int)

############################################
# SPLIT + MERGE
############################################
sapt_k = sapt[sapt["Ion"] == "K"].copy()
sapt_na = sapt[sapt["Ion"] == "Na"].copy()

sapt_k = sapt_k.merge(
    pt_dist[["frame", "distance_COM_Pt_A"]],
    left_on="Frame", right_on="frame", how="left"
)

sapt_na = sapt_na.merge(
    na_dist[["frame", "distance_COM_Na_A"]],
    left_on="Frame", right_on="frame", how="left"
)

sapt_k.sort_values("Frame", inplace=True)
sapt_na.sort_values("Frame", inplace=True)

############################################
# RELABEL FRAMES (35–40 → 1–6)
############################################
for df in (sapt_k, sapt_na):
    df["PlotFrame"] = range(1, len(df) + 1)

############################################
# SAPT COMPONENTS
############################################
components = [
    ("Electrostatics (kJ/mol)", "Electrostatics"),
    ("Exchange (kJ/mol)", "Exchange"),
    ("Induction (kJ/mol)", "Induction"),
    ("Dispersion (kJ/mol)", "Dispersion"),
]

colors = {
    "Electrostatics": "#4B0082",
    "Exchange": "#9370DB",
    "Induction": "#DA70D6",
    "Dispersion": "#C71585",
}

############################################
# GLOBAL Y LIMITS (same scale)
############################################
all_vals = pd.concat([
    sapt_k[[c[0] for c in components]],
    sapt_na[[c[0] for c in components]],
])

ymin = all_vals.min().min() * 1.1
ymax = all_vals.max().max() * 1.1

############################################
# BAR PLOT FUNCTION (NO DISTANCES)
############################################
def plot_components(ax, df, ion_label, panel_label):
    x = np.arange(len(df))
    width = 0.18

    for i, (col, label) in enumerate(components):
        ax.bar(
            x + i * width,
            df[col],
            width=width,
            color=colors[label],
            label=label,
            alpha=0.85
        )

    ax.set_xticks(x + width * 1.5)
    ax.set_xticklabels(df["PlotFrame"], fontsize=14)
    ax.set_ylabel("SAPT Energy (kJ/mol)", fontsize=16)
    ax.set_title(f"{ion_label}–NNQQ SAPT Energy Components", fontsize=18)
    ax.set_ylim(ymin, ymax)
    ax.tick_params(axis="y", labelsize=14)
    ax.legend(frameon=False, fontsize=12, ncol=2, loc="lower left")

    ax.text(
        0.005, 0.93, panel_label,
        transform=ax.transAxes,
        fontsize=24,
        fontweight="bold",
        va="top",
        ha="left"
    )



############################################
# LaTeX TABLE: SAPT COMPONENTS PER CONFORMATION + AVERAGES
############################################

energy_cols = [
    "Electrostatics (kJ/mol)",
    "Exchange (kJ/mol)",
    "Induction (kJ/mol)",
    "Dispersion (kJ/mol)",
    "Total SAPT (kJ/mol)",
]

def build_sapt_table(df, ion_label):
    table = df.copy()

    table["Conformation"] = [f"C{i}" for i in table["PlotFrame"]]

    table = table[["Conformation"] + energy_cols]

    avg_row = table[energy_cols].mean().to_frame().T
    avg_row["Conformation"] = "Average"

    table = pd.concat([table, avg_row], ignore_index=True)

    table.insert(0, "Ion", ion_label)

    return table


table_k = build_sapt_table(sapt_k, "K$^+$")
table_na = build_sapt_table(sapt_na, "Na$^+$")

full_table = pd.concat([table_k, table_na], ignore_index=True)

full_table.rename(columns={
    "Electrostatics (kJ/mol)": r"$E_{\mathrm{elst}}$ (kJ/mol)",
    "Exchange (kJ/mol)": r"$E_{\mathrm{exch}}$ (kJ/mol)",
    "Induction (kJ/mol)": r"$E_{\mathrm{ind}}$ (kJ/mol)",
    "Dispersion (kJ/mol)": r"$E_{\mathrm{disp}}$ (kJ/mol)",
    "Total SAPT (kJ/mol)": r"$E_{\mathrm{tot}}$ (kJ/mol)",
}, inplace=True)

############################################
# EXPORT LaTeX TABLE
############################################
latex_table = full_table.to_latex(
    index=False,
    float_format="%.2f",
    escape=False,
    column_format="llccccc",
    caption=(
        "Symmetry-adapted perturbation theory (SAPT) energy decomposition for "
        "NNQQ--ion interactions. Individual rows correspond to conformations "
        "C1--C6 selected from steered molecular dynamics simulations, followed by "
        "the frame-averaged interaction energies for each ion."
    ),
    label="tab:sapt_components_conformations"
)

with open("SAPT_components_conformations.tex", "w") as f:
    f.write(latex_table)

print("LaTeX table written to SAPT_components_conformations.tex")


############################################
# FIGURE: SAPT BAR PLOTS
############################################
fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

plot_components(axes[0], sapt_k, ion_label="K⁺", panel_label="A")
plot_components(axes[1], sapt_na, ion_label="Na⁺", panel_label="B")

axes[1].set_xlabel("Conformations", fontsize=16)

plt.tight_layout()
plt.savefig("SAPT_components_bars.pdf", dpi=300, bbox_inches="tight")
plt.show()

