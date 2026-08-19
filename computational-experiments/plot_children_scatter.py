#!/usr/bin/env python3
"""Plot scaffold vs child scatter plots for affinity, INTER, and INTRA.

Usage: plot_children_scatter.py <children.tsv> [output_prefix]

Reads the 10-column TSV produced by compare_children.py and saves three
scatter plots (one per metric) as PNG files.
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    sys.exit("Usage: plot_children_scatter.py <children.tsv> [output_prefix]")

tsv_path = sys.argv[1]
prefix   = sys.argv[2] if len(sys.argv) >= 3 else tsv_path.replace(".tsv", "")

# -- load TSV ------------------------------------------------------------------
cols = {"affinity": (3, 7), "inter": (4, 8), "intra": (5, 9)}
data = {k: ([], []) for k in cols}

with open(tsv_path) as f:
    for line in f:
        parts = line.strip().split("\t")
        for metric, (si, ci) in cols.items():
            data[metric][0].append(float(parts[si]))
            data[metric][1].append(float(parts[ci]))

# -- plot ----------------------------------------------------------------------
labels = {
    "affinity": ("Affinity of scaffold (kcal/mol)", "Affinity of child (kcal/mol)", "Affinity"),
    "inter":    ("INTER of scaffold (kcal/mol)",    "INTER of child (kcal/mol)",    "INTER"),
    "intra":    ("INTRA of scaffold (kcal/mol)",    "INTRA of child (kcal/mol)",    "INTRA"),
}

stage = tsv_path  # used in title

for metric, (xlabel, ylabel, title_metric) in labels.items():
    x, y = np.array(data[metric][0]), np.array(data[metric][1])
    frac_better = np.mean(y < x)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x, y, alpha=0.4, s=18, color="steelblue", linewidths=0)

    lim = [min(x.min(), y.min()) - 0.5, max(x.max(), y.max()) + 0.5]
    ax.plot(lim, lim, color="black", linewidth=1, linestyle="--", label="y = x")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")

    ax.set_xlabel(xlabel, fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_title(f"{title_metric}: scaffold vs child\n({tsv_path})", fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, linestyle="--", alpha=0.4)

    ax.text(0.05, 0.95, f"Child better: {frac_better:.1%}",
            transform=ax.transAxes, fontsize=10, va="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray"))

    out = f"{prefix}-{metric}.png"
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"Saved {out}  (n={len(x)}, child better: {frac_better:.1%})")
