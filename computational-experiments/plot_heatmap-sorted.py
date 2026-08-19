import os
import re
import numpy as np
import matplotlib.pyplot as plt

# -- collect data --
combos = sorted(
    [e.name.replace("lig0_", "") for e in os.scandir("/shared")
     if re.match(r"lig0_pocket\d+_sphere\d+$", e.name)],
    key=lambda s: (int(re.search(r"pocket(\d+)", s).group(1)),
                   int(re.search(r"sphere(\d+)", s).group(1)))
)

n_ligs  = 10
matrix  = np.zeros((len(combos), n_ligs))

for row, combo in enumerate(combos):
    for lig in range(n_ligs):
        path = f"/shared/lig{lig}_{combo}/1iep_ligand_vina_out.pdbqt"
        with open(path) as f:
            for line in f:
                m = re.match(r"REMARK VINA RESULT:\s+([-\d.]+)", line)
                if m:
                    matrix[row, lig] = float(m.group(1))
                    break

# -- sort rows by lig9 affinity (most negative first) --
order         = np.argsort(matrix[:, -1])   # ascending = most negative first
matrix_sorted = matrix[order]

canonical     = "pocket0_sphere24"
canonical_row = list(order).index(combos.index(canonical))

# -- plot --
fig_h = max(8, len(combos) * 0.07)
fig, ax = plt.subplots(figsize=(7, fig_h))

vmin, vmax = -13, 5
im = ax.imshow(matrix_sorted, aspect="auto", interpolation="none",
               cmap="RdYlGn", vmin=vmin, vmax=vmax)

ax.set_xticks(range(n_ligs))
ax.set_xticklabels([f"lig{i}" for i in range(n_ligs)], fontsize=10)
ax.set_xlabel("Ligand stage", fontsize=12)
ax.set_ylabel("Spheres (sorted by lig9 affinity)", fontsize=12)
ax.set_title("Docking affinity as imatinib is assembled\n(rows sorted by lig9 affinity, best at top)", fontsize=13)
ax.set_yticks([])

# mark the canonical imatinib sphere
ax.annotate("← pocket0_sphere24\n   (canonical)",
            xy=(n_ligs - 0.5, canonical_row),
            xytext=(n_ligs + 0.15, canonical_row),
            xycoords=("data", "data"),
            fontsize=7, va="center", color="black",
            arrowprops=dict(arrowstyle="-", color="black", lw=1.0))

cbar = fig.colorbar(im, ax=ax, shrink=0.4, pad=0.02)
cbar.set_label("Affinity (kcal/mol)", fontsize=10)
cbar.ax.text(0.5, 1.03, f"(capped at {vmax})", transform=cbar.ax.transAxes,
             ha="center", fontsize=7, color="gray")

plt.tight_layout()
out = "/shared/heatmap-sorted.png"
plt.savefig(out, dpi=180)
print(f"Saved to {out}  ({len(combos)} spheres)")
