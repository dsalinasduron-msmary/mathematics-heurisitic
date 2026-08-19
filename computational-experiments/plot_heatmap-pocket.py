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

# -- pocket boundaries and labels --
pocket_of   = [int(re.search(r"pocket(\d+)", c).group(1)) for c in combos]
boundaries  = [j for j in range(1, len(pocket_of)) if pocket_of[j] != pocket_of[j-1]]
pocket_ids  = sorted(set(pocket_of))

# midpoint row index for each pocket (for y-axis label placement)
pocket_starts = [0] + boundaries
pocket_ends   = boundaries + [len(combos)]
pocket_mids   = {pid: (pocket_starts[i] + pocket_ends[i] - 1) / 2
                 for i, pid in enumerate(pocket_ids)}

# -- plot --
fig_h = max(8, len(combos) * 0.07)
fig, ax = plt.subplots(figsize=(7, fig_h))

vmin, vmax = -13, 5
im = ax.imshow(matrix, aspect="auto", interpolation="none",
               cmap="RdYlGn", vmin=vmin, vmax=vmax)

# solid white lines between pockets
for b in boundaries:
    ax.axhline(b - 0.5, color="white", linewidth=1.5)

ax.set_xticks(range(n_ligs))
ax.set_xticklabels([f"lig{i}" for i in range(n_ligs)], fontsize=10)
ax.set_xlabel("Ligand stage", fontsize=12)
ax.set_ylabel("Pocket", fontsize=12)
ax.set_title("Docking affinity as imatinib is assembled\n(spheres grouped by pocket)", fontsize=13)

# label each pocket at its midpoint row
ax.set_yticks([pocket_mids[pid] for pid in pocket_ids])
ax.set_yticklabels([f"pocket{pid}" for pid in pocket_ids], fontsize=10)

cbar = fig.colorbar(im, ax=ax, shrink=0.4, pad=0.02)
cbar.set_label("Affinity (kcal/mol)", fontsize=10)
cbar.ax.text(0.5, 1.03, f"(capped at {vmax})", transform=cbar.ax.transAxes,
             ha="center", fontsize=7, color="gray")

plt.tight_layout()
out = "/shared/heatmap-pocket.png"
plt.savefig(out, dpi=180)
print(f"Saved to {out}  ({len(combos)} spheres across {len(pocket_ids)} pockets)")
