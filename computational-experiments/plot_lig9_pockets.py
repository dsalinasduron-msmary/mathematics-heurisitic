import os
import re
import matplotlib.pyplot as plt
import numpy as np

data = {}  # pocket_label -> list of affinities

for entry in os.scandir("/shared"):
    m = re.match(r"lig9_(pocket\d+)_sphere\d+$", entry.name)
    if not m:
        continue
    pocket = m.group(1)
    pdbqt = os.path.join(entry.path, "1iep_ligand_vina_out.pdbqt")
    if not os.path.exists(pdbqt):
        continue
    with open(pdbqt) as f:
        for line in f:
            hit = re.match(r"REMARK VINA RESULT:\s+([-\d.]+)", line)
            if hit:
                data.setdefault(pocket, []).append(float(hit.group(1)))
                break

pockets = sorted(data.keys(), key=lambda p: int(p.replace("pocket", "")))
values  = [data[p] for p in pockets]

print("Pocket | n spheres | min affinity | median affinity")
for p, v in zip(pockets, values):
    print(f"  {p}  |  {len(v):>5}  |  {min(v):>8.3f}    |  {np.median(v):>8.3f}")

fig, ax = plt.subplots(figsize=(10, 5))

parts = ax.violinplot(values, positions=range(len(pockets)),
                      showmedians=True, showextrema=True)

parts["cmedians"].set_color("black")
parts["cmedians"].set_linewidth(2)
for body in parts["bodies"]:
    body.set_alpha(0.6)
    body.set_facecolor("steelblue")

ax.set_xticks(range(len(pockets)))
ax.set_xticklabels(pockets, fontsize=12)
ax.set_xlabel("Pocket", fontsize=13)
ax.set_ylabel("Affinity (kcal/mol)", fontsize=13)
ax.set_title("Lig9 (imatinib) docking affinity distribution per pocket\n(one data point per sphere)", fontsize=14)
ax.invert_yaxis()
ax.grid(axis="y", linestyle="--", alpha=0.5)

# annotate n per pocket
for i, v in enumerate(values):
    ax.text(i, ax.get_ylim()[0] - 0.3, f"n={len(v)}", ha="center", fontsize=9, color="gray")

plt.tight_layout()
out = "/shared/lig9_affinity_by_pocket.png"
plt.savefig(out, dpi=150)
print(f"\nPlot saved to {out}")
