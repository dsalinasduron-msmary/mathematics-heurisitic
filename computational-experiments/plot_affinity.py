import os
import re
import matplotlib.pyplot as plt

pocket = "pocket0_sphere24"
ligands = range(10)

lig_nums = []
affinities = []

for i in ligands:
    folder = f"/shared/lig{i}_{pocket}"
    pdbqt = os.path.join(folder, "1iep_ligand_vina_out.pdbqt")
    if not os.path.exists(pdbqt):
        print(f"Missing: {pdbqt}")
        continue
    with open(pdbqt) as f:
        for line in f:
            m = re.match(r"REMARK VINA RESULT:\s+([-\d.]+)", line)
            if m:
                affinity = float(m.group(1))
                lig_nums.append(i)
                affinities.append(affinity)
                break

print("Ligand | Affinity (kcal/mol)")
for n, a in zip(lig_nums, affinities):
    print(f"  lig{n}  |  {a:.3f}")

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(lig_nums, affinities, marker="o", linewidth=2, markersize=8, color="steelblue")
ax.set_xlabel("Ligand number (assembly stage)", fontsize=13)
ax.set_ylabel("Affinity (kcal/mol)", fontsize=13)
ax.set_title("Imatinib assembly: docking affinity vs. build stage\n(pocket0_sphere24)", fontsize=14)
ax.set_xticks(lig_nums)
ax.grid(True, linestyle="--", alpha=0.5)
ax.invert_yaxis()  # more negative = better, so flip so "better" is up

plt.tight_layout()
out = "/shared/affinity_pocket0_sphere24.png"
plt.savefig(out, dpi=150)
print(f"\nPlot saved to {out}")
