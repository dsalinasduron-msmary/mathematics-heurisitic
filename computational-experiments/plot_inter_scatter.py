import matplotlib.pyplot as plt
import numpy as np

x, y = [], []
with open("/shared/playout_stage5-inter.tsv") as f:
    for line in f:
        parts = line.strip().split("\t")
        x.append(float(parts[3]))
        y.append(float(parts[5]))

x, y = np.array(x), np.array(y)

fig, ax = plt.subplots(figsize=(6, 6))

ax.scatter(x, y, alpha=0.4, s=18, color="steelblue", linewidths=0)

# diagonal y = x reference line
lim = [min(x.min(), y.min()) - 1, max(x.max(), y.max()) + 1]
ax.plot(lim, lim, color="black", linewidth=1, linestyle="--", label="y = x")
ax.set_xlim(lim)
ax.set_ylim(lim)

ax.set_xlabel("INTER of scaffold (kcal/mol)", fontsize=13)
ax.set_ylabel("INTER of child (kcal/mol)", fontsize=13)
ax.set_title("Change in INTER force upon fragment addition\n(playout stage 5)", fontsize=13)
ax.legend(fontsize=10)

# annotate fraction where child improves on scaffold
frac_better = np.mean(y < x)
ax.text(0.05, 0.95, f"Child better: {frac_better:.1%}",
        transform=ax.transAxes, fontsize=10, va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray"))

ax.grid(True, linestyle="--", alpha=0.4)
ax.set_aspect("equal")

plt.tight_layout()
out = "/shared/inter_scatter_stage5.png"
plt.savefig(out, dpi=150)
print(f"Saved to {out}  (n={len(x)}, child better: {frac_better:.1%})")
