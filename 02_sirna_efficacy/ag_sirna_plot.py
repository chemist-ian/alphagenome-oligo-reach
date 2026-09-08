import os, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

D = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(D + "/ag_sirna_features.csv")
feats = ["rna_seq", "dnase", "splice_usage"]
labels = {"rna_seq": "predicted RNA-seq coverage",
          "dnase": "predicted DNase (chromatin) — NEG CONTROL",
          "splice_usage": "predicted splice-site usage"}

summary = {"n": int(len(df)), "genes": {}, "pooled": {}}
for f in feats:
    rho, p = spearmanr(df[f], df.efficacy)
    summary["pooled"][f] = {"rho": round(float(rho), 3), "p": float(p)}
for g, sub in df.groupby("gene"):
    d = {"n": int(len(sub))}
    for f in feats:
        if sub[f].nunique() > 1:
            rho, p = spearmanr(sub[f], sub.efficacy)
            d[f] = {"rho": round(float(rho), 3), "p": float(p)}
    summary["genes"][g] = d
json.dump(summary, open(D + "/ag_sirna_summary.json", "w"), indent=2)

# Figure 1: pooled scatter for each feature
fig, axes = plt.subplots(1, 3, figsize=(13, 4))
for ax, f in zip(axes, feats):
    ax.scatter(df[f], df.efficacy, s=10, alpha=0.5)
    rho, p = spearmanr(df[f], df.efficacy)
    ax.set_title(f"{labels[f]}\nSpearman rho={rho:+.2f}, p={p:.1e}", fontsize=9)
    ax.set_xlabel(f); ax.set_ylabel("measured siRNA efficacy")
fig.suptitle(f"AlphaGenome predicted signal vs measured siRNA knockdown  (n={len(df)} siRNAs, 5 human genes)",
             fontsize=11)
fig.tight_layout()
fig.savefig(D + "/ag_sirna_pooled.png", dpi=140)
plt.close(fig)

# Figure 2: per-gene Spearman rho bars (rna_seq vs dnase control)
genes = sorted(summary["genes"])
x = np.arange(len(genes)); w = 0.35
rna = [summary["genes"][g].get("rna_seq", {}).get("rho", 0) for g in genes]
dna = [summary["genes"][g].get("dnase", {}).get("rho", 0) for g in genes]
fig, ax = plt.subplots(figsize=(9, 4))
ax.bar(x - w/2, rna, w, label="RNA-seq coverage")
ax.bar(x + w/2, dna, w, label="DNase (neg control)")
ax.axhline(0, color="k", lw=0.8)
ax.set_xticks(x); ax.set_xticklabels(genes)
ax.set_ylabel("Spearman rho vs efficacy"); ax.set_ylim(-0.6, 0.6)
ax.set_title("Per-gene correlation of AlphaGenome signal with siRNA efficacy")
ax.legend()
fig.tight_layout()
fig.savefig(D + "/ag_sirna_pergene.png", dpi=140)
plt.close(fig)

print(json.dumps(summary, indent=2))
