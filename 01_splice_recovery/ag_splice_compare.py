import os, json, time, requests
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from alphagenome.data import genome
from alphagenome.models import dna_client

OUTDIR = os.path.dirname(os.path.abspath(__file__))
LENGTHS = sorted([dna_client.SEQUENCE_LENGTH_16KB, dna_client.SEQUENCE_LENGTH_100KB,
                  dna_client.SEQUENCE_LENGTH_500KB, dna_client.SEQUENCE_LENGTH_1MB])
GENES = ["SOD1", "SMN2", "TTR"]
TOL = 3          # bp tolerance matching predicted peak to annotated boundary
THRESH = 0.5     # predicted probability call threshold

client = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])

def ensembl_canonical(sym):
    r = requests.get(
        f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{sym}?expand=1",
        headers={"Content-Type": "application/json"}, timeout=30)
    r.raise_for_status()
    d = r.json()
    ts = d["Transcript"]
    can = [t for t in ts if t.get("is_canonical")] or ts
    t = can[0]
    exons = sorted(t["Exon"], key=lambda e: e["start"])
    return {"chrom": "chr" + str(d["seq_region_name"]), "start": d["start"],
            "end": d["end"], "strand": d["strand"], "tx": t["id"], "exons": exons}

def pick_len(span):
    for L in LENGTHS:
        if L >= span + 4000:
            return L
    return LENGTHS[-1]

def junctions(g):
    # returns donor and acceptor genomic positions on the gene's coding strand
    ex = g["exons"]
    starts = [e["start"] for e in ex]  # 1-based
    ends = [e["end"] for e in ex]
    if g["strand"] == 1:
        donors = ends[:-1]            # 5' splice site = exon end (+ strand)
        acceptors = starts[1:]        # 3' splice site = exon start (+ strand)
    else:
        donors = starts[1:]           # on - strand, donor is at exon start (lower coord) of downstream exon
        acceptors = ends[:-1]
    return donors, acceptors

summary = []
for sym in GENES:
    g = ensembl_canonical(sym)
    span = g["end"] - g["start"]
    L = pick_len(span)
    mid = (g["start"] + g["end"]) // 2
    istart = mid - L // 2
    iv = genome.Interval(chromosome=g["chrom"], start=istart, end=istart + L)
    out = client.predict_interval(
        interval=iv, organism=dna_client.Organism.HOMO_SAPIENS,
        requested_outputs=[dna_client.OutputType.SPLICE_SITES], ontology_terms=None)
    ss = out.splice_sites
    vals = np.asarray(ss.values)          # (L, 4)
    md = ss.metadata.reset_index(drop=True)
    # find track columns
    def col(name, strand):
        for i, row in md.iterrows():
            if row["name"] == name and row["strand"] == strand:
                return i
        return None
    plus = g["strand"] == 1
    d_strand = "+" if plus else "-"
    di = col("donor", d_strand); ai = col("acceptor", d_strand)
    donor_tr = vals[:, di]; acc_tr = vals[:, ai]

    donors, acceptors = junctions(g)
    def recall(sites, track):
        got = 0; preds = []
        for pos in sites:
            idx = pos - istart - 1     # 1-based -> 0-based index
            lo, hi = max(0, idx - TOL), min(len(track), idx + TOL + 1)
            peak = float(track[lo:hi].max()) if hi > lo else 0.0
            preds.append(peak)
            if peak >= THRESH:
                got += 1
        return got, len(sites), preds
    dg, dn, dpk = recall(donors, donor_tr)
    ag, an, apk = recall(acceptors, acc_tr)
    # background: mean predicted prob away from any annotated site
    all_sites = set()
    for pos in list(donors) + list(acceptors):
        for k in range(-TOL, TOL + 1):
            all_sites.add(pos - istart - 1 + k)
    mask = np.ones(len(donor_tr), bool)
    for k in all_sites:
        if 0 <= k < len(mask):
            mask[k] = False
    bg = float(np.concatenate([donor_tr[mask], acc_tr[mask]]).mean())

    summary.append({"gene": sym, "tx": g["tx"], "chrom": g["chrom"],
                    "strand": "+" if plus else "-", "window_bp": L, "n_exons": len(g["exons"]),
                    "donor_recall": f"{dg}/{dn}", "acceptor_recall": f"{ag}/{an}",
                    "mean_prob_true_donors": round(np.mean(dpk), 3),
                    "mean_prob_true_acceptors": round(np.mean(apk), 3),
                    "background_mean_prob": round(bg, 5)})

    # figure: zoom to gene span
    xs = np.arange(istart, istart + L)
    gsel = (xs >= g["start"] - 300) & (xs <= g["end"] + 300)
    fig, ax = plt.subplots(figsize=(11, 3.2))
    ax.plot(xs[gsel], donor_tr[gsel], lw=0.8, label="AlphaGenome donor P")
    ax.plot(xs[gsel], acc_tr[gsel], lw=0.8, label="AlphaGenome acceptor P")
    for pos in donors:
        ax.axvline(pos, color="tab:blue", ls="--", lw=0.6, alpha=0.6)
    for pos in acceptors:
        ax.axvline(pos, color="tab:orange", ls="--", lw=0.6, alpha=0.6)
    ax.set_title(f"{sym} ({g['tx']}, {g['chrom']} {'+' if plus else '-'}) "
                 f"canonical GENCODE junctions (dashed) vs AlphaGenome predicted splice sites")
    ax.set_xlabel("genomic position (GRCh38)"); ax.set_ylabel("predicted probability")
    ax.set_ylim(-0.02, 1.05); ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, f"splice_{sym}.png"), dpi=140)
    plt.close(fig)
    print("done", sym, summary[-1])
    time.sleep(0.5)

with open(os.path.join(OUTDIR, "splice_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("\nSUMMARY")
for s in summary:
    print(s)
