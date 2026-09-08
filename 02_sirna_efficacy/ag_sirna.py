import os, sys, gc, time, json, requests
import numpy as np, pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client

D = os.path.dirname(os.path.abspath(__file__))
LENGTHS = sorted([dna_client.SEQUENCE_LENGTH_16KB, dna_client.SEQUENCE_LENGTH_100KB,
                  dna_client.SEQUENCE_LENGTH_500KB, dna_client.SEQUENCE_LENGTH_1MB])
client = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])

# gene -> (mRNA accession, chrom, gstart, gend)
GENES = {
    "P2RX3": ("NM_002559", "chr11", 57338346, 57374566),
    "UBE2B": ("NM_003337", "chr5", 134371174, 134394447),
    "UBE2S": ("NM_014501", "chr19", 55385311, 55408164),
    "UBE2N": ("NM_003348", "chr12", 93405673, 93442263),
    "UBE2K": ("NM_005339", "chr4",  39697879, 39782792),
}
ONLY = sys.argv[1] if len(sys.argv) > 1 else None

comp = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(comp)[::-1]

def load_fasta(path):
    d, cur = {}, None
    for ln in open(path):
        ln = ln.strip()
        if ln.startswith(">"): cur = ln[1:]; d[cur] = []
        elif cur: d[cur].append(ln)
    return {k: "".join(v).upper() for k, v in d.items()}

mrna = load_fasta(D + "/mRNA_HUVK.fas")
sirna = load_fasta(D + "/siRNA_HUVK.fas")
eff = pd.read_csv(D + "/eff.csv")

def pick_len(span):
    for L in LENGTHS:
        if L >= span + 4000: return L
    return LENGTHS[-1]

def ensembl_seq(chrom, start, end):
    c = chrom.replace("chr", "")
    r = requests.get(f"https://rest.ensembl.org/sequence/region/human/{c}:{start}..{end}",
                     headers={"Content-Type": "text/plain"}, timeout=90)
    r.raise_for_status()
    return r.text.strip().upper()

def reduce_mean(track_data):
    v = np.asarray(track_data.values, dtype=np.float32)
    return v.mean(axis=1)   # 1D over positions

rows = []
for gene, (acc, chrom, gs, ge) in GENES.items():
    if ONLY and gene != ONLY: continue
    t0 = time.time()
    span = ge - gs
    L = pick_len(span)
    mid = (gs + ge) // 2
    istart = mid - L // 2
    iend = istart + L
    gseq = ensembl_seq(chrom, istart + 1, iend)   # L bases, index o -> AG position istart+o
    if len(gseq) != L:
        print(f"[{gene}] seq len {len(gseq)} != {L}, trimming"); gseq = gseq[:L].ljust(L, "N")
    m = mrna.get(acc, "").replace("U", "T")
    ids = eff.loc[eff.mRNA == acc, "siRNA"].tolist()

    # map each siRNA -> mRNA target pos -> genomic offset in gseq
    recs = []
    for sid in ids:
        s = sirna.get(sid)
        if not s: continue
        d = s.replace("U", "T")
        # find sense (mRNA-matching) form
        sense = None
        if d in m: sense = d
        elif rc(d) in m: sense = rc(d)
        else: continue
        # locate in genome (either strand)
        o = gseq.find(sense)
        if o < 0:
            o = gseq.find(rc(sense))
        if o < 0:
            continue   # spans an exon junction -> not contiguous in genome
        e = float(eff.loc[(eff.mRNA == acc) & (eff.siRNA == sid), "efficacy"].iloc[0])
        recs.append((sid, o, len(sense), e))

    n_tot = len(ids); n_map = len(recs)
    if n_map < 10:
        print(f"[{gene}] only {n_map}/{n_tot} mapped, skipping"); continue

    out = client.predict_interval(
        interval=genome.Interval(chromosome=chrom, start=istart, end=iend),
        organism=dna_client.Organism.HOMO_SAPIENS,
        requested_outputs=[dna_client.OutputType.RNA_SEQ,
                           dna_client.OutputType.DNASE,
                           dna_client.OutputType.SPLICE_SITE_USAGE],
        ontology_terms=None)
    rna = reduce_mean(out.rna_seq)
    dnase = reduce_mean(out.dnase)
    try: splice = reduce_mean(out.splice_site_usage)
    except Exception: splice = np.zeros(L, np.float32)
    del out; gc.collect()

    for sid, o, ln, e in recs:
        w = slice(max(0, o), min(L, o + ln))
        rows.append({"gene": gene, "siRNA": sid, "efficacy": e,
                     "rna_seq": float(rna[w].mean()),
                     "dnase": float(dnase[w].mean()),
                     "splice_usage": float(splice[w].mean())})
    print(f"[{gene}] mapped {n_map}/{n_tot}  L={L}  {time.time()-t0:.1f}s")
    del rna, dnase, splice; gc.collect()

df = pd.DataFrame(rows)
df.to_csv(D + "/ag_sirna_features.csv", index=False)
print("\nTOTAL siRNAs with features:", len(df))
if len(df):
    from scipy.stats import spearmanr
    print("\n== pooled Spearman (feature vs efficacy) ==")
    for f in ["rna_seq", "dnase", "splice_usage"]:
        rho, p = spearmanr(df[f], df.efficacy)
        print(f"  {f:14} rho={rho:+.3f}  p={p:.2e}")
    print("\n== per-gene Spearman (rna_seq vs efficacy) ==")
    for g, sub in df.groupby("gene"):
        if len(sub) > 8:
            rho, p = spearmanr(sub.rna_seq, sub.efficacy)
            print(f"  {g:8} n={len(sub):3} rho={rho:+.3f} p={p:.2e}")
