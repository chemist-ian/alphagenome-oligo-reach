import os, requests, numpy as np, pandas as pd, time
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from alphagenome.data import genome
from alphagenome.models import dna_client

OUT = os.path.dirname(os.path.abspath(__file__))
client = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])
CHROM="chrX"; ES=31773960; EE=31774192   # DMD exon51, minus strand: donor=ES(5'SS), acceptor=EE(3'SS)
L=dna_client.SEQUENCE_LENGTH_100KB
center=(ES+EE)//2; istart=center-L//2; iend=istart+L
iv=genome.Interval(chromosome=CHROM,start=istart,end=iend)
comp=str.maketrans("ACGT","TGCA"); rc=lambda s:s.translate(comp)[::-1]
TRANS=str.maketrans("ACGT","GTAC")

def ens(a,b):
    r=requests.get(f"https://rest.ensembl.org/sequence/region/human/X:{a}..{b}",
       headers={"Content-Type":"text/plain"},timeout=120); r.raise_for_status(); return r.text.strip().upper()
ref=ens(istart+1,iend); assert len(ref)==L

def cols(md):
    dc=ac=None
    for i,r in md.reset_index(drop=True).iterrows():
        if r["name"]=="donor" and r["strand"]=="-": dc=i
        if r["name"]=="acceptor" and r["strand"]=="-": ac=i
    return dc,ac
def inclusion(seq):
    out=client.predict_sequence(sequence=seq,organism=dna_client.Organism.HOMO_SAPIENS,
        requested_outputs=[dna_client.OutputType.SPLICE_SITE_USAGE],ontology_terms=None,interval=iv)
    ss=out.splice_site_usage; v=np.asarray(ss.values); dc,ac=cols(ss.metadata)
    ae=EE-istart; de=ES-istart
    return (v[max(0,ae-3):ae+4,ac].max()+v[max(0,de-3):de+4,dc].max())/2

base=inclusion(ref)

# locate eteplirsen target in exon (annotation)
etep="CTCCAACATCAAGGAAGATGGCATTTCTAG"
et_span=None
for cand in [etep, rc(etep)]:
    o=ref.find(cand)
    if o>=0: et_span=(istart+o, istart+o+len(cand)); break

# sweep masks across exon +/- 60 nt flanks
WIN=20; STEP=4
lo, hi = ES-60, EE+60
rows=[]
for gpos in range(lo, hi-WIN, STEP):
    o1=gpos-istart; o2=o1+WIN
    mut=ref[:o1]+ref[o1:o2].translate(TRANS)+ref[o2:]
    inc=inclusion(mut)
    ctr=gpos+WIN//2
    in_exon = ES<=ctr<=EE
    rows.append({"center":ctr, "rel_exon5p": EE-ctr, "in_exon":in_exon,
                 "inclusion":inc, "delta":inc-base})
    df=pd.DataFrame(rows); df.to_csv(OUT+"/dmd_sweep.csv",index=False)
print("base inclusion",round(base,4),"n_masks",len(rows),"etep_span",et_span)

df=pd.DataFrame(rows)
fig,ax=plt.subplots(figsize=(12,4))
ax.axhline(0,color="k",lw=.6)
ax.axvspan(ES,EE,color="tab:blue",alpha=.08,label="exon 51 body")
ax.axvline(ES,color="tab:blue",ls="--",lw=.8); ax.axvline(EE,color="tab:blue",ls="--",lw=.8)
if et_span: ax.axvspan(et_span[0],et_span[1],color="tab:red",alpha=.30,label="eteplirsen target")
ax.plot(df.center, df.delta, lw=1.1, color="tab:green")
ax.scatter(df.center, df.delta, s=10, color="tab:green")
ax.set_title("AlphaGenome in-silico ASO-mask sweep across DMD exon 51\n(negative delta = masking this window is predicted to skip the exon)")
ax.set_xlabel("genomic position (chrX, GRCh38)"); ax.set_ylabel("delta exon-51 inclusion (mask - ref)")
ax.legend(loc="lower right", fontsize=8); fig.tight_layout()
fig.savefig(OUT+"/dmd_exon51_map.png",dpi=140)
print("SAVED figure")
