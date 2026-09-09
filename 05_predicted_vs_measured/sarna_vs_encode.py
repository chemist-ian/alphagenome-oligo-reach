"""saRNA target selection: AlphaGenome PREDICTED tracks vs the free MEASURED ENCODE
tracks (HepG2), both against measured LHPP activation (Bi et al. 2024).
For a known human promoter in a common cell line the prediction reproduces the public
track, so it adds nothing a free ENCODE download does not already give you.
"""
import os, requests, numpy as np, pandas as pd, pyBigWig
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from alphagenome.data import genome
from alphagenome.models import dna_client

OUT=os.path.dirname(os.path.abspath(__file__))
CHROM="chr10"; TSS=124461823
LEADS=[("RAG7-133",-185,3.65),("RAG7-162",-214,2.98),("RAG7-694",-746,2.70),
       ("RAG7-892",-944,2.70),("RAG7-177",-229,2.48),("RAG7-132",-184,2.35),
       ("RAG7-178",-230,2.08),("RAG7-846",-898,1.65),("RAG7-139",-191,2.25),
       ("RAG7-707",-759,1.43)]
# measured ENCODE HepG2 bigWigs (GRCh38)
MEAS={"DNase":"https://www.encodeproject.org/files/ENCFF995ZMK/@@download/ENCFF995ZMK.bigWig",
      "ATAC": "https://www.encodeproject.org/files/ENCFF645OHB/@@download/ENCFF645OHB.bigWig",
      "CAGE": "https://www.encodeproject.org/files/ENCFF452THO/@@download/ENCFF452THO.bigWig"}
OT={"DNase":dna_client.OutputType.DNASE,"ATAC":dna_client.OutputType.ATAC,
    "CAGE":dna_client.OutputType.CAGE}
def open_bw(u): return pyBigWig.open(requests.head(u,allow_redirects=True,timeout=60).url)
bw={k:open_bw(u) for k,u in MEAS.items()}

client=dna_client.create(os.environ["ALPHAGENOME_API_KEY"])
L=dna_client.SEQUENCE_LENGTH_100KB; istart=TSS-L//2; iend=istart+L
iv=genome.Interval(chromosome=CHROM,start=istart,end=iend)
out=client.predict_interval(interval=iv,organism=dna_client.Organism.HOMO_SAPIENS,
    requested_outputs=list(OT.values()),ontology_terms=None)
pred={"DNase":np.asarray(out.dnase.values,np.float32).mean(1),
      "ATAC": np.asarray(out.atac.values,np.float32).mean(1),
      "CAGE": np.asarray(out.cage.values,np.float32).mean(1)}

def wmax_bw(b,g,ln=19,pad=10):
    v=b.stats(CHROM,g-pad,g+ln+pad,type="max")[0]; return float(v) if v is not None else 0.0
def wmax_pred(a,g,ln=19,pad=10):
    o=g-istart; return float(a[max(0,o-pad):o+ln+pad].max())
def sp(x,y):
    if np.ptp(x)==0 or np.ptp(y)==0: return None,None
    return spearmanr(x,y)

rows=[]
for name,rel,act in LEADS:
    g=TSS+rel; r={"name":name,"rel_TSS":rel,"activation":act}
    for k in MEAS:
        r[f"meas_{k}"]=wmax_bw(bw[k],g); r[f"pred_{k}"]=wmax_pred(pred[k],g)
    rows.append(r)
df=pd.DataFrame(rows).sort_values("activation",ascending=False)
df.to_csv(OUT+"/sarna_vs_encode.csv",index=False)
print(df.to_string(index=False))
print("\nSpearman vs LHPP activation (n=10):")
print(f"{'track':6} {'measured ENCODE':>18} {'AlphaGenome pred':>18} {'agreement':>10}")
for k in MEAS:
    rM,pM=sp(df[f'meas_{k}'].values,df.activation.values)
    rP,pP=sp(df[f'pred_{k}'].values,df.activation.values)
    rMP,_=sp(df[f'meas_{k}'].values,df[f'pred_{k}'].values)
    mtxt="flat (no signal upstream)" if rM is None else f"{rM:+.3f} (p={pM:.2f})"
    ptxt="flat" if rP is None else f"{rP:+.3f} (p={pP:.2f})"
    atxt="n/a" if rMP is None else f"rho={rMP:+.2f}"
    print(f"{k:6} {mtxt:>18} {ptxt:>18} {atxt:>10}")

lo,hi=TSS-1200,TSS+300; xs=np.arange(lo,hi)-TSS
fig,axes=plt.subplots(3,1,figsize=(12,8),sharex=True)
for ax,k in zip(axes,["DNase","ATAC","CAGE"]):
    m=np.nan_to_num(np.array(bw[k].values(CHROM,lo,hi))); p=pred[k][(lo-istart):(hi-istart)]
    ax.plot(xs,m/(m.max()+1e-9),color="black",lw=1.0,label=f"measured ENCODE HepG2 {k}")
    ax.plot(xs,p/(p.max()+1e-9),color="tab:blue",lw=1.0,alpha=.8,label=f"AlphaGenome predicted {k}")
    for _,rel,_ in LEADS: ax.axvline(rel,color="crimson",alpha=.30,lw=1)
    ax.axvline(0,color="k",lw=.6); ax.set_ylabel(f"norm {k}"); ax.legend(loc="upper left",fontsize=8)
axes[0].set_title("LHPP promoter: AlphaGenome prediction vs free measured ENCODE track (HepG2), saRNA sites in red")
axes[2].set_xlabel("position relative to LHPP TSS (bp)")
fig.tight_layout(); fig.savefig(OUT+"/sarna_vs_encode.png",dpi=140)
for b in bw.values(): b.close()
print("SAVED")
