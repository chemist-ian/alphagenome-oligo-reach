import os, requests, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr
from alphagenome.data import genome
from alphagenome.models import dna_client
OUT=os.path.dirname(os.path.abspath(__file__))
client=dna_client.create(os.environ["ALPHAGENOME_API_KEY"])

CHROM="chr10"; STRAND=1; TSS=124461823
# lead saRNAs: sense sequence (U->T), and activation = mean(Huh7,HepG2) read from Fig 3
LEADS={
 "RAG7-133":("GCTCTTTGTCCGCTGATCT",3.65),"RAG7-162":("TTCTTAGGGACTTGTTTTC",2.98),
 "RAG7-694":("AGGTCCTATGCATCCTCAT",2.70),"RAG7-892":("TGTTGGACCAGAAGTAAAG",2.70),
 "RAG7-177":("ATTTGCCTTTGACCTTTCT",2.48),"RAG7-132":("CTCTTTGTCCGCTGATCTC",2.35),
 "RAG7-178":("AATTTGCCTTTGACCTTTC",2.08),"RAG7-846":("AAGGTTCCGAGGGGCCATT",1.65),
 "RAG7-139":("TTTCCTGCTCTTTGTCCGC",2.25),"RAG7-707":("TTCTTCTCAGCCCAGGTCC",1.43)}
comp=str.maketrans("ACGT","TGCA"); rc=lambda s:s.translate(comp)[::-1]

def ens(chrom,a,b):
    r=requests.get(f"https://rest.ensembl.org/sequence/region/human/{chrom.replace('chr','')}:{a}..{b}",
       headers={"Content-Type":"text/plain"},timeout=120); r.raise_for_status(); return r.text.strip().upper()

# AlphaGenome interval centered on TSS
L=dna_client.SEQUENCE_LENGTH_100KB
istart=TSS-L//2; iend=istart+L
iv=genome.Interval(chromosome=CHROM,start=istart,end=iend)
prom=ens(CHROM,istart+1,iend)  # index o -> genomic istart+o

# map each lead to genomic offset
sites={}
for name,(seq,act) in LEADS.items():
    o=prom.find(seq)
    strand_hit="+"
    if o<0: o=prom.find(rc(seq)); strand_hit="-"
    if o<0:
        # fallback: position label
        N=int(name.split("-")[1]); g=TSS-N; o=g-istart; strand_hit="label"
    g=istart+o
    rel=g-TSS  # relative to TSS (negative=upstream)
    sites[name]=dict(offset=o,genomic=g,rel=rel,act=act,strand=strand_hit,ln=len(seq))

out=client.predict_interval(interval=iv,organism=dna_client.Organism.HOMO_SAPIENS,
    requested_outputs=[dna_client.OutputType.CAGE,dna_client.OutputType.PROCAP,
                       dna_client.OutputType.DNASE,dna_client.OutputType.ATAC],ontology_terms=None)
def red(td):
    return np.asarray(td.values,dtype=np.float32).mean(axis=1)
sig={"CAGE":red(out.cage),"PROCAP":red(out.procap),"DNASE":red(out.dnase),"ATAC":red(out.atac)}

rows=[]
for name,d in sites.items():
    o=d["offset"]; w=slice(max(0,o-10),min(L,o+d["ln"]+10))
    row=dict(name=name,rel_TSS=d["rel"],activation=d["act"],strand=d["strand"])
    for k,s in sig.items(): row[k]=float(s[w].max())  # peak signal over target window+/-10
    rows.append(row)
df=pd.DataFrame(rows).sort_values("activation",ascending=False)
df.to_csv(OUT+"/sarna_lhpp_features.csv",index=False)
print(df.to_string(index=False))
print("\n== Spearman (feature vs activation), n=%d ==" % len(df))
for k in ["CAGE","PROCAP","DNASE","ATAC"]:
    rho,p=spearmanr(df[k],df.activation); print(f"  {k:7} rho={rho:+.3f} p={p:.3f}")

# figure: promoter regulatory profile + lead positions
xs=np.arange(istart,iend)-TSS
sel=(xs>=-1200)&(xs<=300)
fig,ax=plt.subplots(figsize=(12,4))
for k,c in [("DNASE","tab:blue"),("CAGE","tab:green")]:
    y=sig[k][sel]; y=y/ (y.max()+1e-9)
    ax.plot(xs[sel],y,lw=1.0,color=c,alpha=.8,label=f"AlphaGenome {k} (norm)")
for name,d in sites.items():
    ax.axvline(d["rel"],color="crimson",alpha=.35,lw=1)
    ax.scatter([d["rel"]],[d["act"]/df.activation.max()],s=25+30*(d["act"]-1.4),color="crimson",zorder=5)
    ax.annotate(name.replace("RAG7-",""),(d["rel"],d["act"]/df.activation.max()),fontsize=7,rotation=90,ha="center",va="bottom")
ax.axvline(0,color="k",lw=.8); ax.text(5,0.95,"TSS",fontsize=8)
ax.set_xlabel("position relative to LHPP TSS (bp)"); ax.set_ylabel("normalized signal / activation")
ax.set_title("LHPP saRNA leads (Bi et al. 2024) vs AlphaGenome promoter regulatory signal")
ax.legend(loc="upper left",fontsize=8); fig.tight_layout(); fig.savefig(OUT+"/sarna_lhpp.png",dpi=140)
print("\nSAVED")
