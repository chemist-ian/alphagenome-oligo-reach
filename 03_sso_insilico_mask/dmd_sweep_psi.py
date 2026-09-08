import os, requests, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from alphagenome.data import genome
from alphagenome.models import dna_client
OUT=os.path.dirname(os.path.abspath(__file__))
client=dna_client.create(os.environ["ALPHAGENOME_API_KEY"])
CHROM="chrX"; ES=31773960; EE=31774192
L=dna_client.SEQUENCE_LENGTH_100KB
center=(ES+EE)//2; istart=center-L//2; iend=istart+L
iv=genome.Interval(chromosome=CHROM,start=istart,end=iend)
TRANS=str.maketrans("ACGT","GTAC"); comp=str.maketrans("ACGT","TGCA"); rc=lambda s:s.translate(comp)[::-1]
def ens(a,b):
    r=requests.get(f"https://rest.ensembl.org/sequence/region/human/X:{a}..{b}",
       headers={"Content-Type":"text/plain"},timeout=120); r.raise_for_status(); return r.text.strip().upper()
ref=ens(istart+1,iend); assert len(ref)==L

INCL=[(31729748,31773959),(31774192,31819974)]; SKIP=(31729748,31819974)
def find_j(juncs, tgt):
    # predict_sequence returns interval-RELATIVE coords; convert to genomic with istart
    for i,j in enumerate(juncs):
        if abs((j.start+istart)-tgt[0])<=2 and abs((j.end+istart)-tgt[1])<=2: return i
    return None
def psi(seq):
    out=client.predict_sequence(sequence=seq,organism=dna_client.Organism.HOMO_SAPIENS,
        requested_outputs=[dna_client.OutputType.SPLICE_JUNCTIONS],ontology_terms=None,interval=iv)
    sj=out.splice_junctions; v=np.asarray(sj.values); j=list(sj.junctions)
    tot=np.nan_to_num(v).sum(axis=1)  # sum across tracks per junction
    i1=find_j(j,INCL[0]); i2=find_j(j,INCL[1]); isk=find_j(j,SKIP)
    incl=(tot[i1] if i1 is not None else 0)+(tot[i2] if i2 is not None else 0)
    skip=(tot[isk] if isk is not None else 0)
    return incl/(incl+skip+1e-9)

base=psi(ref)
etep="CTCCAACATCAAGGAAGATGGCATTTCTAG"; et_span=None
for cand in [etep, rc(etep)]:
    o=ref.find(cand)
    if o>=0: et_span=(istart+o,istart+o+len(cand)); break

WIN=20; STEP=4; lo,hi=ES-60,EE+60; rows=[]
for gpos in range(lo,hi-WIN,STEP):
    o1=gpos-istart; o2=o1+WIN
    mut=ref[:o1]+ref[o1:o2].translate(TRANS)+ref[o2:]
    p=psi(mut); ctr=gpos+WIN//2
    rows.append({"center":ctr,"in_exon":ES<=ctr<=EE,"psi":p,"delta":p-base})
    pd.DataFrame(rows).to_csv(OUT+"/dmd_sweep_psi.csv",index=False)
print("base PSI",round(base,4),"n",len(rows),"etep",et_span)
df=pd.DataFrame(rows)
fig,ax=plt.subplots(figsize=(12,4)); ax.axhline(0,color="k",lw=.6)
ax.axvspan(ES,EE,color="tab:blue",alpha=.08,label="exon 51 body")
ax.axvline(ES,color="tab:blue",ls="--",lw=.8); ax.axvline(EE,color="tab:blue",ls="--",lw=.8)
if et_span: ax.axvspan(*et_span,color="tab:red",alpha=.30,label="eteplirsen target")
ax.plot(df.center,df.delta,lw=1.1,color="darkorange"); ax.scatter(df.center,df.delta,s=10,color="darkorange")
ax.set_title("DMD exon 51: in-silico ASO-mask sweep, PSI (inclusion/skip junction) readout\n(negative = masking predicted to skip exon 51)")
ax.set_xlabel("genomic position (chrX, GRCh38)"); ax.set_ylabel("delta PSI (mask - ref)")
ax.legend(loc="lower right",fontsize=8); fig.tight_layout(); fig.savefig(OUT+"/dmd_exon51_psi.png",dpi=140)
# region stats
def st(m,lab):
    s=df[m]; print(f"{lab:30} mean={s.delta.mean():+.3f} min={s.delta.min():+.3f} n={len(s)}")
st((df.center>=ES-30)&(df.center<=ES+15),"donor 5'SS flank")
st((df.center>=EE-15)&(df.center<=EE+30),"acceptor 3'SS flank")
st((df.center>=ES+25)&(df.center<=EE-25),"exon INTERIOR")
if et_span: st((df.center>=et_span[0])&(df.center<=et_span[1]),"eteplirsen span")
print("PSI_SWEEP_DONE")
