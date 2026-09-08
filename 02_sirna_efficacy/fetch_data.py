"""Download the siRNA inputs (siRNADiscovery HUVK compilation) next to this script.
Run once before ag_sirna.py. Source: Long et al., Brief Bioinform 2024.
"""
import os, requests

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "https://raw.githubusercontent.com/BertramLoong/siRNADiscovery/main/Data/Dataset_HUVK"
FILES = {
    "eff.csv": "siRNA_mRNA_Efficacy.csv",
    "mRNA_HUVK.fas": "mRNA_HUVK.fas",
    "siRNA_HUVK.fas": "siRNA_HUVK.fas",
}
for local, remote in FILES.items():
    r = requests.get(f"{BASE}/{remote}", timeout=60)
    r.raise_for_status()
    open(os.path.join(HERE, local), "wb").write(r.content)
    print("saved", local, len(r.content), "bytes")
print("done")
