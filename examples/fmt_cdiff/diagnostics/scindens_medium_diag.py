import cobra, logging, warnings
warnings.filterwarnings("ignore"); logging.disable(logging.CRITICAL)
from pathlib import Path
GEMS=Path("examples/fmt_cdiff/gems")
avail={}
for l in open(GEMS/"western_gut_modelseed.csv"):
    if l.startswith(("#","metabolite")): continue
    p=l.strip().split(",")
    if len(p)>=4:
        try: mu=float(p[3]) if p[3] else 0.0
        except: mu=0.0
        if mu>0: avail[p[0]]=mu
NAME={"cpd00023":"Glu","cpd00033":"Gly","cpd00035":"Ala","cpd00039":"Lys","cpd00041":"Asp",
 "cpd00051":"Arg","cpd00053":"Gln","cpd00054":"Ser","cpd00060":"Met","cpd00065":"Trp",
 "cpd00066":"Phe","cpd00069":"Tyr","cpd00084":"Cys","cpd00107":"Leu","cpd00119":"His",
 "cpd00129":"Pro","cpd00156":"Val","cpd00161":"Thr","cpd00322":"Ile","cpd00027":"Glc",
 "cpd00013":"NH3","cpd00048":"SO4","cpd00007":"O2","cpd00029":"acetate","cpd00211":"butyrate"}
nm=lambda c: NAME.get(c.split("_")[0], c.split("_")[0])
def setmed(m, mode, extra_relax=None):
    for r in m.exchanges:
        met=list(r.metabolites)[0].id
        if mode=="open": r.lower_bound=-1000.0
        elif mode=="medium": r.lower_bound=-avail.get(met,0.0)
        elif mode=="relaxed": r.lower_bound=-10.0 if met in avail else 0.0
        if extra_relax and met==extra_relax: r.lower_bound=-10.0
m=cobra.io.read_sbml_model(str(GEMS/"C_scindens_ATCC35704.xml.gz"))
setmed(m,"medium");  g_med=m.slim_optimize() or 0
setmed(m,"relaxed"); g_rel=m.slim_optimize() or 0
setmed(m,"open");    g_open=m.slim_optimize() or 0
print(f"g_medium={g_med:.4f}  g_medium_relaxed(all avail -10)={g_rel:.3f}  g_open={g_open:.3f}")
print("-> if relaxed>>medium: the medium's tight max_uptake bounds are the limit\n")
# which single medium substrate, relaxed to -10, helps most?
print("single-substrate relaxation (medium bounds, one substrate -> -10):")
res=[]
for met in avail:
    setmed(m,"medium",extra_relax=met)
    g=m.slim_optimize() or 0
    if g>g_med*1.5+0.001: res.append((round(g,3), nm(met), round(avail[met],3)))
res.sort(reverse=True)
for g,name,mu in res[:12]: print(f"    +{name:10} (max_uptake in medium={mu}) -> growth {g}")
