import io, pyrodigal, pyhmmer
from pathlib import Path
from pyhmmer.easel import Alphabet, SequenceFile, TextSequence, DigitalSequenceBlock

GEN = Path("/mnt/isilon/w_gmi/blanked2lab/cumbof/git/muODE/data/clark2021/genomes")
OUT = Path("/home/cumbof/muODE/data/clark2021/cds"); OUT.mkdir(parents=True, exist_ok=True)
REF = "/home/cumbof/muODE/data/clark2021/ecoli_riboprot.faa"
abc = Alphabet.amino()
with SequenceFile(REF, digital=True, alphabet=abc) as sf:
    refs = sf.read_block()

def read_fna(p):
    seqs=[]; name=None; buf=[]
    for line in open(p):
        if line.startswith(">"):
            if name: seqs.append((name,"".join(buf)))
            name=line[1:].split()[0]; buf=[]
        else: buf.append(line.strip())
    if name: seqs.append((name,"".join(buf)))
    return seqs

def parse_fasta(text):
    recs=[]; gid=None; buf=[]
    for line in text.splitlines():
        if line.startswith(">"):
            if gid is not None: recs.append((gid,"".join(buf)))
            gid=line[1:].split()[0]; buf=[]
        else: buf.append(line.strip())
    if gid is not None: recs.append((gid,"".join(buf)))
    return recs

for fna in sorted(GEN.glob("*.fna")):
    code=fna.stem
    seqs=read_fna(fna)
    gf=pyrodigal.GeneFinder(meta=False); gf.train(*[s for _,s in seqs])
    nt_io=io.StringIO(); aa_io=io.StringIO()
    for sid,seq in seqs:
        g=gf.find_genes(seq); g.write_genes(nt_io, sid); g.write_translations(aa_io, sid)
    nt=parse_fasta(nt_io.getvalue()); aa=parse_fasta(aa_io.getvalue())
    # digital block of proteins (strip trailing *)
    prots=[TextSequence(name=gid.encode(), sequence=s.rstrip("*")).digitize(abc)
           for gid,s in aa if s.rstrip("*")]
    block=DigitalSequenceBlock(abc, prots)
    ribo=set()
    for hits in pyhmmer.phmmer(refs, block, E=1e-8, cpus=8):
        for h in hits:
            if h.included:
                _n=h.name; ribo.add(_n.decode() if isinstance(_n,bytes) else _n)
    with open(OUT/f"{code}.ffn","w") as fh:
        for gid,s in nt:
            tag=" ribosomal_protein" if gid in ribo else ""
            fh.write(f">{gid}{tag}\n{s}\n")
    print(f"{code}: {len(nt)} genes, {len(ribo)} ribosomal", flush=True)
print("ANNOTATE_DONE", flush=True)
