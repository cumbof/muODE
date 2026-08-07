import urllib.request, urllib.parse
q = '(proteome:UP000000625) AND (protein_name:"ribosomal protein")'
url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(
    {"query": q, "format": "fasta"})
out = "/home/cumbof/muODE/data/clark2021/ecoli_riboprot.faa"
urllib.request.urlretrieve(url, out)
n = sum(1 for l in open(out) if l.startswith(">"))
print(f"fetched {n} E. coli ribosomal protein sequences -> {out}")
