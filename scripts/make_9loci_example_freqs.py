"""Regenerate data/freqs_9loci/CAU.freqs.gz, the 9-locus example frequency file.

Run from the repo root:  python scripts/make_9loci_example_freqs.py

The data is illustrative, not a reference dataset - it is built from the
5-locus CAU sample that ships with the repo.

The 5-locus sample carries A~C~B~DRB1~DQB1. The four extra loci are filled in
from the linkage they actually have with what is already there: DRB3/4/5 and
DQA1 follow DRB1 and DQB1 closely enough to be derived, while DP sits across a
recombination hotspot and is spread over a few common DPA1~DPB1 haplotypes.
"""
import gzip

SRC = "data/freqs/CAU.freqs.gz"
DST = "data/freqs_9loci/CAU.freqs.gz"
N_SOURCE_HAPS = 200

# DRB1 group -> the DRB3/4/5 allele carried on that haplotype. Groups absent
# here (DRB1*01, *08, *10) carry no second DRB locus at all, so they are left
# out of the example rather than modelled with a null placeholder.
DRB345 = {
    "03": "DRB3*01:01", "11": "DRB3*02:02", "12": "DRB3*02:02",
    "13": "DRB3*01:01", "14": "DRB3*02:02",
    "04": "DRB4*01:03", "07": "DRB4*01:03", "09": "DRB4*01:03",
    "15": "DRB5*01:01", "16": "DRB5*02:02",
}
# A few DRB1 alleles break the group default.
DRB345_EXACT = {
    "DRB1*13:02": "DRB3*03:01",
    "DRB1*15:02": "DRB5*01:02",
    "DRB1*16:01": "DRB5*02:02",
}

DQA1 = {
    "DQB1*02:01": "DQA1*05:01", "DQB1*03:01": "DQA1*05:05",
    "DQB1*03:02": "DQA1*03:01", "DQB1*03:03": "DQA1*03:02",
    "DQB1*03:04": "DQA1*04:01", "DQB1*03:05": "DQA1*05:05",
    "DQB1*04:02": "DQA1*04:01", "DQB1*05:01": "DQA1*01:01",
    "DQB1*05:02": "DQA1*01:02", "DQB1*05:03": "DQA1*01:04",
    "DQB1*05:04": "DQA1*01:02", "DQB1*06:01": "DQA1*01:03",
    "DQB1*06:02": "DQA1*01:02", "DQB1*06:03": "DQA1*01:03",
    "DQB1*06:04": "DQA1*01:02", "DQB1*06:09": "DQA1*01:02",
}

# DPA1~DPB1 haplotype and the share of each 5-locus haplotype it takes.
DP = [
    ("DPA1*01:03", "DPB1*04:01", 0.45),
    ("DPA1*01:03", "DPB1*02:01", 0.25),
    ("DPA1*01:03", "DPB1*04:02", 0.18),
    ("DPA1*02:01", "DPB1*01:01", 0.12),
]

rows = []
with gzip.open(SRC, "rt") as f:
    for line in f:
        haplotype, count, freq = line.strip().split(",")
        if haplotype == "Haplo":
            continue
        a, c, b, drb1, dqb1 = haplotype.split("~")
        group = drb1.split("*")[1].split(":")[0]
        if group not in DRB345:
            continue
        rows.append((float(freq), a, b, c, drb1, dqb1))

rows.sort(reverse=True)
rows = rows[:N_SOURCE_HAPS]

out = []
for freq, a, b, c, drb1, dqb1 in rows:
    drbx = DRB345_EXACT.get(drb1, DRB345[drb1.split("*")[1].split(":")[0]])
    dqa1 = DQA1[dqb1]
    for dpa1, dpb1, share in DP:
        hap = "~".join([a, b, c, dpa1, dpb1, dqa1, dqb1, drb1, drbx])
        out.append((hap, freq * share))

with gzip.open(DST, "wt") as f:
    for hap, freq in out:
        f.write("{},1,{:.3e}\n".format(hap, freq))

print("haplotypes:", len(out))
print("min freq:  ", min(f for _, f in out))
print("max freq:  ", max(f for _, f in out))
print("trim threshold must stay below:", min(f for _, f in out) * len(out))
