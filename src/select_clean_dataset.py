from Bio import SeqIO
from collections import Counter

src = "data/proteins.fasta"
dst = "data/proteins_clean.fasta"

length_ranges = {
    "cytochrome_c": (70, 180),
    "ferredoxin": (70, 180),
    "lysozyme": (100, 220),
    "ubiquitin": (70, 120),
    "triosephosphate_isomerase": (220, 280),
}

records = list(SeqIO.parse(src, "fasta"))
selected = []

for family, (low, high) in length_ranges.items():

    candidates = [
        r for r in records
        if r.id.startswith(f"{family}|")
        and low <= len(r.seq) <= high
    ]

    print(f"{family}: {len(candidates)} valid candidates")

    if len(candidates) < 6:
        raise RuntimeError(
            f"Not enough valid sequences for {family}: {len(candidates)}"
        )

    selected.extend(candidates[:6])

with open(dst, "w") as f:
    for record in selected:
        f.write(f">{record.id}\n")
        f.write(f"{record.seq}\n")

print("\nSelected proteins:", len(selected))
print("Saved:", dst)

print("\nFamily counts:")
print(Counter(r.id.split("|")[0] for r in selected))
