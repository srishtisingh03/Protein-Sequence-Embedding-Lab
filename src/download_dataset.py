import requests
from pathlib import Path

families = {
    "cytochrome_c": "family:Cytochrome c AND reviewed:true",
    "ferredoxin": "family:Ferredoxin AND reviewed:true",
    "lysozyme": "family:Lysozyme AND reviewed:true",
    "ubiquitin": "family:Ubiquitin AND reviewed:true",
    "triosephosphate_isomerase": "family:Triosephosphate isomerase AND reviewed:true",
}

out = Path("data/proteins.fasta")

with out.open("w") as f:
    for family, query in families.items():

        params = {
            "query": query,
            "format": "fasta",
            "size": 6,
        }

        print(f"\nDownloading {family}...")

        response = requests.get(
            "https://rest.uniprot.org/uniprotkb/stream",
            params=params,
            timeout=60,
        )

        response.raise_for_status()

        records = response.text.strip().split(">")[1:]

        for record in records:
            lines = record.splitlines()
            header = lines[0]
            sequence = "".join(lines[1:])

            accession = header.split("|")[1]

            f.write(
                f">{family}|{accession}\n"
                f"{sequence}\n"
            )

            print(f"  {accession}")

print(f"\nDataset saved to: {out}")
