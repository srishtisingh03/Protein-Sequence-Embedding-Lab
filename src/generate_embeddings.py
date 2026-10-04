import torch
import numpy as np
from pathlib import Path
from Bio import SeqIO
from transformers import AutoTokenizer, AutoModel

MODEL_NAME = "facebook/esm2_t6_8M_UR50D"
FASTA = "data/proteins_clean.fasta"
OUTPUT = "results/embeddings.npz"

device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
print(f"Using device: {device}")

tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
model = AutoModel.from_pretrained(MODEL_NAME)
model = model.to(device)
model.eval()

names = []
embeddings = []

records = list(SeqIO.parse(FASTA, "fasta"))

for i, record in enumerate(records, 1):
    sequence = str(record.seq)

    print(f"[{i}/{len(records)}] {record.id} ({len(sequence)} aa)")

    inputs = tokenizer(
        sequence,
        return_tensors="pt",
        truncation=True,
        max_length=1022,
    )

    inputs = {k: v.to(device) for k, v in inputs.items()}

    with torch.no_grad():
        outputs = model(**inputs)

    # Remove special tokens and mean-pool residue embeddings
    attention_mask = inputs["attention_mask"]
    hidden = outputs.last_hidden_state

    mask = attention_mask.unsqueeze(-1).expand(hidden.size()).float()

    embedding = (hidden * mask).sum(dim=1) / mask.sum(dim=1)

    embeddings.append(embedding.cpu().numpy()[0])
    names.append(record.id)

embeddings = np.vstack(embeddings)

Path("results").mkdir(exist_ok=True)

np.savez(
    OUTPUT,
    embeddings=embeddings,
    names=np.array(names),
)

print("\nEmbedding generation complete!")
print("Shape:", embeddings.shape)
print("Saved:", OUTPUT)
