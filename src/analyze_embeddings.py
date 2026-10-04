import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import cosine_similarity

INPUT = "results/embeddings.npz"
OUTDIR = Path("results")
OUTDIR.mkdir(exist_ok=True)

# Load embeddings
data = np.load(INPUT)
embeddings = data["embeddings"]
names = data["names"].astype(str)

# Extract protein groups
families = np.array([
    name.split("|")[0].replace("_", " ")
    for name in names
])

# -------------------------
# PCA
# -------------------------
pca = PCA(n_components=2)
coords = pca.fit_transform(embeddings)

pca_df = pd.DataFrame({
    "PC1": coords[:, 0],
    "PC2": coords[:, 1],
    "Family": families,
    "Protein": names
})

pca_df.to_csv(OUTDIR / "pca_coordinates.csv", index=False)

plt.figure(figsize=(10, 7))

sns.scatterplot(
    data=pca_df,
    x="PC1",
    y="PC2",
    hue="Family",
    s=120
)

plt.title("ESM-2 Protein Embedding Space")
plt.xlabel(
    f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}% variance)"
)
plt.ylabel(
    f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}% variance)"
)

plt.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left"
)

plt.tight_layout()
plt.savefig(
    OUTDIR / "pca_embeddings.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# -------------------------
# Cosine similarity
# -------------------------
similarity = cosine_similarity(embeddings)

similarity_df = pd.DataFrame(
    similarity,
    index=names,
    columns=names
)

similarity_df.to_csv(
    OUTDIR / "similarity_matrix.csv"
)

plt.figure(figsize=(11, 9))

sns.heatmap(
    similarity,
    xticklabels=False,
    yticklabels=False
)

plt.title("Protein Embedding Cosine Similarity")
plt.xlabel("Protein")
plt.ylabel("Protein")

plt.tight_layout()

plt.savefig(
    OUTDIR / "similarity_heatmap.png",
    dpi=300
)

plt.close()

# -------------------------
# Within vs between family
# -------------------------
within = []
between = []

for i in range(len(names)):
    for j in range(i + 1, len(names)):

        if families[i] == families[j]:
            within.append(similarity[i, j])
        else:
            between.append(similarity[i, j])

within_mean = np.mean(within)
between_mean = np.mean(between)

# -------------------------
# Results
# -------------------------
print()
print("=" * 55)
print("PROTEIN SEQUENCE EMBEDDING ANALYSIS")
print("=" * 55)

print(f"Proteins: {len(names)}")
print(f"Embedding dimensions: {embeddings.shape[1]}")

print()
print("Protein groups:")
for family in sorted(set(families)):
    print(f"  {family}: {np.sum(families == family)}")

print()
print("Similarity analysis")
print("-------------------")
print(f"Mean within-group similarity : {within_mean:.4f}")
print(f"Mean between-group similarity: {between_mean:.4f}")

difference = within_mean - between_mean

print(f"Difference                  : {difference:.4f}")

print()
print("PCA")
print("---")
print(
    f"PC1 variance explained: "
    f"{pca.explained_variance_ratio_[0] * 100:.2f}%"
)

print(
    f"PC2 variance explained: "
    f"{pca.explained_variance_ratio_[1] * 100:.2f}%"
)

print(
    f"Total variance explained: "
    f"{pca.explained_variance_ratio_.sum() * 100:.2f}%"
)

print()
print("Results saved:")
print("  results/pca_embeddings.png")
print("  results/pca_coordinates.csv")
print("  results/similarity_heatmap.png")
print("  results/similarity_matrix.csv")

print()
print("Analysis complete!")
print("=" * 55)
