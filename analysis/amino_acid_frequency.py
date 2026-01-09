"""
Amino acid frequency analysis for influenza A protein sequences.
"""

import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

df = pd.read_excel("data/processed/FASTA_with_length.xlsx")

# Combine all sequences into one long string
all_sequences = "".join(df["Sequence"].tolist())

# Count residues
counts = Counter(all_sequences)

# Standard 20 amino acids
amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
frequencies = [counts[aa] for aa in amino_acids]

# Plot
plt.figure(figsize=(12, 6))
plt.bar(amino_acids, frequencies, edgecolor="black")
plt.xlabel("Amino Acid")
plt.ylabel("Frequency")
plt.title("Amino Acid Frequency – Influenza A Sequences")
plt.show()
