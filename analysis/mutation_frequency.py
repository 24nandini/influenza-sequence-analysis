# HA 
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# Load aligned FASTA
alignment_file = r"C:\Users\24nan\Downloads\HA.fasta"
alignment = list(SeqIO.parse(alignment_file, "fasta"))

alignment_length = len(alignment[0].seq)
num_sequences = len(alignment)

mutation_freq = []

for i in range(alignment_length):
    column = [record.seq[i] for record in alignment if record.seq[i] != '-']  # ignore gaps
    unique_residues = set(column)
    freq = 1 - (column.count(max(unique_residues, key=column.count)) / len(column))  # fraction of sequences that differ from the most common residue
    mutation_freq.append(freq)

mutation_freq = np.array(mutation_freq)
plt.figure(figsize=(20,5))
plt.bar(range(1, alignment_length+1), mutation_freq, color='orange')
plt.xlabel("Alignment Position")
plt.ylabel("Mutation Frequency")
plt.title("Mutation Frequency Across Alignment in HA protein")
plt.ylim(0, 1)
plt.show()

# NA
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# Load aligned FASTA
alignment_file = r"C:\Users\24nan\Downloads\NA.fasta"
alignment = list(SeqIO.parse(alignment_file, "fasta"))

alignment_length = len(alignment[0].seq)
num_sequences = len(alignment)

mutation_freq = []

for i in range(alignment_length):
    column = [record.seq[i] for record in alignment if record.seq[i] != '-']  # ignore gaps
    unique_residues = set(column)
    freq = 1 - (column.count(max(unique_residues, key=column.count)) / len(column))  # fraction of sequences that differ from the most common residue
    mutation_freq.append(freq)

mutation_freq = np.array(mutation_freq)
plt.figure(figsize=(20,5))
plt.bar(range(1, alignment_length+1), mutation_freq, color='orange')
plt.xlabel("Alignment Position")
plt.ylabel("Mutation Frequency")
plt.title("Mutation Frequency Across Alignment in NA protein")
plt.ylim(0, 1)
plt.show()

# M1 
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# Load aligned FASTA
alignment_file = r"C:\Users\24nan\Downloads\M1.fasta"
alignment = list(SeqIO.parse(alignment_file, "fasta"))

alignment_length = len(alignment[0].seq)
num_sequences = len(alignment)

mutation_freq = []

for i in range(alignment_length):
    column = [record.seq[i] for record in alignment if record.seq[i] != '-']  # ignore gaps
    unique_residues = set(column)
    freq = 1 - (column.count(max(unique_residues, key=column.count)) / len(column))  # fraction of sequences that differ from the most common residue
    mutation_freq.append(freq)

mutation_freq = np.array(mutation_freq)
plt.figure(figsize=(20,5))
plt.bar(range(1, alignment_length+1), mutation_freq, color='orange')
plt.xlabel("Alignment Position")
plt.ylabel("Mutation Frequency")
plt.title("Mutation Frequency Across Alignment in M1 protein")
plt.ylim(0, 1)
plt.show()


# NP
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# Load aligned FASTA
alignment_file = r"C:\Users\24nan\Downloads\NP.fasta"
alignment = list(SeqIO.parse(alignment_file, "fasta"))

alignment_length = len(alignment[0].seq)
num_sequences = len(alignment)

mutation_freq = []

for i in range(alignment_length):
    column = [record.seq[i] for record in alignment if record.seq[i] != '-']  # ignore gaps
    unique_residues = set(column)
    freq = 1 - (column.count(max(unique_residues, key=column.count)) / len(column))  # fraction of sequences that differ from the most common residue
    mutation_freq.append(freq)

mutation_freq = np.array(mutation_freq)
plt.figure(figsize=(20,5))
plt.bar(range(1, alignment_length+1), mutation_freq, color='orange')
plt.xlabel("Alignment Position")
plt.ylabel("Mutation Frequency")
plt.title("Mutation Frequency Across Alignment in NP protein")
plt.ylim(0, 1)
plt.show()

# NS1 
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# Load aligned FASTA
alignment_file = r"C:\Users\24nan\Downloads\NS1.fasta"
alignment = list(SeqIO.parse(alignment_file, "fasta"))

alignment_length = len(alignment[0].seq)
num_sequences = len(alignment)

mutation_freq = []

for i in range(alignment_length):
    column = [record.seq[i] for record in alignment if record.seq[i] != '-']  # ignore gaps
    unique_residues = set(column)
    freq = 1 - (column.count(max(unique_residues, key=column.count)) / len(column))  # fraction of sequences that differ from the most common residue
    mutation_freq.append(freq)

mutation_freq = np.array(mutation_freq)
plt.figure(figsize=(20,5))
plt.bar(range(1, alignment_length+1), mutation_freq, color='orange')
plt.xlabel("Alignment Position")
plt.ylabel("Mutation Frequency")
plt.title("Mutation Frequency Across Alignment in NS1 protein")
plt.ylim(0, 1)
plt.show()




