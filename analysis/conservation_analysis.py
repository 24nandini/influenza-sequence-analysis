# HA Serotype
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# Path to aligned FASTA
alignment_file = r"C:\Users\24nan\Downloads\HA.fasta"
alignment = list(SeqIO.parse(alignment_file, "fasta"))

alignment_length = len(alignment[0].seq)
num_sequences = len(alignment)

# Build a 2D array: rows=sequences, columns=positions
data = np.zeros((num_sequences, alignment_length))

for col in range(alignment_length):
    column_residues = [record.seq[col] for record in alignment if record.seq[col] != '-']
    # Determine if position is conserved
    if len(set(column_residues)) == 1 and column_residues:
        data[:, col] = 1  # conserved
    else:
        data[:, col] = 0  # variable

# Plot
plt.figure(figsize=(20, 4))
plt.imshow(data, aspect='auto', cmap='plasma', interpolation='none')
plt.colorbar(label='Conserved (1) / Variable (0)')
plt.xlabel("Alignment Position")
plt.ylabel("Sequences")
plt.title("Conserved vs Variable Positions in Alignment for HA Serotype")
plt.show()

# NA Serotype
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# Path to aligned FASTA
alignment_file = r"C:\Users\24nan\Downloads\NA.fasta"
alignment = list(SeqIO.parse(alignment_file, "fasta"))

alignment_length = len(alignment[0].seq)
num_sequences = len(alignment)

# Build a 2D array: rows=sequences, columns=positions
data = np.zeros((num_sequences, alignment_length))

for col in range(alignment_length):
    column_residues = [record.seq[col] for record in alignment if record.seq[col] != '-']
    # Determine if position is conserved
    if len(set(column_residues)) == 1 and column_residues:
        data[:, col] = 1  # conserved
    else:
        data[:, col] = 0  # variable

# Plot
plt.figure(figsize=(20, 4))
plt.imshow(data, aspect='auto', cmap='plasma', interpolation='none')
plt.colorbar(label='Conserved (1) / Variable (0)')
plt.xlabel("Alignment Position")
plt.ylabel("Sequences")
plt.title("Conserved vs Variable Positions in Alignment for NA Serotype")
plt.show()

# NS1 Serotype
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# Path to aligned FASTA
alignment_file = r"C:\Users\24nan\Downloads\NS1.fasta"
alignment = list(SeqIO.parse(alignment_file, "fasta"))

alignment_length = len(alignment[0].seq)
num_sequences = len(alignment)

# Build a 2D array: rows=sequences, columns=positions
data = np.zeros((num_sequences, alignment_length))

for col in range(alignment_length):
    column_residues = [record.seq[col] for record in alignment if record.seq[col] != '-']
    # Determine if position is conserved
    if len(set(column_residues)) == 1 and column_residues:
        data[:, col] = 1  # conserved
    else:
        data[:, col] = 0  # variable

# Plot
plt.figure(figsize=(20, 4))
plt.imshow(data, aspect='auto', cmap='plasma', interpolation='none')
plt.colorbar(label='Conserved (1) / Variable (0)')
plt.xlabel("Alignment Position")
plt.ylabel("Sequences")
plt.title("Conserved vs Variable Positions in Alignment for NS1 Serotype")
plt.show()

# M1 Serotype
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# Path to aligned FASTA
alignment_file = r"C:\Users\24nan\Downloads\M1.fasta"
alignment = list(SeqIO.parse(alignment_file, "fasta"))

alignment_length = len(alignment[0].seq)
num_sequences = len(alignment)

# Build a 2D array: rows=sequences, columns=positions
data = np.zeros((num_sequences, alignment_length))

for col in range(alignment_length):
    column_residues = [record.seq[col] for record in alignment if record.seq[col] != '-']
    # Determine if position is conserved
    if len(set(column_residues)) == 1 and column_residues:
        data[:, col] = 1  # conserved
    else:
        data[:, col] = 0  # variable

# Plot
plt.figure(figsize=(20, 4))
plt.imshow(data, aspect='auto', cmap='plasma', interpolation='none')
plt.colorbar(label='Conserved (1) / Variable (0)')
plt.xlabel("Alignment Position")
plt.ylabel("Sequences")
plt.title("Conserved vs Variable Positions in Alignment for M1 Serotype")
plt.show()

# M2 Serotype
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# Path to aligned FASTA
alignment_file = r"C:\Users\24nan\Downloads\M2.fasta"
alignment = list(SeqIO.parse(alignment_file, "fasta"))

alignment_length = len(alignment[0].seq)
num_sequences = len(alignment)

# Build a 2D array: rows=sequences, columns=positions
data = np.zeros((num_sequences, alignment_length))

for col in range(alignment_length):
    column_residues = [record.seq[col] for record in alignment if record.seq[col] != '-']
    # Determine if position is conserved
    if len(set(column_residues)) == 1 and column_residues:
        data[:, col] = 1  # conserved
    else:
        data[:, col] = 0  # variable

# Plot
plt.figure(figsize=(20, 4))
plt.imshow(data, aspect='auto', cmap='plasma', interpolation='none')
plt.colorbar(label='Conserved (1) / Variable (0)')
plt.xlabel("Alignment Position")
plt.ylabel("Sequences")
plt.title("Conserved vs Variable Positions in Alignment for M2 Serotype")
plt.show()

# NP Serotype
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# Path to aligned FASTA
alignment_file = r"C:\Users\24nan\Downloads\NP.fasta"
alignment = list(SeqIO.parse(alignment_file, "fasta"))

alignment_length = len(alignment[0].seq)
num_sequences = len(alignment)

# Build a 2D array: rows=sequences, columns=positions
data = np.zeros((num_sequences, alignment_length))

for col in range(alignment_length):
    column_residues = [record.seq[col] for record in alignment if record.seq[col] != '-']
    # Determine if position is conserved
    if len(set(column_residues)) == 1 and column_residues:
        data[:, col] = 1  # conserved
    else:
        data[:, col] = 0  # variable

# Plot
plt.figure(figsize=(20, 4))
plt.imshow(data, aspect='auto', cmap='plasma', interpolation='none')
plt.colorbar(label='Conserved (1) / Variable (0)')
plt.xlabel("Alignment Position")
plt.ylabel("Sequences")
plt.title("Conserved vs Variable Positions in Alignment for NP Serotype")
plt.show()
