# HA 
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

alignment_file = r"C:\Users\24nan\Downloads\HA.fasta"
alignment = list(SeqIO.parse(alignment_file, "fasta"))

alignment_length = len(alignment[0].seq)
num_sequences = len(alignment)

data = np.zeros((num_sequences, alignment_length))

for col in range(alignment_length):
    column_residues = [record.seq[col] for record in alignment if record.seq[col] != '-']
    if len(set(column_residues)) == 1 and column_residues:
        data[:, col] = 1  # conserved
    else:
        data[:, col] = 0  # variable

plt.figure(figsize=(20, 4))
plt.imshow(data, aspect='auto', cmap='plasma', interpolation='none')
plt.colorbar(label='Conserved (1) / Variable (0)')
plt.xlabel("Alignment Position")
plt.ylabel("Sequences")
plt.title("Conserved vs Variable Positions in Alignment for HA protein")
plt.show()

# NA 
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

alignment_file = r"C:\Users\24nan\Downloads\NA.fasta"
alignment = list(SeqIO.parse(alignment_file, "fasta"))

alignment_length = len(alignment[0].seq)
num_sequences = len(alignment)

data = np.zeros((num_sequences, alignment_length))

for col in range(alignment_length):
    column_residues = [record.seq[col] for record in alignment if record.seq[col] != '-']
    if len(set(column_residues)) == 1 and column_residues:
        data[:, col] = 1  # conserved
    else:
        data[:, col] = 0  # variable
        
plt.figure(figsize=(20, 4))
plt.imshow(data, aspect='auto', cmap='plasma', interpolation='none')
plt.colorbar(label='Conserved (1) / Variable (0)')
plt.xlabel("Alignment Position")
plt.ylabel("Sequences")
plt.title("Conserved vs Variable Positions in Alignment for NA protein")
plt.show()

# NS1 
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

alignment_file = r"C:\Users\24nan\Downloads\NS1.fasta"
alignment = list(SeqIO.parse(alignment_file, "fasta"))

alignment_length = len(alignment[0].seq)
num_sequences = len(alignment)

data = np.zeros((num_sequences, alignment_length))

for col in range(alignment_length):
    column_residues = [record.seq[col] for record in alignment if record.seq[col] != '-']
    if len(set(column_residues)) == 1 and column_residues:
        data[:, col] = 1  # conserved
    else:
        data[:, col] = 0  # variable

plt.figure(figsize=(20, 4))
plt.imshow(data, aspect='auto', cmap='plasma', interpolation='none')
plt.colorbar(label='Conserved (1) / Variable (0)')
plt.xlabel("Alignment Position")
plt.ylabel("Sequences")
plt.title("Conserved vs Variable Positions in Alignment for NS1 ptotein")
plt.show()

# M1 
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

alignment_file = r"C:\Users\24nan\Downloads\M1.fasta"
alignment = list(SeqIO.parse(alignment_file, "fasta"))

alignment_length = len(alignment[0].seq)
num_sequences = len(alignment)

data = np.zeros((num_sequences, alignment_length))

for col in range(alignment_length):
    column_residues = [record.seq[col] for record in alignment if record.seq[col] != '-']
    if len(set(column_residues)) == 1 and column_residues:
        data[:, col] = 1  # conserved
    else:
        data[:, col] = 0  # variable

plt.figure(figsize=(20, 4))
plt.imshow(data, aspect='auto', cmap='plasma', interpolation='none')
plt.colorbar(label='Conserved (1) / Variable (0)')
plt.xlabel("Alignment Position")
plt.ylabel("Sequences")
plt.title("Conserved vs Variable Positions in Alignment for M1 protein")
plt.show()

# M2 
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

alignment_file = r"C:\Users\24nan\Downloads\M2.fasta"
alignment = list(SeqIO.parse(alignment_file, "fasta"))

alignment_length = len(alignment[0].seq)
num_sequences = len(alignment)

data = np.zeros((num_sequences, alignment_length))

for col in range(alignment_length):
    column_residues = [record.seq[col] for record in alignment if record.seq[col] != '-']
    if len(set(column_residues)) == 1 and column_residues:
        data[:, col] = 1  # conserved
    else:
        data[:, col] = 0  # variable

plt.figure(figsize=(20, 4))
plt.imshow(data, aspect='auto', cmap='plasma', interpolation='none')
plt.colorbar(label='Conserved (1) / Variable (0)')
plt.xlabel("Alignment Position")
plt.ylabel("Sequences")
plt.title("Conserved vs Variable Positions in Alignment for M2 protein")
plt.show()

# NP 
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

alignment_file = r"C:\Users\24nan\Downloads\NP.fasta"
alignment = list(SeqIO.parse(alignment_file, "fasta"))

alignment_length = len(alignment[0].seq)
num_sequences = len(alignment)

data = np.zeros((num_sequences, alignment_length))

for col in range(alignment_length):
    column_residues = [record.seq[col] for record in alignment if record.seq[col] != '-']
    if len(set(column_residues)) == 1 and column_residues:
        data[:, col] = 1  # conserved
    else:
        data[:, col] = 0  # variable

plt.figure(figsize=(20, 4))
plt.imshow(data, aspect='auto', cmap='plasma', interpolation='none')
plt.colorbar(label='Conserved (1) / Variable (0)')
plt.xlabel("Alignment Position")
plt.ylabel("Sequences")
plt.title("Conserved vs Variable Positions in Alignment for NP protein")
plt.show()
