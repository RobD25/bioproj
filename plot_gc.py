from Bio import Entrez, SeqIO
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt

Entrez.email = "your.email@example.com"

# Fetch Human Insulin sequence
handle = Entrez.efetch(db="nucleotide", id="NM_000207.3", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()

seq = record.seq
window_size = 30
positions = []
gc_values = []

# Sliding window calculation
for i in range(len(seq) - window_size + 1):
    subseq = seq[i : i + window_size]
    gc_pct = gc_fraction(subseq) * 100
    positions.append(i + window_size // 2)
    gc_values.append(gc_pct)

# Plotting the data
plt.figure(figsize=(10, 4))
plt.plot(positions, gc_values, color='teal', linewidth=2)
plt.axhline(y=50, color='gray', linestyle='--', alpha=0.7, label='50% GC baseline')

plt.title(f"Sliding Window GC Content ({record.id} - {record.description})")
plt.xlabel("Nucleotide Position (bp)")
plt.ylabel("GC Content (%)")
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.tight_layout()

# Save the plot to an image file
output_file = "gc_content_plot.png"
plt.savefig(output_file)
print(f"Plot successfully saved to {output_file}!")
