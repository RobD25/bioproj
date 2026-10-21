from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

print("DNA:", dna)
print("Reverse complement:", dna.reverse_complement())
print("RNA:", dna.transcribe())
print("Protein:", dna.translate(to_stop=True))
print("GC content:", round(gc_fraction(dna) * 100, 1), "%")
