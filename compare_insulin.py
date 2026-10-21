from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner

Entrez.email = "your.email@example.com"

# RefSeq Accessions for Insulin mRNA
species_ids = {
    "Human": "NM_000207.3",
    "Mouse (Ins2)": "NM_008387.5"
}

proteins = {}

# 1. Fetch CDS translation for each species
for name, acc_id in species_ids.items():
    handle = Entrez.efetch(db="nucleotide", id=acc_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    
    for feature in record.features:
        if feature.type == "CDS":
            cds_seq = feature.extract(record.seq)
            proteins[name] = cds_seq.translate(to_stop=True)
            break

print("--- Extracted Protein Sequences ---")
for species, prot in proteins.items():
    print(f"{species} ({len(prot)} aa): {prot}\n")

# 2. Perform Pairwise Sequence Alignment
aligner = PairwiseAligner()
aligner.mode = 'global'

alignment = aligner.align(proteins["Human"], proteins["Mouse (Ins2)"])[0]

print("--- Protein Sequence Alignment ---")
print(alignment)

# Calculate percent identity
matches = alignment.counts().identities
alignment_len = alignment.length
identity = (matches / alignment_len) * 100

print(f"Sequence Identity: {identity:.2f}%")
