from Bio import Entrez
from Bio import SeqIO

Entrez.email = "your.email@example.com"

# Fetch record
handle = Entrez.efetch(db="nucleotide", id="NM_000207.3", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()

# Loop through annotations to find the Coding Sequence (CDS)
print("--- Features Found in Record ---")
for feature in record.features:
    print(f"Feature type: {feature.type} at location {feature.location}")
    if feature.type == "CDS":
        # Extract CDS sequence directly using the feature's location slice
        cds_seq = feature.extract(record.seq)
        
        print("\n--- Coding Sequence (CDS) ---")
        print("CDS Length:", len(cds_seq), "bp")
        print("CDS Nucleotides:", cds_seq)
        
        print("\n--- Translated Protein ---")
        print("Protein from sequence:", cds_seq.translate())
        if "translation" in feature.qualifiers:
            print("NCBI Annotated Protein:", feature.qualifiers["translation"][0])
