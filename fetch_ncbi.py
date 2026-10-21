from Bio import Entrez
from Bio import SeqIO

# 1. NCBI requires an email address for automated API requests
Entrez.email = "your.email@example.com"

# 2. Fetch the Human Insulin mRNA sequence (RefSeq accession: NM_000207.3)
print("Fetching Human Insulin sequence from NCBI...")
handle = Entrez.efetch(
    db="nucleotide",
    id="NM_000207.3",
    rettype="gb",      # GenBank format (includes annotations + sequence)
    retmode="text"
)

# 3. Parse the record using SeqIO
record = SeqIO.read(handle, "genbank")
handle.close()

# 4. Inspect the downloaded data
print("\n--- Record Summary ---")
print("ID:", record.id)
print("Name:", record.name)
print("Description:", record.description)
print("Sequence Length:", len(record.seq), "bp")

# Print the first 100 nucleotides
print("\nFirst 100 bases:")
print(record.seq[:100])
