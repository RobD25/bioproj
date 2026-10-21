import sys
from Bio import Entrez, SeqIO

# Set NCBI Email
Entrez.email = "your.email@example.com"

# Check if the user passed an accession ID in the terminal
if len(sys.argv) < 2:
    print("Error: Please provide an NCBI accession ID.")
    print("Usage: python3 gene_analyzer.py <ACCESSION_ID>")
    sys.exit(1)

# sys.argv[1] takes the ID you type after the script name in terminal
accession_id = sys.argv[1]

print(f"Fetching record for accession: {accession_id}...")

try:
    # 1. Fetch record from NCBI
    handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    # 2. Print basic summary
    print("\n================ SUCCESS ================")
    print(f"ID:          {record.id}")
    print(f"Description: {record.description}")
    print(f"Length:      {len(record.seq)} base pairs")

    # 3. Look for CDS and translate
    cds_found = False
    for feature in record.features:
        if feature.type == "CDS":
            cds_found = True
            cds_seq = feature.extract(record.seq)
            protein = cds_seq.translate(to_stop=True)
            print(f"\nCoding Sequence Length: {len(cds_seq)} bp")
            print(f"Translated Protein ({len(protein)} aa):")
            print(protein)
            break

    if not cds_found:
        print("\nNo CDS annotation found in this sequence record.")

except Exception as e:
    print(f"\nAn error occurred: {e}")
