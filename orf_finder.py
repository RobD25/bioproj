from Bio.Seq import Seq

def find_orfs(dna_sequence, min_protein_len=20):
    """
    Scans a DNA sequence across all 6 reading frames to find open reading frames.
    """
    orfs = []
    
    # Store both forward and reverse complement strands
    strands = [
        ("+", dna_sequence),
        ("-", dna_sequence.reverse_complement())
    ]
    
    for strand_name, seq in strands:
        # Check all 3 frame offsets (0, 1, 2)
        for frame in range(3):
            # Translate from offset to end of sequence
            protein = seq[frame:].translate()
            
            # Find all methionine (M) start locations in this protein frame
            start_pos = 0
            while True:
                start = protein.find("M", start_pos)
                if start == -1:
                    break  # No more M found
                
                # Look for the first stop codon (*) following this start codon
                stop = protein.find("*", start)
                if stop != -1:
                    orf_protein = protein[start:stop]
                    
                    # Filter out tiny spurious peptides below min_protein_len
                    if len(orf_protein) >= min_protein_len:
                        orfs.append({
                            "strand": strand_name,
                            "frame": frame + 1,
                            "length_aa": len(orf_protein),
                            "protein": str(orf_protein)
                        })
                    start_pos = stop + 1
                else:
                    break
                    
    return orfs

# Test our function on human insulin mRNA sequence
# (We hardcode sequence for testing logic)
test_dna = Seq("AGCCCTCCAGGACAGGCTGCATCAGAAGAGCCATCAAGCAGATCACTGTCCTTCTGCCATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAGACGCAGCCCGCAGGCAGCCCCACACCCGCCGCCTCCTGCACCGAGAGAGATGGAATAAAGCCCTTGAACCAGC")

results = find_orfs(test_dna, min_protein_len=30)

print(f"Found {len(results)} ORF(s) with length >= 30 amino acids:\n")
for i, orf in enumerate(results, 1):
    print(f"ORF #{i} | Strand: {orf['strand']} | Frame: {orf['frame']} | Length: {orf['length_aa']} aa")
    print(f"Protein: {orf['protein']}\n")
