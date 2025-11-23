#!/usr/bin/env python3
"""
Python implementation of tBLASTn - search protein against DNA database
"""

def translate_dna(dna_seq):
    """Translate DNA to protein in all 6 frames"""
    
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    
    def translate_frame(seq):
        protein = []
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if 'N' in codon:
                protein.append('X')
            else:
                protein.append(codon_table.get(codon, 'X'))
        return ''.join(protein)
    
    # Reverse complement
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    rev_comp = ''.join(complement.get(base, 'N') for base in reversed(dna_seq))
    
    # Translate all 6 frames
    frames = []
    for frame in range(3):
        frames.append(translate_frame(dna_seq[frame:]))
        frames.append(translate_frame(rev_comp[frame:]))
    
    return frames

def read_fasta(filename):
    """Read FASTA file"""
    with open(filename, 'r') as f:
        header = ""
        sequence = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    yield header, ''.join(sequence)
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
        if sequence:
            yield header, ''.join(sequence)

def simple_alignment_score(query, target, match=1, mismatch=-1):
    """Quick local alignment score (simplified Smith-Waterman)"""
    if len(query) > len(target):
        return 0
    
    best_score = 0
    for i in range(len(target) - len(query) + 1):
        score = sum(match if q == target[i+j] else mismatch 
                   for j, q in enumerate(query))
        best_score = max(best_score, score)
    
    return best_score

def search_protein_in_dna(protein_query, dna_file, chunk_size=50000):
    """Search for protein in DNA by translating in all 6 frames"""
    
    print(f"🔍 Searching for chitinase in Metrosideros genome...")
    print(f"   Query length: {len(protein_query)} amino acids")
    print(f"   This will take a few minutes for 283 Mb genome...\n")
    
    hits = []
    total_bases = 0
    
    for header, dna_seq in read_fasta(dna_file):
        print(f"Scanning {header}...")
        total_bases += len(dna_seq)
        
        # Process in chunks to save memory
        for start in range(0, len(dna_seq), chunk_size):
            chunk = dna_seq[start:start + chunk_size + 1000]  # Overlap
            
            # Translate all 6 frames
            translations = translate_dna(chunk)
            
            # Score each frame
            for frame_num, translation in enumerate(translations):
                score = simple_alignment_score(protein_query, translation)
                
                if score > len(protein_query) * 0.3:  # 30% match threshold
                    strand = '+' if frame_num < 3 else '-'
                    frame = (frame_num % 3) + 1
                    
                    hits.append({
                        'chromosome': header,
                        'position': start,
                        'strand': strand,
                        'frame': frame,
                        'score': score
                    })
    
    print(f"\n✅ Scanned {total_bases:,} bp")
    return hits

def main():
    query_file = "eucalyptus_chitinase.fasta"
    genome_file = "metrosideros_chr1.fasta"
    
    # Read query
    print("📖 Reading query protein...")
    query_header, query_seq = next(read_fasta(query_file))
    print(f"   {query_header}")
    print(f"   Length: {len(query_seq)} aa\n")
    
    # Search
    hits = search_protein_in_dna(query_seq, genome_file)
    
    # Report results
    print(f"\n{'='*60}")
    print(f"🎯 Found {len(hits)} potential chitinase regions!")
    print(f"{'='*60}\n")
    
    if hits:
        # Sort by score
        hits.sort(key=lambda x: x['score'], reverse=True)
        
        print("Top 10 hits:")
        for i, hit in enumerate(hits[:10], 1):
            print(f"\n{i}. Chromosome: {hit['chromosome']}")
            print(f"   Position: {hit['position']:,} bp")
            print(f"   Strand: {hit['strand']}")
            print(f"   Frame: {hit['frame']}")
            print(f"   Score: {hit['score']}")
    else:
        print("❌ No significant hits found")
        print("   Try lowering the threshold or check your query sequence")

if __name__ == "__main__":
    main()
    