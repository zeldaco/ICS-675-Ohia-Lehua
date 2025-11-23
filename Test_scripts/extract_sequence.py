#!/usr/bin/env python3
"""
Extract DNA sequence from Metrosideros .gbff file
FIXED VERSION - saves in current directory
"""

def extract_sequence_from_gbff(gbff_file, output_fasta):
    """Extract DNA sequence and save as FASTA"""
    
    print(f"Reading: {gbff_file}")
    
    sequence = []
    in_origin = False
    
    with open(gbff_file, 'r') as f:
        for line in f:
            if line.startswith('LOCUS'):
                print(f"Found: {line.strip()}")
            elif line.startswith('ORIGIN'):
                in_origin = True
                print("✓ Extracting sequence...")
            elif line.startswith('//'):
                in_origin = False
            elif in_origin:
                seq_part = ''.join(line.split()[1:])
                sequence.append(seq_part)
    
    full_sequence = ''.join(sequence).upper()
    
    print(f"Extracted {len(full_sequence):,} bp")
    
    with open(output_fasta, 'w') as f:
        f.write(f">Metrosideros_polymorpha_chr1\n")
        for i in range(0, len(full_sequence), 60):
            f.write(full_sequence[i:i+60] + '\n')
    
    import os
    print(f"Saved to: {os.path.abspath(output_fasta)}")
    print(f"   File size: {os.path.getsize(output_fasta):,} bytes")

def main():
    # Use relative paths from where script runs
    gbff_file = r"Genome_Data\Metrosideros\metrosideros.gbff"
    output_fasta = "metrosideros_chr1.fasta"  # Save in current directory!
    
    extract_sequence_from_gbff(gbff_file, output_fasta)

if __name__ == "__main__":
    main()