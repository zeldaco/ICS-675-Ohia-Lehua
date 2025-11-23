#!/usr/bin/env python3
"""
ICS 675 Final Project - Chitinase Sequence Alignment
Smith-Waterman Local Alignment Implementation
Author: Zelda Cole
Date: November 2025

This implementation compares two chitinase sequences from Eucalyptus grandis:
1. Basic endochitinase (Class I) - XP_010040215.1
2. Endochitinase EP3 (Class IV) - XP_010036821.1
"""

import numpy as np
import time
import sys
import argparse
from typing import Tuple, List
import matplotlib.pyplot as plt
import seaborn as sns

class SmithWatermanAligner:
    def __init__(self, match_score=2, mismatch_score=-1, gap_penalty=-1):
        """
        Initialize the Smith-Waterman aligner with scoring parameters.
        
        Args:
            match_score: Score for matching residues
            mismatch_score: Score for mismatching residues
            gap_penalty: Penalty for gaps
        """
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty
        
        # For benchmarking
        self.runtime = 0
        self.memory_used = 0
        
    def read_fasta(self, filename: str) -> Tuple[str, str]:
        """
        Read a FASTA file and return header and sequence.
        
        Args:
            filename: Path to FASTA file
            
        Returns:
            Tuple of (header, sequence)
        """
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        header = lines[0].strip()
        sequence = ''.join([line.strip() for line in lines[1:]])
        return header, sequence
    
    def scoring_function(self, a: str, b: str) -> int:
        """
        Score function for amino acid pairs.
        
        Args:
            a, b: Amino acid characters
            
        Returns:
            Score for the pair
        """
        if a == b:
            return self.match_score
        else:
            return self.mismatch_score
    
    def smith_waterman(self, seq1: str, seq2: str) -> Tuple[int, np.ndarray]:
        """
        Perform Smith-Waterman local alignment.
        
        Args:
            seq1, seq2: Input sequences
            
        Returns:
            Tuple of (max_score, scoring_matrix)
        """
        start_time = time.time()
        
        m, n = len(seq1), len(seq2)
        
        # Initialize scoring matrix with zeros
        V = np.zeros((m + 1, n + 1), dtype=int)
        
        # Track memory usage (approximate)
        self.memory_used = V.nbytes
        
        # Fill the matrix using Smith-Waterman recurrence
        max_score = 0
        max_pos = (0, 0)
        
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                # Calculate scores for all possible moves
                diagonal = V[i-1, j-1] + self.scoring_function(seq1[i-1], seq2[j-1])
                vertical = V[i-1, j] + self.gap_penalty
                horizontal = V[i, j-1] + self.gap_penalty
                
                # Smith-Waterman: take max of all scores and 0
                V[i, j] = max(diagonal, vertical, horizontal, 0)
                
                # Track maximum score and position
                if V[i, j] > max_score:
                    max_score = V[i, j]
                    max_pos = (i, j)
        
        self.runtime = time.time() - start_time
        return max_score, V
    
    def traceback(self, V: np.ndarray, seq1: str, seq2: str, max_pos: Tuple[int, int]) -> Tuple[str, str]:
        """
        Trace back through the matrix to find optimal alignment.
        
        Args:
            V: Scoring matrix
            seq1, seq2: Input sequences
            max_pos: Position of maximum score
            
        Returns:
            Tuple of aligned sequences
        """
        i, j = max_pos
        aligned_seq1 = []
        aligned_seq2 = []
        
        # Trace back until we hit a zero
        while i > 0 and j > 0 and V[i, j] > 0:
            current_score = V[i, j]
            diagonal = V[i-1, j-1] + self.scoring_function(seq1[i-1], seq2[j-1])
            vertical = V[i-1, j] + self.gap_penalty
            horizontal = V[i, j-1] + self.gap_penalty
            
            if current_score == diagonal:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append(seq2[j-1])
                i -= 1
                j -= 1
            elif current_score == vertical:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append('-')
                i -= 1
            else:  # horizontal
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j-1])
                j -= 1
        
        return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))
    
    def align_sequences(self, file1: str, file2: str) -> dict:
        """
        Complete alignment workflow.
        
        Args:
            file1, file2: FASTA file paths
            
        Returns:
            Dictionary with alignment results
        """
        # Read sequences
        header1, seq1 = self.read_fasta(file1)
        header2, seq2 = self.read_fasta(file2)
        
        print(f"Aligning sequences:")
        print(f"Seq1: {header1}")
        print(f"Seq2: {header2}")
        print(f"Sequence lengths: {len(seq1)} x {len(seq2)}")
        
        # Perform alignment
        max_score, matrix = self.smith_waterman(seq1, seq2)
        
        # Find max position for traceback
        max_pos = np.unravel_index(np.argmax(matrix), matrix.shape)
        
        # Traceback to get alignment
        aligned_seq1, aligned_seq2 = self.traceback(matrix, seq1, seq2, max_pos)
        
        # Calculate statistics
        identity = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')
        alignment_length = len(aligned_seq1)
        percent_identity = (identity / alignment_length) * 100 if alignment_length > 0 else 0
        
        return {
            'seq1_header': header1,
            'seq2_header': header2,
            'seq1_length': len(seq1),
            'seq2_length': len(seq2),
            'max_score': max_score,
            'matrix': matrix,
            'aligned_seq1': aligned_seq1,
            'aligned_seq2': aligned_seq2,
            'alignment_length': alignment_length,
            'identity': identity,
            'percent_identity': percent_identity,
            'runtime': self.runtime,
            'memory_used': self.memory_used
        }

def benchmark_alignment(file1: str, file2: str) -> dict:
    """Run benchmarking analysis."""
    aligner = SmithWatermanAligner()
    results = aligner.align_sequences(file1, file2)
    
    print("\n" + "="*60)
    print("BENCHMARKING RESULTS")
    print("="*60)
    print(f"Input sequences: {results['seq1_length']} x {results['seq2_length']} amino acids")
    print(f"Matrix size: {(results['seq1_length']+1) * (results['seq2_length']+1):,} cells")
    print(f"Memory usage: {results['memory_used']:,} bytes ({results['memory_used']/1024:.1f} KB)")
    print(f"Runtime: {results['runtime']:.4f} seconds")
    print(f"Max alignment score: {results['max_score']}")
    print(f"Alignment length: {results['alignment_length']}")
    print(f"Percent identity: {results['percent_identity']:.1f}%")
    
    return results

def create_alignment_visualization(results: dict, output_file: str = "alignment_heatmap.png"):
    """Create visualization of alignment matrix."""
    plt.figure(figsize=(12, 8))
    
    # Create heatmap of scoring matrix
    matrix = results['matrix']
    plt.subplot(2, 1, 1)
    sns.heatmap(matrix[:50, :50], cmap='viridis', cbar_kws={'label': 'Alignment Score'})
    plt.title(f"Smith-Waterman Scoring Matrix (50x50 subset)\nMax Score: {results['max_score']}")
    plt.xlabel("Sequence 2 Position")
    plt.ylabel("Sequence 1 Position")
    
    # Show alignment snippet
    plt.subplot(2, 1, 2)
    alignment_snippet = min(60, len(results['aligned_seq1']))
    seq1_snippet = results['aligned_seq1'][:alignment_snippet]
    seq2_snippet = results['aligned_seq2'][:alignment_snippet]
    
    # Create alignment visualization
    match_line = ''.join(['|' if a == b and a != '-' else ' ' for a, b in zip(seq1_snippet, seq2_snippet)])
    
    plt.text(0.02, 0.8, f"Seq1: {seq1_snippet}", fontfamily='monospace', transform=plt.gca().transAxes)
    plt.text(0.02, 0.6, f"      {match_line}", fontfamily='monospace', transform=plt.gca().transAxes)
    plt.text(0.02, 0.4, f"Seq2: {seq2_snippet}", fontfamily='monospace', transform=plt.gca().transAxes)
    plt.text(0.02, 0.2, f"Identity: {results['percent_identity']:.1f}%", transform=plt.gca().transAxes)
    
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.axis('off')
    plt.title("Local Alignment Result (First 60 positions)")
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"\nVisualization saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Smith-Waterman Chitinase Alignment')
    parser.add_argument('seq1', help='First FASTA file')
    parser.add_argument('seq2', help='Second FASTA file')
    parser.add_argument('--output', '-o', default='alignment_results.png', help='Output visualization file')
    
    args = parser.parse_args()
    
    print("ICS 675 Final Project - Chitinase Sequence Alignment")
    print("Smith-Waterman Local Alignment Implementation")
    print("="*60)
    
    # Run benchmark
    results = benchmark_alignment(args.seq1, args.seq2)
    
    # Create visualization
    create_alignment_visualization(results, args.output)
    
    print(f"\nAlignment completed successfully!")
    print(f"Results visualization: {args.output}")

if __name__ == "__main__":
    main()
