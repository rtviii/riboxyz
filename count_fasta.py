import os
from pathlib import Path
import time
from typing import List, Tuple

from ribctl.lib.libmsa import Fasta, muscle_align_N_seq

def align_fasta_files(directory: str, num_files: int = 10) -> None:
    """
    Create multiple sequence alignments for FASTA files using MUSCLE (dry run)
    
    Args:
        directory (str): Path to directory containing FASTA files
        num_files (int): Number of files to process as a test
    """
    # Get all files with common FASTA extensions
    fasta_files = []
    for ext in ['.fasta', '.fa', '.fna', '.faa']:
        fasta_files.extend(Path(directory).glob(f'*{ext}'))
    
    if not fasta_files:
        print(f"No FASTA files found in {directory}")
        return
    
    # Process only the specified number of files
    for fasta_path in sorted(fasta_files)[:num_files]:
        try:
            print(f"\nProcessing {fasta_path.name}...")
            
            # Time the loading of sequences
            start_time = time.time()
            fasta = Fasta(str(fasta_path))
            load_time = time.time() - start_time
            print(f"Loading {len(fasta.records)} sequences took {load_time:.2f} seconds")
            
            # Time the alignment
            start_time = time.time()
            aligned_records = muscle_align_N_seq(fasta.records[:10])
            print(aligned_records)
            align_time = time.time() - start_time
            
            # Calculate output path (not saving, just for demonstration)
            output_path = fasta_path.with_suffix('.afasta')
            
            print(f"Alignment completed in {align_time:.2f} seconds")
            print(f"Would save aligned sequences to: {output_path}")
            
        except Exception as e:
            print(f"Error processing {fasta_path.name}: {str(e)}")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 2:
        print("Usage: python script.py <fasta_directory>")
        sys.exit(1)
    
    directory = sys.argv[1]
    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory")
        sys.exit(1)
        
    align_fasta_files(directory)