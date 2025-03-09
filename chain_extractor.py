import os
import json
from Bio.PDB import *
import argparse
import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict, Optional, Any
import sys

class ChainExtractor:
    def __init__(self, data_dir, output_dir):
        """
        Initialize the ChainExtractor.
        
        Args:
            data_dir (str): Directory containing the structure folders
            output_dir (str): Directory where extracted chains will be saved
        """
        self.data_dir = Path(data_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.parser = MMCIFParser(QUIET=True)
        self.io = MMCIFIO()
        
    def get_all_structures(self):
        """Get a list of all structure IDs from the data directory."""
        return [d.name for d in self.data_dir.iterdir() if d.is_dir()]
        
    def get_polymer_classes_from_structure(self, rcsb_id):
        """
        Get a mapping of chain IDs to polymer classes for a structure.
        
        Args:
            rcsb_id (str): The RCSB ID of the structure
            
        Returns:
            dict: A dictionary mapping chain IDs to lists of polymer classes
        """
        try:
            json_file = self.data_dir / rcsb_id / f"{rcsb_id}.json"
            if not json_file.exists():
                print(f"JSON file not found for {rcsb_id}")
                return {}
                
            with open(json_file, "r") as f:
                ribosome_structure = json.load(f)
            
            # Create a dictionary to map chain IDs to polymer classes
            chain_to_polymer_classes = {}
            
            # Process proteins
            if "proteins" in ribosome_structure:
                for protein in ribosome_structure["proteins"]:
                    auth_asym_id = protein["auth_asym_id"]
                    if "nomenclature" in protein and len(protein["nomenclature"]) > 0:
                        chain_to_polymer_classes[auth_asym_id] = protein["nomenclature"]
            
            # Process RNAs
            if "rnas" in ribosome_structure:
                for rna in ribosome_structure["rnas"]:
                    auth_asym_id = rna["auth_asym_id"]
                    if "nomenclature" in rna and len(rna["nomenclature"]) > 0:
                        chain_to_polymer_classes[auth_asym_id] = rna["nomenclature"]
            
            # Process other polymers
            if "other_polymers" in ribosome_structure:
                for polymer in ribosome_structure["other_polymers"]:
                    auth_asym_id = polymer["auth_asym_id"]
                    if "nomenclature" in polymer and len(polymer["nomenclature"]) > 0:
                        chain_to_polymer_classes[auth_asym_id] = polymer["nomenclature"]
            
            return chain_to_polymer_classes
            
        except Exception as e:
            print(f"Error getting polymer classes for {rcsb_id}: {e}")
            return {}
    
    def has_target_polymer_class(self, polymer_classes, target_polymer_class):
        """
        Check if a list of polymer classes contains the target class.
        
        Args:
            polymer_classes (list): List of polymer classes
            target_polymer_class (str): The target polymer class to search for
            
        Returns:
            bool: True if the target class is found, False otherwise
        """
        # Check if the target class is in the list of polymer classes
        for class_entry in polymer_classes:
            # Handle either string or dict representation
            if isinstance(class_entry, str) and class_entry == target_polymer_class:
                return True
            elif isinstance(class_entry, dict) and class_entry.get("name") == target_polymer_class:
                return True
        return False
            
    def extract_chain(self, rcsb_id, chain_id, polymer_class):
        """
        Extract a specific chain from a structure and save it to the output directory.
        
        Args:
            rcsb_id (str): The RCSB ID of the structure
            chain_id (str): The chain ID to extract
            polymer_class (str): The polymer class of the chain (for directory organization)
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            structure_dir = self.data_dir / rcsb_id
            cif_file = structure_dir / f"{rcsb_id}.cif"
            
            if not cif_file.exists():
                print(f"CIF file not found for {rcsb_id}")
                return False
                
            # Parse the structure
            structure = self.parser.get_structure(rcsb_id, cif_file)
            
            # Create a new structure with just the requested chain
            new_structure = Structure.Structure(rcsb_id)
            model = Model.Model(0)
            new_structure.add(model)
            
            # Get the first model
            original_model = structure[0]
            
            # Check if the chain exists
            if chain_id not in original_model:
                print(f"Chain {chain_id} not found in structure {rcsb_id}")
                return False
                
            # Add the chain to the new model
            chain = original_model[chain_id]
            model.add(chain)
            
            # Create output directory for polymer class
            polymer_dir = self.output_dir / polymer_class
            polymer_dir.mkdir(exist_ok=True)
            
            # Save the extracted chain
            output_file = polymer_dir / f"{rcsb_id}_{chain_id}.cif"
            self.io.set_structure(new_structure)
            self.io.save(str(output_file))
            
            print(f"Extracted chain {chain_id} from {rcsb_id} to {output_file}")
            return True
            
        except Exception as e:
            print(f"Error extracting chain {chain_id} from {rcsb_id}: {e}")
            return False
            
    def extract_chains_by_polymer_class(self, target_polymer_class, max_workers=4):
        """
        Extract all chains of a specific polymer class from all structures.
        
        Args:
            target_polymer_class (str): The polymer class to extract
            max_workers (int): Maximum number of parallel workers
            
        Returns:
            int: Number of chains extracted
        """
        extracted_count = 0
        structures = self.get_all_structures()
        print(f"Processing {len(structures)} structures for polymer class: {target_polymer_class}")
        
        # Process structures sequentially to avoid memory issues
        for rcsb_id in structures:
            try:
                # Get the mapping of chains to polymer classes
                chain_to_polymer_classes = self.get_polymer_classes_from_structure(rcsb_id)
                
                # If no polymer classes found, skip this structure
                if not chain_to_polymer_classes:
                    continue
                
                # Find chains that match the target polymer class
                matching_chains = []
                for chain_id, polymer_classes in chain_to_polymer_classes.items():
                    if self.has_target_polymer_class(polymer_classes, target_polymer_class):
                        matching_chains.append(chain_id)
                
                if not matching_chains:
                    continue
                    
                print(f"Found {len(matching_chains)} matching chains in {rcsb_id}")
                
                # Process matching chains in parallel
                tasks = []
                with ProcessPoolExecutor(max_workers=max_workers) as executor:
                    for chain_id in matching_chains:
                        task = executor.submit(
                            self.extract_chain, 
                            rcsb_id, 
                            chain_id,
                            target_polymer_class
                        )
                        tasks.append(task)
                        
                    # Process results
                    for future in as_completed(tasks):
                        if future.result():
                            extracted_count += 1
                            
            except Exception as e:
                print(f"Error processing structure {rcsb_id}: {e}")
                continue
                    
        print(f"Extracted {extracted_count} chains of polymer class {target_polymer_class}")
        return extracted_count
        
    def list_polymer_classes(self, sample_size=10):
        """
        List all unique polymer classes found in the first sample_size structures.
        
        Args:
            sample_size (int): Number of structures to sample
            
        Returns:
            set: Set of unique polymer classes
        """
        structures = self.get_all_structures()[:sample_size]
        print(f"Sampling {len(structures)} structures to find polymer classes")
        
        unique_classes = set()
        
        for rcsb_id in structures:
            try:
                chain_to_polymer_classes = self.get_polymer_classes_from_structure(rcsb_id)
                
                for polymer_classes in chain_to_polymer_classes.values():
                    for class_entry in polymer_classes:
                        if isinstance(class_entry, str):
                            unique_classes.add(class_entry)
                        elif isinstance(class_entry, dict) and "name" in class_entry:
                            unique_classes.add(class_entry["name"])
            except Exception as e:
                print(f"Error processing structure {rcsb_id}: {e}")
                continue
                
        return unique_classes

def main():
    parser = argparse.ArgumentParser(description="Extract chains by polymer class from structure files")
    parser.add_argument("--data_dir", required=True, help="Directory containing structure folders")
    parser.add_argument("--output_dir", required=True, help="Directory where extracted chains will be saved")
    parser.add_argument("--polymer_class", help="Polymer class to extract")
    parser.add_argument("--list_classes", action="store_true", help="List all polymer classes found in the dataset")
    parser.add_argument("--sample_size", type=int, default=10, help="Number of structures to sample when listing classes")
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel workers")
    
    args = parser.parse_args()
    
    extractor = ChainExtractor(args.data_dir, args.output_dir)
    
    if args.list_classes:
        unique_classes = extractor.list_polymer_classes(args.sample_size)
        print("Found the following polymer classes:")
        for cls in sorted(unique_classes):
            print(f"  - {cls}")
    elif args.polymer_class:
        extractor.extract_chains_by_polymer_class(args.polymer_class, args.workers)
    else:
        parser.print_help()
        print("\nError: Either --polymer_class or --list_classes must be specified")
        sys.exit(1)

if __name__ == "__main__":
    main()