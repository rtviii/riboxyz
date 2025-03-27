from pathlib import Path
import json
from typing import List, Dict, Set
import argparse
import os

from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.npet.npet_pipeline import create_npet_mesh

def load_structure_ids(filepath: str) -> List[str]:
    """Load structure IDs from a JSON file, handling UTF-8 BOM if present."""
    with open(filepath, 'r', encoding='utf-8-sig') as f:
        return json.load(f)

def find_missing_meshes(structure_ids: List[str]) -> List[str]:
    """Find structures that are missing NPET_MESH artifacts."""
    missing_meshes = []
    
    for rcsb_id in structure_ids:
        mesh_path = AssetType.NPET_MESH.get_path(rcsb_id)
        if not mesh_path.exists():
            missing_meshes.append(rcsb_id)
    
    return missing_meshes

def process_missing_meshes(missing_ids: List[str], force: bool = False, log_dir: Path = None) -> Dict[str, bool]:
    """Process the missing meshes using create_npet_mesh function."""
    results = {}
    
    for rcsb_id in missing_ids:
        print(f"Processing {rcsb_id}...")
        try:
            tracker = create_npet_mesh(rcsb_id, log_dir=log_dir, force=force)
            results[rcsb_id] = tracker.success
        except Exception as e:
            print(f"Error processing {rcsb_id}: {str(e)}")
            results[rcsb_id] = False
    
    return results

def report_results(category_results: Dict[str, Dict[str, List[str]]]) -> None:
    """Generate a report of the processing results."""
    print("\n=== Processing Report ===")
    
    total_missing = 0
    total_processed = 0
    total_successful = 0
    
    for category, results in category_results.items():
        missing = len(results['missing'])
        processed = len(results['processed'])
        successful = len(results['successful'])
        failed = len(results['failed'])
        
        print(f"\n{category} structures:")
        print(f"  - Missing meshes: {missing}")
        print(f"  - Processed: {processed}")
        print(f"  - Successful: {successful}")
        print(f"  - Failed: {failed}")
        
        if failed > 0:
            print(f"  - Failed IDs: {', '.join(results['failed'])}")
        
        total_missing += missing
        total_processed += processed
        total_successful += successful
    
    print("\nOverall summary:")
    print(f"  - Total missing meshes: {total_missing}")
    print(f"  - Total processed: {total_processed}")
    print(f"  - Total successful: {total_successful}")

def main():
    parser = argparse.ArgumentParser(description="Process missing NPET meshes for different categories of structures")
    parser.add_argument("--archaeal", type=str, help="Path to archaeal structure IDs JSON file")
    parser.add_argument("--mitochondrial", type=str, help="Path to mitochondrial structure IDs JSON file")
    parser.add_argument("--eukaryotic", type=str, help="Path to eukaryotic structure IDs JSON file")
    parser.add_argument("--bacterial", type=str, help="Path to bacterial structure IDs JSON file")
    parser.add_argument("--process", action="store_true", help="Process the missing meshes")
    parser.add_argument("--force", action="store_true", help="Force regeneration of meshes")
    parser.add_argument("--log-dir", type=str, default="logs", help="Directory to store logs")
    args = parser.parse_args()
    
    category_files = {
        "Archaeal": args.archaeal,
        "Mitochondrial": args.mitochondrial,
        "Eukaryotic": args.eukaryotic,
        "Bacterial": args.bacterial
    }
    
    # Filter out None values
    category_files = {k: v for k, v in category_files.items() if v is not None}
    
    if not category_files:
        print("Error: At least one category file must be provided.")
        parser.print_help()
        return
    
    log_dir = Path(args.log_dir)
    log_dir.mkdir(exist_ok=True, parents=True)
    
    category_results = {}
    
    for category, file_path in category_files.items():
        print(f"\nProcessing {category} structures...")
        
        structure_ids = load_structure_ids(file_path)
        print(f"Loaded {len(structure_ids)} {category.lower()} structure IDs")
        
        missing_ids = find_missing_meshes(structure_ids)
        print(f"Found {len(missing_ids)} {category.lower()} structures missing NPET meshes")
        
        category_results[category] = {
            "missing": missing_ids,
            "processed": [],
            "successful": [],
            "failed": []
        }
        
        if args.process and missing_ids:
            print(f"Processing {len(missing_ids)} missing {category.lower()} meshes...")
            results = process_missing_meshes(missing_ids, args.force, log_dir)
            
            category_results[category]["processed"] = list(results.keys())
            category_results[category]["successful"] = [id for id, success in results.items() if success]
            category_results[category]["failed"] = [id for id, success in results.items() if not success]
    
    report_results(category_results)
    
    # Export results to JSON for later analysis
    results_path = log_dir / "missing_mesh_results.json"
    with open(results_path, 'w') as f:
        json.dump(category_results, f, indent=2)
    print(f"\nDetailed results saved to {results_path}")

if __name__ == "__main__":
    main()