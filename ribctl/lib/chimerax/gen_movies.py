import os
import time
from chimerax.core.commands import run

def verify_structure_loaded(session):
    """Check if any models are loaded in the session"""
    return len(session.models.list()) > 0

def process_structure(session, struct_id, base_dir):
    """Process a single structure with verification steps between commands"""
    print(f"Starting {struct_id} at {time.strftime('%H:%M:%S')}")
    
    # Step 1: Run ribetl to open the structure
    run(session, "view")
    run(session, f"ribetl {struct_id}")
    
    # Quick check for structure loading
    for _ in range(5):  # Wait up to 5 seconds
        if verify_structure_loaded(session):
            break
        time.sleep(0.2)
    
    if not verify_structure_loaded(session):
        raise Exception("Structure failed to load")
    
    print(f"{struct_id}: Structure loaded successfully")
    
    # Step 2: Apply symmetry
    run(session, "sym #1 assembly 1")
    time.sleep(0.2)
    
    # Step 3: Run representation
    run(session, "ribrep #2")
    time.sleep(0.2)
    
    # Step 4: Record and rotate
    print(f"{struct_id}: Starting movie recording")
    run(session, "movie record")
    run(session, "turn y 2 180")
    run(session, "wait 180")
    
    # Step 5: Encode movie
    print(f"{struct_id}: Encoding movie")
    run(session, f"movie encode {base_dir}/{struct_id}/{struct_id}.mp4")
    
    # Quick check for movie file
    mp4_path = os.path.join(base_dir, struct_id, f"{struct_id}.mp4")
    for _ in range(5):  # Wait up to 5 seconds
        if os.path.exists(mp4_path) and os.path.getsize(mp4_path) > 0:
            break
        time.sleep(0.2)
    
    # Clean up
    run(session, "close all")
    time.sleep(0.2)
    
    # Final verification
    if not os.path.exists(mp4_path) or os.path.getsize(mp4_path) == 0:
        raise Exception("Movie file not created or empty")
    
    print(f"Successfully completed {struct_id} at {time.strftime('%H:%M:%S')}")

def get_structures_without_movies(base_dir):
    """Returns a list of structure IDs that don't have corresponding .mp4 files"""
    structures = []
    for item in sorted(os.listdir(base_dir)):
        dir_path = os.path.join(base_dir, item)
        if os.path.isdir(dir_path):
            mp4_path = os.path.join(dir_path, f"{item}.mp4")
            if not os.path.exists(mp4_path):
                structures.append(item)
    return structures

def process_missing_movies(session):
    """Main loop to process structures that don't have movies"""
    while True:
        try:
            base_dir = os.environ.get("RIBETL_DATA")
            if not base_dir:
                print("ERROR: RIBETL_DATA environment variable not set")
                return
            
            structures = get_structures_without_movies(base_dir)
            
            if not structures:
                print("All structures processed! Checking again in 10 seconds...")
                time.sleep(10)
                continue
            
            print(f"Found {len(structures)} structures needing processing")
            
            for struct_id in structures:
                try:
                    process_structure(session, struct_id, base_dir)
                    time.sleep(0.5)  # Small delay between structures
                    
                except Exception as e:
                    print(f"Error processing {struct_id}: {str(e)}")
                    run(session, "close all")  # Cleanup on error
                    time.sleep(1)  # Brief delay after errors
                    continue
            
        except Exception as e:
            print(f"Error in main loop: {str(e)}")
            time.sleep(1)
            continue

# Start processing
process_missing_movies(session)
