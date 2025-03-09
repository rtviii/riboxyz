import os
import sys
from chimerax.core.commands import run

# ================================================================
# CONFIGURATION - MODIFY THIS SECTION
# ================================================================

# Define the directory where your CIF files are located
homolog_dir = "/Users/rtviii/dev/riboxyz/mL45_polymers/mL45/"

# Set your reference/pivot structure here (without .cif extension)
# Examples: "6ZSB" or "7OI8" or "3J7Y"
# Leave empty ("") to use the first structure in the directory as reference
REFERENCE_STRUCTURE = "6ZSB"  # <-- SET YOUR REFERENCE STRUCTURE HERE

# ================================================================
# END OF CONFIGURATION
# ================================================================

# Get list of all CIF files in the directory
cif_files = [f for f in os.listdir(homolog_dir) if f.endswith(".cif")]

if not cif_files:
    print("No CIF files found in directory:", homolog_dir)
    exit()

# Process arguments if any were passed to the script
if len(sys.argv) > 1:
    arg_ref = sys.argv[1]
    if arg_ref.endswith('.cif'):
        REFERENCE_STRUCTURE = arg_ref[:-4]  # Remove .cif extension
    else:
        REFERENCE_STRUCTURE = arg_ref

# Find the reference structure
ref_idx = 0  # Default to first structure
ref_file = ""

if REFERENCE_STRUCTURE:
    # Try to find the specified reference structure
    matching_files = [f for f in cif_files if REFERENCE_STRUCTURE in f]
    if matching_files:
        ref_file = matching_files[0]
        print(f"Using specified reference structure: {ref_file}")
        
        # Move the reference file to the beginning of the list
        cif_files.remove(ref_file)
        cif_files.insert(0, ref_file)
    else:
        print(f"Warning: Specified reference structure '{REFERENCE_STRUCTURE}' not found.")
        print(f"Using first structure as reference instead: {cif_files[0]}")
        ref_file = cif_files[0]
else:
    # Use the first structure as reference
    ref_file = cif_files[0]
    print(f"Using first structure as reference: {ref_file}")

# Command to open all structures
open_cmd = "open " + " ".join([os.path.join(homolog_dir, cif) for cif in cif_files])
run(session, open_cmd)

# Set all structures to the same representation
run(session, "cartoon")
run(session, "color white")  # Start with all white for simplicity

# Make the reference structure (first one) green for easy identification
run(session, "color #1 green")

print(f"Reference structure: {ref_file}")
print(f"Total structures loaded: {len(cif_files)}")

# Try better matching approach using the correct matchmaker syntax
# Using "to" instead of direct model specifications per documentation
print("Starting superposition using documented matchmaker syntax...")

for i in range(2, len(cif_files) + 1):
    try:
        # Using proper syntax: matchmaker matchstruct to refstruct [options]
        # Adding multiple options for more robust alignment
        match_cmd = f"matchmaker #{i} to #1 alg nw matrix blosum-62 cutoffDistance 5.0 verbose false"
        run(session, match_cmd)
        print(f"Successfully aligned structure #{i}")
    except Exception as e:
        print(f"Warning: Could not align structure #{i}: {str(e)}")
        # Try with alternative parameters
        try:
            alt_cmd = f"matchmaker #{i} to #1 alg sw matrix blosum-45 cutoffDistance 10.0"
            run(session, alt_cmd)
            print(f"Successfully aligned structure #{i} with alternative parameters")
        except Exception as e2:
            print(f"Failed to align structure #{i} with both methods")

# Center the view on all models
run(session, "view")

print("Superposition complete!")