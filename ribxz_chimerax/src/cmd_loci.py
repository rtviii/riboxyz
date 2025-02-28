import urllib.request
import tempfile
import os
import json
from chimerax.core.commands import run

# Base URL for the backend API
BASE_URL = "http://localhost:8000"

def loci(session, rcsb_id, loci_name):
    """Render specific loci (tunnel, PTC, constriction) for a given RCSB ID in ChimeraX.

    Arguments:
        session: The ChimeraX session object.
        rcsb_id (str): The RCSB ID of the structure (e.g., "4UG0").
        loci_name (str): The name of the loci ("ptc", "constriction", or "tunnel").
    """
    # Normalize inputs
    rcsb_id = rcsb_id.upper()
    loci_name = loci_name.lower()

    # Validate loci_name
    if loci_name not in ['ptc', 'constriction', 'tunnel']:
        session.logger.error(f"Invalid loci name: {loci_name}")
        return

    if loci_name == 'tunnel':
        # Handle tunnel visualization as a mesh
        url = f"{BASE_URL}/loci/tunnel_geometry?rcsb_id={rcsb_id}&is_ascii=false&format=vtk"
        try:
            with urllib.request.urlopen(url) as response:
                if response.status == 200:
                    # Save VTK file to a temporary location
                    with tempfile.NamedTemporaryFile(delete=False, suffix='.vtk') as tmp:
                        tmp.write(response.read())
                        tmp_path = tmp.name
                    try:
                        # Load the VTK file as a surface model
                        models = run(session, f"open {tmp_path}")
                        if models:
                            tunnel_model = models[0]
                            # Style the tunnel: gray and semi-transparent
                            run(session, f"color #{tunnel_model.id_string} gray")
                            run(session, f"transparency #{tunnel_model.id_string} 50")
                            session.logger.info(f"Added tunnel mesh for {rcsb_id}")
                    finally:
                        os.remove(tmp_path)  # Clean up temporary file
                else:
                    session.logger.error(f"Failed to fetch tunnel geometry: HTTP {response.status}")
        except urllib.error.URLError as e:
            session.logger.error(f"Network error for tunnel {rcsb_id}: {e}")
    else:
        # Handle PTC or constriction as markers with labels
        endpoint = 'ptc' if loci_name == 'ptc' else 'constriction_site'
        url = f"{BASE_URL}/loci/{endpoint}?rcsb_id={rcsb_id}"
        try:
            with urllib.request.urlopen(url) as response:
                if response.status == 200:
                    data = json.loads(response.read().decode())
                    if 'location' in data:
                        x, y, z = data['location']
                        
                        # Get all existing model IDs to find next available ID
                        models = session.models
                        existing_ids = {m.id[0] for m in models.list() if m.id}
                        
                        # Find the next available ID
                        next_id = 1
                        while next_id in existing_ids:
                            next_id += 1
                        
                        # Create a marker with the next available ID - removing the 'name' parameter
                        marker_cmd = f"marker #{next_id} position {x},{y},{z} radius 3 color red"
                        marker_set = run(session, marker_cmd)
                        
                        if marker_set:
                            # Add a label to the marker using separate label command
                            run(session, f"label #{next_id} text {loci_name.upper()} height 4 color blue")
                            session.logger.info(f"Added {loci_name} landmark for {rcsb_id}")
                        else:
                            session.logger.error(f"Failed to create marker for {loci_name}")
                    else:
                        session.logger.error(f"Invalid data for {loci_name}: 'location' key missing")
                else:
                    session.logger.error(f"Failed to fetch {loci_name} data: HTTP {response.status}")
        except urllib.error.URLError as e:
            session.logger.error(f"Network error for {loci_name} {rcsb_id}: {e}")
        except ValueError as e:
            session.logger.error(f"Invalid coordinates for {loci_name}: {e}")

def register_loci_command(logger):
    """Register the 'loci' command with ChimeraX."""
    from chimerax.core.commands import CmdDesc, register, StringArg
    desc = CmdDesc(
        required=[("rcsb_id", StringArg), ("loci_name", StringArg)],
        synopsis="Render specific loci (tunnel, PTC, constriction) for a given RCSB ID"
    )
    register("loci", desc, loci, logger=logger)

# To use this script in ChimeraX, save it (e.g., as 'loci.py') and run:
#   runscript loci.py
# Then register the command with:
register_loci_command(session.logger)
# Example usage in ChimeraX:
#   loci 4UG0 ptc
#   loci 4UG0 constriction
#   loci 4UG0 tunnel   loci 4UG0 tunnel