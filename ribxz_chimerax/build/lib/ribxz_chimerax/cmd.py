# ribxz_chimerax/cmd.py

from chimerax.core.commands import CmdDesc, StringArg
from .io import fetch_ribosome_data
import asyncio

def ribxz_repr(session, pdb_id):
    """
    Command to fetch ribosome metadata and apply representation based on it.
    
    Parameters
    ----------
    session : chimerax.core.session.Session
        The current ChimeraX session
    pdb_id : str
        The PDB ID to fetch metadata for (e.g. "4UG0")
    """
    # Log start of operation
    session.logger.info(f"Fetching ribosome metadata for {pdb_id}...")
    
    # Get the event loop to run our async function
    loop = asyncio.get_event_loop()
    metadata = loop.run_until_complete(fetch_ribosome_data(pdb_id))
    
    if "error" in metadata:
        session.logger.error(metadata["error"])
        return
    
    # Example of using the metadata to apply representation
    # This will depend on the structure of your metadata
    try:
        # First, open the structure if not already open
        from chimerax.atomic import AtomicStructure
        structures = session.models.list(type=AtomicStructure)
        structure_loaded = False
        
        for s in structures:
            if s.pdb_id.lower() == pdb_id.lower():
                structure_loaded = True
                structure = s
                break
        
        if not structure_loaded:
            # Open the structure from PDB
            from chimerax.core.commands import run
            run(session, f"open {pdb_id}")
            # Get the newly opened structure
            structures = session.models.list(type=AtomicStructure)
            for s in structures:
                if s.pdb_id.lower() == pdb_id.lower():
                    structure = s
                    break
        
        # Apply visualization based on metadata
        # This is just an example - customize based on your metadata structure
        if "rRNA_chains" in metadata:
            # Color rRNA chains
            chain_ids = ",".join(metadata["rRNA_chains"])
            run(session, f"color bychain #{structure.id}/{chain_ids}")
        
        if "functional_sites" in metadata:
            for site in metadata["functional_sites"]:
                # Highlight functional sites
                residues = site["residues"]
                res_spec = "+".join(residues)
                run(session, f"color red #{structure.id}:{res_spec}")
                run(session, f"style ball #{structure.id}:{res_spec}")
        
        # Display any additional information in the log
        session.logger.info(f"Applied ribosome representation for {pdb_id}")
        if "description" in metadata:
            session.logger.info(f"Description: {metadata['description']}")
            
    except Exception