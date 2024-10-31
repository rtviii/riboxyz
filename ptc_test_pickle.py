from pprint import pprint
from ribctl.lib.landmarks.ptc_via_trna import PTC_reference_residues
import pickle
from Bio.PDB.Residue import Residue
import copy

def pickle_residue_array(residues:list[Residue], output_file:str):
    """
    Pickle an array of BioPython Residue instances.
    
    Args:
        residues (list): List of Bio.PDB.Residue objects
        output_file (str): Path to save the pickled data
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # # Create a deep copy to avoid modifying original objects
        # residue_copies = []
        
        # for res in residues:
        #     if not isinstance(res, Residue):
        #         raise TypeError(f"Expected Residue object, got {type(res)}")
                
        #     # Deep copy the residue
        #     res_copy = copy.deepcopy(res)
        #     residue_copies.append(res_copy)
            
        # Pickle the copied residues
        with open(output_file, 'wb') as f:
            pickle.dump(residues, f, protocol=pickle.HIGHEST_PROTOCOL)
        return True
        
    except Exception as e:
        print(f"Error pickling residues: {str(e)}")
        return False

def unpickle_residue_array(input_file:str):
    """
    Unpickle an array of BioPython Residue instances.
    
    Args:
        input_file (str): Path to the pickled data file
        
    Returns:
        list: List of unpickled Residue objects, or None if error occurs
    """
    try:
        with open(input_file, 'rb') as f:
            residues = pickle.load(f)
            
        # Verify we got back Residue objects
        for res in residues:
            if not isinstance(res, Residue):
                raise TypeError(f"Unpickled object is not a Residue: {type(res)}")
                
        return residues
        
    except Exception as e:
        print(f"Error unpickling residues: {str(e)}")
        return None


# c = PTC_reference_residues('euk')
cached_name  = 'ptc_reference_residues_EUK.pickle'
# pickle_residue_array(c, )
c=  unpickle_residue_array(cached_name)

