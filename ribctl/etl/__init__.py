from enum import Enum
import os
from ribctl import RIBETL_DATA


class AssetType(Enum):
    PTC          = "PTC"
    PROFILE      = "PROFILE"
    CIF          = "CIF"
    LIGAND       = "LIGAND"
    MODIFIED_CIF = "MODIFIED_CIF"
    CHAINS       = "CHAINS"



class AssetFile:
    
    
    def __init__(self, rcsb_id):
        self.rcsb_id = rcsb_id.upper()

    def profile(self):
        return "{}/{}/{}.json".format(RIBETL_DATA,self.rcsb_id, self.rcsb_id)

    def ptc_coords(self):
        return "{}/{}/{}_PTC_COORDINATES.json".format(RIBETL_DATA,self.rcsb_id, self.rcsb_id)

    def cif(self):
        return "{}/{}/{}.cif".format(RIBETL_DATA,self.rcsb_id, self.rcsb_id)
    
    def chain(self, auth_asym_id):
        return "{}/{}/CHAINS/{}.cif".format(RIBETL_DATA,self.rcsb_id, auth_asym_id)

    def ligand(self,ligand_id):
        return "{}/{}/LIGAND_{}.json".format(RIBETL_DATA,self.rcsb_id, ligand_id)
    
    def status(self)->dict[str, bool]:
        return {
           AssetType.PTC.name    : os.path.exists(self.ptc_coords()),
           AssetType.PROFILE.name: os.path.exists(self.profile()),
           AssetType.CIF.name    : os.path.exists(self.cif()),
        }


    @staticmethod
    def list_all_structs():
        return os.listdir(RIBETL_DATA)
    
    @staticmethod
    def status_all()->list[tuple[str, dict[str, bool]]]:
        _ = []
        for struct in AssetFile.list_all_structs():
            _.append([struct, AssetFile(struct).status()])
        return _

