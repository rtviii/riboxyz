import json
from pprint import pprint
import sys
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.ribosome_types.types_ribosome import MitochondrialRNAClass, Polymer



def struct_stats(ra:RibosomeAssets):
    profile     = ra.profile()
    struct_stat = {}

    has_factors = False
    has_trna    = False
    n_dbank_compounds  = 0


    rnas  = profile.rnas
    prots = profile.proteins
    ligs  = profile.nonpolymeric_ligands


    for lig in ligs:
        if lig.nonpolymer_comp.drugbank != None:
            n_dbank_compounds += 1
            print(lig.nonpolymer_comp.drugbank)
        else:
            print("No drugbank")
    
    
    def lsu_ssu_presence(rnas:list[Polymer]):
        has_lsu     = 0
        has_ssu     = 0
        for rna in rnas:
            if profile.mitochondrial :
                if  "mt12SrRNA" in rna.nomenclature:
                    has_ssu = 1
                elif  "mt16SrRNA" in rna.nomenclature:
                    has_lsu = 2

        match has_ssu + has_lsu:
            case 1:
                return "ssu"
            case 2:
                return "lsu"
            case 3:
                return "both"

        
            
            


    



def get_stats():

    global_stats = {
        "mitochondrial"     : 0,
        "cytosolic"         : 0,
        "lsu_only"          : 0,
        "ssu_only"          : 0,
        "ssu_lsu"           : 0,
        "with_trna"         : 0,
        "with_factor"       : 0,
        "drugbank_compounds": 0
    }
    chain_classes={}
       
    for struct in RibosomeAssets.list_all_structs():





struct_stats(RibosomeAssets(sys.argv[1].upper()))