import json
import os
from pprint import pprint
from typing import Literal
import typing

from pydantic import BaseModel
from ribctl import ASSETS, ASSETS_PATH
from ribctl.etl.etl_assets_ops import Assets, RibosomeOps, Structure
from ribctl.lib.libtax import PhylogenyNode
from ribctl.lib.schema.types_ribosome import (
    RNA,
    CytosolicRNAClass,
    ElongationFactorClass,
    InitiationFactorClass,
    LifecycleFactorClass,
    MitochondrialRNAClass,
    Polymer,
    RibosomeStructureMetadata,
    tRNA,
)
from ribctl.lib.libmsa import Taxid, ncbi

def lsu_ssu_presence(rnas:list[RNA], is_mitochondrial :bool) -> list[Literal[ "ssu", "lsu"]]:
    has_lsu = 0
    has_ssu = 0
    for rna in rnas:
        # print("Inspecting rnas:", rna)
        if is_mitochondrial:
            if MitochondrialRNAClass.mtrRNA12S in rna.nomenclature:
                has_ssu = 1
            elif MitochondrialRNAClass.mtrRNA16S in rna.nomenclature:
                has_lsu = 2
        else:
            if ( (CytosolicRNAClass.rRNA_5_8S in rna.nomenclature)
                or (CytosolicRNAClass.rRNA_5S in rna.nomenclature)
                or (CytosolicRNAClass.rRNA_28S in rna.nomenclature)
                or (CytosolicRNAClass.rRNA_25S in rna.nomenclature)
                or (CytosolicRNAClass.rRNA_23S in rna.nomenclature) ):
                has_lsu = 2
            elif (CytosolicRNAClass.rRNA_16S in rna.nomenclature) or ( CytosolicRNAClass.rRNA_18S in rna.nomenclature ):
                has_ssu = 1
    match has_ssu + has_lsu:
        case 1:
            return [ "ssu" ]
        case 2:
            return [ "lsu" ]
        case 3:
            return ['ssu','lsu']
        case 0:
            return []
        case _:
            raise ValueError("Invalid case")

def struct_stats(ra: RibosomeOps):
    profile = ra.profile()
    struct_stat = {}
    lig_compounds = {}
    n_dbank_compounds = 0

    ligs = profile.nonpolymeric_ligands

    for lig in ligs:
        # If this particualr ligand is a drugbank compound
        chem_name = lig.chemicalName

        if lig.nonpolymer_comp != None and lig.nonpolymer_comp.drugbank != None:
            n_dbank_compounds += 1
            if chem_name not in lig_compounds:
                lig_compounds[chem_name] = 1
            else:
                lig_compounds[chem_name] = lig_compounds[chem_name] + 1
        else:
            ...
            print("No drugbank")


    struct_stat["mitochondrial"]       = profile.mitochondrial
    struct_stat["subunit_composition"] = lsu_ssu_presence(profile.rnas, profile.mitochondrial)


    return [
        struct_stat,
        lig_compounds,
        n_dbank_compounds,
        Taxid.superkingdom(profile.src_organism_ids[0]),
    ]

class CompositionStats(BaseModel):
    lsu_only          :int
    ssu_only          :int
    ssu_lsu           :int
    drugbank_compounds:int
    mitochondrial     :int

class StructureCompositionStats(BaseModel):

    archaea  : CompositionStats
    bacteria : CompositionStats
    eukaryota: CompositionStats

def get_stats():
    global_stats = {
        "bacteria": {
            "lsu_only": 0,
            "ssu_only": 0,
            "ssu_lsu": 0,
            "drugbank_compounds": 0,
            "mitochondrial": 0,
        },
        "eukaryota": {
            "lsu_only": 0,
            "ssu_only": 0,
            "ssu_lsu": 0,
            "drugbank_compounds": 0,
            "mitochondrial": 0,
        },
        "archaea": {
            "lsu_only": 0,
            "ssu_only": 0,
            "ssu_lsu": 0,
            "drugbank_compounds": 0,
            "mitochondrial": 0,
        },
    }
    lig_global = {
        "archaea": {},
        "eukaryota": {},
        "bacteria": {},
    }

    for struct in Assets.list_all_structs():
        try:
            [struct_stat, lig_compounds, n_dbank_compounds, superkingdom] = (struct_stats(RibosomeOps(struct)))

            if superkingdom not in list(global_stats.keys()):
                continue
            # ---------- SUBUNUTIS
            if struct_stat["subunit_composition"] == "lsu":
                global_stats[superkingdom]["lsu_only"] += 1

            elif struct_stat["subunit_composition"] == "ssu":
                global_stats[superkingdom]["ssu_only"] += 1

            elif struct_stat["subunit_composition"] == "both":
                global_stats[superkingdom]["ssu_lsu"] += 1

            if struct_stat["mitochondrial"] == True:
                global_stats[superkingdom]["mitochondrial"] += 1

            global_stats[superkingdom]["drugbank_compounds"] += n_dbank_compounds

            for k, v in lig_compounds.items():
                if k not in lig_global:
                    lig_global[superkingdom][k] = v
                else:
                    lig_global[superkingdom][k] += v

        except Exception as e:
            print("Error --->", e)

    return StructureCompositionStats.model_validate({
        **global_stats,
        "ligands": lig_global,
        # "chain_classes": {key: chain_classes[key] for key in sorted(chain_classes)},
    }).model_dump()

def run_composition_stats():
    with open(os.path.join(ASSETS_PATH,"structure_composition_stats.json"), "w") as of:
        d = get_stats()
        json.dump(d, of, indent=4)
        print("Saved structure composition_stats: ", os.path.join(ASSETS_PATH,"structure_composition_stats.json"))

if __name__ == "__main__":
    run_composition_stats()