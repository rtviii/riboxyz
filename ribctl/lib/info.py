import json
import os
from pprint import pprint
from typing import Literal
import typing

from pydantic import BaseModel
from ribctl import ASSETS, ASSETS_PATH
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.libtax import PhylogenyNode
from ribctl.lib.schema.types_ribosome import (
    RNA,
    CytosolicRNAClass,
    ElongationFactorClass,
    InitiationFactorClass,
    LifecycleFactorClass,
    MitochondrialRNAClass,
    Polymer,
    tRNA,
)
from ribctl.lib.libmsa import Taxid, ncbi


def struct_stats(ra: RibosomeAssets):
    profile = ra.profile()
    struct_stat = {}
    lig_compounds = {}
    n_dbank_compounds = 0

    rnas = profile.rnas
    ligs = profile.nonpolymeric_ligands

    for lig in ligs:
        # If this particualr ligand is a drugbank compound
        chem_name = lig.chemicalName

        if lig.nonpolymer_comp.drugbank != None:
            n_dbank_compounds += 1
            if chem_name not in lig_compounds:
                lig_compounds[chem_name] = 1
            else:
                lig_compounds[chem_name] = lig_compounds[chem_name] + 1
        else:
            ...
            print("No drugbank")

    def lsu_ssu_presence(rnas: list[RNA]) -> Literal["both", "ssu", "lsu"]:
        has_lsu = 0
        has_ssu = 0

        for rna in rnas:
            if profile.mitochondrial:
                if MitochondrialRNAClass.mtrRNA12S in rna.nomenclature:
                    has_ssu = 1
                elif MitochondrialRNAClass.mtrRNA16S in rna.nomenclature:
                    has_lsu = 2
            else:
                if (
                    (CytosolicRNAClass.rRNA_5_8S in rna.nomenclature)
                    or (CytosolicRNAClass.rRNA_5S in rna.nomenclature)
                    or (CytosolicRNAClass.rRNA_28S in rna.nomenclature)
                    or (CytosolicRNAClass.rRNA_25S in rna.nomenclature)
                    or (CytosolicRNAClass.rRNA_23S in rna.nomenclature)
                ):
                    has_lsu = 2
                elif (CytosolicRNAClass.rRNA_16S in rna.nomenclature) or (
                    CytosolicRNAClass.rRNA_18S in rna.nomenclature
                ):
                    has_ssu = 1

        match has_ssu + has_lsu:
            case 1:
                return "ssu"
            case 2:
                return "lsu"
            case 3:
                return "both"
            case _:
                raise Exception("No ssu or lsu found")

    struct_stat["subunit_composition"] = lsu_ssu_presence(rnas)
    struct_stat["mitochondrial"] = profile.mitochondrial

    # for chain in [*rnas, *prots,*other_polys ]:
    #     if chain.assembly_id != 0:
    #         continue
    #     if len(chain.nomenclature) > 0:
    #         cls = chain.nomenclature[0].value

    #         #! class acc
    #         if cls not in nomenclature_classes:
    #             nomenclature_classes[cls] = 1
    #         else:
    #             nomenclature_classes[cls] = nomenclature_classes[cls] + 1

    #         #! subunit determinatiaon
    #         if cls in [k.value for k in list(LifecycleFactorClass)]:
    #             struct_stat["has_factors"] = True

    #         if cls in [k.value for k in list(ElongationFactorClass)]:
    #             struct_stat["has_elongation_factors"] = True

    #         if cls in [k.value for k in list(InitiationFactorClass)]:
    #             struct_stat["has_initiation_factors"] = True

    #         if cls in [k.value for k in list(tRNA)]:
    # #             struct_stat["has_trna"] = True

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

def get_stats()->StructureCompositionStats:
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

    for struct in RibosomeAssets.list_all_structs():
        try:
            [struct_stat, lig_compounds, n_dbank_compounds, superkingdom] = (
                struct_stats(RibosomeAssets(struct))
            )

            # for k, v in nomenclature_classes.items():
            #     if k not in chain_classes:
            #         chain_classes[k] = v
            #     else:
            #         chain_classes[k] += v
            if superkingdom not in list(global_stats.keys()):
                continue
            # ---------- SUBUNUTIS
            if struct_stat["subunit_composition"] == "lsu":
                global_stats[superkingdom]["lsu_only"] += 1

            elif struct_stat["subunit_composition"] == "ssu":
                global_stats[superkingdom]["ssu_only"] += 1

            elif struct_stat["subunit_composition"] == "both":
                global_stats[superkingdom]["ssu_lsu"] += 1

            # if struct_stat["has_trna"] == True:
            #     global_stats["with_trna"] += 1

            # if struct_stat["has_elongation_factors"] == True:
            #     global_stats["with_elongation_factor"] += 1
            # if struct_stat["has_initiation_factors"] == True:
            #     global_stats["with_initiation_factor"] += 1

            # if struct_stat["has_factors"] == True:
            #     global_stats["with_factor"] += 1

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
    })



def run_composition_stats():
    with open(os.path.join(ASSETS_PATH,"structure_composition_stats.json"), "w") as of:
        json.dump(get_stats().model_dump(), of, indent=4)
        print("Saved structure composition_stats: ", os.path.join(ASSETS_PATH,"structure_composition_stats.json"))

if __name__ == "__main__":
    run_composition_stats()