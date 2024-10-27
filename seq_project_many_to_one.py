from pprint import pprint
import sys
from typing import NewType, TypeVar
from Bio.PDB.Residue import Residue
from ribctl.etl.etl_assets_ops import RibosomeOps
from ribctl.lib.libbsite import map_motifs, bsite_ligand, bsite_transpose
from ribctl.lib.libseq import BiopythonChain
from ribctl.lib.schema.types_binding_site import (
    PredictionTarget,
    ResiduesMapping,
    ResidueSummary,
)
from ribctl.lib.schema.types_ribosome import PolymerClass
from functools import partial

sys.path.append("/home/rtviii/dev/riboxyz")
from ribctl.lib.libmsa import Fasta
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Pool
from functools import partial


# ? Ligand Projections
# Ok let's not complicate things. The prototype is thus:
# take a motif in 7K00 protein, say uL10

# 1. To Record:
# - locus name and its anchor structure/chain/residue
# - project anchor into the MSA of homologs of other structure
# - project anchor into each other structure PAIRWISE

# 2. Backtrack to structural files indices for each structure, record as a row in a table.

s = Fasta.poly_class_all_seq(PolymerClass("uL4"))

type ResidiueIndices = list[int]
type RCSB_ID = str


tetracycline_structs = [
    (
        {
            "InChI": "InChI=1S/C21H25ClN2O7/c1-24(2)15-6-5-7-11(17(27)10(6)19(29)14(20(15)30)21(23)31)18(28)13-9(25)4-3-8(22)12(13)16(7)26/h3-4,6-7,10-11,14-18,25-28H,5H2,1-2H3,(H2,23,31)/t6-,7+,10-,11+,14+,15+,16+,17+,18-/m1/s1",
            "InChIKey": "DMSXRSJJYGANOR-VAVXRSHRSA-N",
            "SMILES": "CN(C)C1C2CC3C(c4c(ccc(c4C(C3C(C2C(=O)C(C1=O)C(=O)N)O)O)O)Cl)O",
            "SMILES_stereo": "CN(C)[C@H]1[C@@H]2C[C@@H]3[C@@H](c4c(ccc(c4C(C3C([C@@H]2C(=O)C(C1=O)C(=O)N)O)O)O)Cl)O",
            "chemicalId": "D2C",
            "chemicalName": "(2S,4S,4AR,5AS,6S,11R,11AS,12R,12AR)-7-CHLORO-4-(DIMETHYLAMINO)-6,10,11,12-TETRAHYDROXY-1,3-DIOXO-1,2,3,4,4A,5,5A,6,11,11A,12,12A-DODECAHYDROTETRACENE-2-CARBOXAMIDE",
            "formula_weight": 0.453,
            "number_of_instances": 1,
            "pdbx_description": "(2S,4S,4AR,5AS,6S,11R,11AS,12R,12AR)-7-CHLORO-4-(DIMETHYLAMINO)-6,10,11,12-TETRAHYDROXY-1,3-DIOXO-1,2,3,4,4A,5,5A,6,11,11A,12,12A-DODECAHYDROTETRACENE-2-CARBOXAMIDE",
        },
        [
            {
                "rcsb_id": "2F4V",
                "tax_node": {
                    "ncbi_tax_id": 274,
                    "rank": "species",
                    "scientific_name": "Thermus thermophilus",
                },
            }
        ],
    ),
    (
        {
            "InChI": "InChI=1S/C22H24N2O8/c1-21(31)8-5-4-6-11(25)12(8)16(26)13-9(21)7-10-15(24(2)3)17(27)14(20(23)30)19(29)22(10,32)18(13)28/h4-6,9-10,15,25,27-28,31-32H,7H2,1-3H3,(H2,23,30)/t9-,10-,15-,21+,22-/m0/s1",
            "InChIKey": "OFVLGDICTFRJMM-WESIUVDSSA-N",
            "SMILES": "CC1(c2cccc(c2C(=O)C3=C(C4(C(CC31)C(C(=C(C4=O)C(=O)N)O)N(C)C)O)O)O)O",
            "SMILES_stereo": "C[C@]1(c2cccc(c2C(=O)C3=C([C@]4([C@@H](C[C@@H]31)C(C(=C(C4=O)C(=O)N)O)N(C)C)O)O)O)O",
            "chemicalId": "TAC",
            "chemicalName": "TETRACYCLINE",
            "drugbank_description": "Tetracycline is a broad spectrum polyketide "
            "antibiotic produced by the Streptomyces genus of "
            "Actinobacteria. It exerts a bacteriostatic effect "
            "on bacteria by binding reversible to the bacterial "
            "30S ribosomal subunit and blocking incoming "
            "aminoacyl tRNA from binding to the ribosome "
            "acceptor site. It also binds to some extent to the "
            "bacterial 50S ribosomal subunit and may alter the "
            "cytoplasmic membrane causing intracellular "
            "components to leak from bacterial cells.  The FDA "
            "withdrew its approval for the use of all liquid "
            "oral drug products formulated for pediatric use "
            "containing tetracycline in a concentration greater "
            "than 25 mg/ml.[L43942] Other formulations of "
            "tetracycline continue to be used.",
            "drugbank_id": "DB00759",
            "formula_weight": 0.444,
            "number_of_instances": 2,
            "pdbx_description": "TETRACYCLINE",
        },
        [
            {
                "rcsb_id": "1HNW",
                "tax_node": {
                    "ncbi_tax_id": 274,
                    "rank": "species",
                    "scientific_name": "Thermus thermophilus",
                },
            },
            {
                "rcsb_id": "5J5B",
                "tax_node": {
                    "ncbi_tax_id": 83333,
                    "rank": "strain",
                    "scientific_name": "Escherichia coli K-12",
                },
            },
            {
                "rcsb_id": "1I97",
                "tax_node": {
                    "ncbi_tax_id": 274,
                    "rank": "species",
                    "scientific_name": "Thermus thermophilus",
                },
            },
            {
                "rcsb_id": "4V9A",
                "tax_node": {
                    "ncbi_tax_id": 300852,
                    "rank": "strain",
                    "scientific_name": "Thermus thermophilus HB8",
                },
            },
            {
                "rcsb_id": "8CGJ",
                "tax_node": {
                    "ncbi_tax_id": 679895,
                    "rank": "no rank",
                    "scientific_name": "Escherichia coli BW25113",
                },
            },
            {
                "rcsb_id": "5J7L",
                "tax_node": {
                    "ncbi_tax_id": 83333,
                    "rank": "strain",
                    "scientific_name": "Escherichia coli K-12",
                },
            },
        ],
    ),
    (
        {
            "InChI": "InChI=1S/C27H31FN4O8/c1-31(2)20-13-8-11-7-12-14(28)9-15(30-16(33)10-32-5-3-4-6-32)21(34)18(12)22(35)17(11)24(37)27(13,40)25(38)19(23(20)36)26(29)39/h9,11,13,20,34-35,38,40H,3-8,10H2,1-2H3,(H2,29,39)(H,30,33)/t11-,13-,20-,27-/m0/s1",
            "InChIKey": "AKLMFDDQCHURPW-ISIOAQNYSA-N",
            "SMILES": "CN(C)C1C2CC3Cc4c(cc(c(c4C(=C3C(=O)C2(C(=C(C1=O)C(=O)N)O)O)O)O)NC(=O)CN5CCCC5)F",
            "SMILES_stereo": "CN(C)[C@H]1[C@@H]2C[C@@H]3Cc4c(cc(c(c4C(=C3C(=O)[C@@]2(C(=C(C1=O)C(=O)N)O)O)O)O)NC(=O)CN5CCCC5)F",
            "chemicalId": "YQM",
            "chemicalName": "Eravacycline",
            "formula_weight": 0.559,
            "number_of_instances": 1,
            "pdbx_description": "Eravacycline",
        },
        [
            {
                "rcsb_id": "7M4U",
                "tax_node": {
                    "ncbi_tax_id": 480119,
                    "rank": "strain",
                    "scientific_name": "Acinetobacter baumannii AB0057",
                },
            },
            {
                "rcsb_id": "7M4Y",
                "tax_node": {
                    "ncbi_tax_id": 480119,
                    "rank": "strain",
                    "scientific_name": "Acinetobacter baumannii AB0057",
                },
            },
            {
                "rcsb_id": "7M4Z",
                "tax_node": {
                    "ncbi_tax_id": 480119,
                    "rank": "strain",
                    "scientific_name": "Acinetobacter baumannii AB0057",
                },
            },
            {
                "rcsb_id": "7M4V",
                "tax_node": {
                    "ncbi_tax_id": 480119,
                    "rank": "strain",
                    "scientific_name": "Acinetobacter baumannii AB0057",
                },
            },
            {
                "rcsb_id": "7M4W",
                "tax_node": {
                    "ncbi_tax_id": 480119,
                    "rank": "strain",
                    "scientific_name": "Acinetobacter baumannii AB0057",
                },
            },
            {
                "rcsb_id": "7M4X",
                "tax_node": {
                    "ncbi_tax_id": 480119,
                    "rank": "strain",
                    "scientific_name": "Acinetobacter baumannii AB0057",
                },
            },
        ],
    ),
    (
        {
            "InChI": "InChI=1S/C28H32F3N3O7/c1-3-34(4-2)21-14-9-11-8-13-18(16(35)10-12(15-6-5-7-33-15)20(13)28(29,30)31)22(36)17(11)24(38)27(14,41)25(39)19(23(21)37)26(32)40/h10-11,14-15,21,33,35,37-38,41H,3-9H2,1-2H3,(H2,32,40)/t11-,14-,15-,21-,27-/m0/s1",
            "InChIKey": "IDWTZCZXXFOLNV-DOYYSQEVSA-N",
            "SMILES": "CCN(CC)C1C2CC3Cc4c(c(cc(c4C(F)(F)F)C5CCCN5)O)C(=O)C3=C(C2(C(=O)C(=C1O)C(=O)N)O)O",
            "SMILES_stereo": "CCN(CC)[C@H]1[C@@H]2C[C@@H]3Cc4c(c(cc(c4C(F)(F)F)[C@@H]5CCCN5)O)C(=O)C3=C([C@@]2(C(=O)C(=C1O)C(=O)N)O)O",
            "chemicalId": "80P",
            "chemicalName": "(4S,4aS,5aR,12aS)-4-(diethylamino)-3,10,12,12a-tetrahydroxy-1,11-dioxo-8-[(2S)-pyrrolidin-2-yl]-7-(trifluoromethyl)-1,4,4a,5,5a,6,11,12a-octahydrotetracene-2-carboxamide",
            "formula_weight": 0.58,
            "number_of_instances": 3,
            "pdbx_description": "(4S,4aS,5aR,12aS)-4-(diethylamino)-3,10,12,12a-tetrahydroxy-1,11-dioxo-8-[(2S)-pyrrolidin-2-yl]-7-(trifluoromethyl)-1,4,4a,5,5a,6,11,12a-octahydrotetracene-2-carboxamide",
        },
        [
            {
                "rcsb_id": "7RYG",
                "tax_node": {
                    "ncbi_tax_id": 480119,
                    "rank": "strain",
                    "scientific_name": "Acinetobacter baumannii AB0057",
                },
            },
            {
                "rcsb_id": "7RYF",
                "tax_node": {
                    "ncbi_tax_id": 480119,
                    "rank": "strain",
                    "scientific_name": "Acinetobacter baumannii AB0057",
                },
            },
            {
                "rcsb_id": "7RYH",
                "tax_node": {
                    "ncbi_tax_id": 480119,
                    "rank": "strain",
                    "scientific_name": "Acinetobacter baumannii AB0057",
                },
            },
        ],
    ),
    (
        {
            "InChI": "InChI=1S/C29H39N5O8/c1-28(2,3)31-11-17(35)32-15-10-16(33(4)5)13-8-12-9-14-21(34(6)7)24(38)20(27(30)41)26(40)29(14,42)25(39)18(12)23(37)19(13)22(15)36/h10,12,14,21,31,36,38-39,42H,8-9,11H2,1-7H3,(H2,30,41)(H,32,35)/p+2/t12-,14-,21-,29-/m0/s1",
            "InChIKey": "FPZLLRFZJZRHSY-HJYUBDRYSA-P",
            "SMILES": "CC(C)(C)NCC(=O)Nc1cc(c2c(c1O)C(=O)C3=C(C4(C(CC3C2)C(C(=C(C4=O)C(=O)N)O)[NH+](C)C)O)O)[NH+](C)C",
            "SMILES_stereo": "CC(C)(C)NCC(=O)Nc1cc(c2c(c1O)C(=O)C3=C([C@]4([C@@H](C[C@@H]3C2)[C@@H](C(=C(C4=O)C(=O)N)O)[NH+](C)C)O)O)[NH+](C)C",
            "chemicalId": "T1C",
            "chemicalName": "TIGECYCLINE",
            "formula_weight": 0.588,
            "number_of_instances": 2,
            "pdbx_description": "TIGECYCLINE",
        },
        [
            {
                "rcsb_id": "5J91",
                "tax_node": {
                    "ncbi_tax_id": 83333,
                    "rank": "strain",
                    "scientific_name": "Escherichia coli K-12",
                },
            },
            {
                "rcsb_id": "6YSI",
                "tax_node": {
                    "ncbi_tax_id": 575584,
                    "rank": "strain",
                    "scientific_name": "Acinetobacter baumannii ATCC 19606 = CIP "
                    "70.34 = JCM 6841",
                },
            },
            {
                "rcsb_id": "5J8A",
                "tax_node": {
                    "ncbi_tax_id": 83333,
                    "rank": "strain",
                    "scientific_name": "Escherichia coli K-12",
                },
            },
            {
                "rcsb_id": "4V9B",
                "tax_node": {
                    "ncbi_tax_id": 300852,
                    "rank": "strain",
                    "scientific_name": "Thermus thermophilus HB8",
                },
            },
            {
                "rcsb_id": "4YHH",
                "tax_node": {
                    "ncbi_tax_id": 300852,
                    "rank": "strain",
                    "scientific_name": "Thermus thermophilus HB8",
                },
            },
        ],
    ),
    (
        {
            "InChI": "InChI=1S/C29H30FN3O7/c1-32(2)22-18-8-15-7-14-6-13-4-3-12(9-33-10-16(30)11-33)5-17(13)23(34)19(14)24(35)20(15)26(37)29(18,40)27(38)21(25(22)36)28(31)39/h3-6,15-16,18,22,34,36-37,40H,7-11H2,1-2H3,(H2,31,39)/t15-,18-,22-,29-/m0/s1",
            "InChIKey": "DAUIQSOIFWOKJQ-MGVDVOGZSA-N",
            "SMILES": "CN(C)C1C2CC3Cc4cc5ccc(cc5c(c4C(=O)C3=C(C2(C(=O)C(=C1O)C(=O)N)O)O)O)CN6CC(C6)F",
            "SMILES_stereo": "CN(C)[C@H]1[C@@H]2C[C@@H]3Cc4cc5ccc(cc5c(c4C(=O)C3=C([C@@]2(C(=O)C(=C1O)C(=O)N)O)O)O)CN6CC(C6)F",
            "chemicalId": "P8F",
            "chemicalName": "Pentacycline",
            "formula_weight": 0.552,
            "number_of_instances": 2,
            "pdbx_description": "Pentacycline",
        },
        [
            {
                "rcsb_id": "8CGV",
                "tax_node": {
                    "ncbi_tax_id": 679895,
                    "rank": "no rank",
                    "scientific_name": "Escherichia coli BW25113",
                },
            }
        ],
    ),
    (
        {
            "InChI": "InChI=1S/C24H29N3O8/c1-26(2)18-13-8-11-7-12-10(9-27(3)35-4)5-6-14(28)16(12)19(29)15(11)21(31)24(13,34)22(32)17(20(18)30)23(25)33/h5-6,11,13,18,28,30-31,34H,7-9H2,1-4H3,(H2,25,33)/t11-,13-,18-,24-/m0/s1",
            "InChIKey": "PQJQFLNBMSCUSH-SBAJWEJLSA-N",
            "SMILES": "CN(C)C1C2CC3Cc4c(ccc(c4C(=O)C3=C(C2(C(=O)C(=C1O)C(=O)N)O)O)O)CN(C)OC",
            "SMILES_stereo": "CN(C)[C@H]1[C@@H]2C[C@@H]3Cc4c(ccc(c4C(=O)C3=C([C@@]2(C(=O)C(=C1O)C(=O)N)O)O)O)CN(C)OC",
            "chemicalId": "V7A",
            "chemicalName": "Sarecycline",
            "drugbank_description": "Sarecycline is a semi-synthetic derivative of "
            "tetracycline that was initially discovered by "
            "Paratek Pharmaceuticals from Boston, MA but then "
            "licensed to Warner Chilcott of Rockaway, NJ in "
            "July of 2007 [A40005]. After completing various "
            "phase-II and phase-III trials demonstrating its "
            "effectiveness in treating moderate to severe "
            "facial acne vulgaris [A39993, A39994] the US Food "
            "and Drug Administration approved Barcelona based "
            "Almirall, S.A.'s Seysara (sarecylcine) as a new "
            "first in class narrow spectrum tetracycline "
            "derived oral antibiotic for the treatment of "
            "inflammatory lesions of non-nodular moderate to "
            "severe acne vulgaris in patients nine years of age "
            "and older [L4814]. Seysara (sarecycline) was "
            "originally part of Allergan's US Medical "
            "Dermatology portfolio, before Almirall acquired "
            "the portfolio in the second half of 2018 as a "
            "means of consolidating and reinforcing the "
            "dermatology-focused pharmaceutical company's "
            "presence in the United States [L4815].  Acne "
            "vulgaris itself is a common chronic skin condition "
            "associated with the blockage and/or inflammation "
            "of hair follicles and their accompanying sebaceous "
            "glands [L4814]. The acne often presents physically "
            "as a mixture of non-inflammatory and inflammatory "
            "lesions mainly on the face but on the back and "
            "chest as well [L4814]. Based upon data from Global "
            "Burden of Disease studies, the acne vulgaris "
            "condition affects up to 85% of young adults aged "
            "12 to 25 years globally - with the possibility of "
            "permanent physical and mental scarring resulting "
            "from cases of severe acne [L4814].  Subsequently, "
            "while a number of first line tetracycline "
            "therapies like doxycycline and minocycline do "
            "exist for treating acne vulgaris, sarecycline "
            "presents a new and innovative therapy choice "
            "because it exhibits the necessary antibacterial "
            "activity against relevant pathogens that cause "
            "acne vulgaris but also possesses a low propensity "
            "for resistance development in such pathogens and a "
            "narrower, more specific spectrum of antibacterial "
            "activity, resulting in fewer off-target "
            "antibacterial effects on endogenous intestinal "
            "flora and consequently fewer resultant adverse "
            "effects associated with diarrhea, fungal "
            "overgrowth, etc.",
            "drugbank_id": "DB12035",
            "formula_weight": 0.488,
            "number_of_instances": 2,
            "pdbx_description": "Sarecycline",
        },
        [
            {
                "rcsb_id": "6XQD",
                "tax_node": {
                    "ncbi_tax_id": 300852,
                    "rank": "strain",
                    "scientific_name": "Thermus thermophilus HB8",
                },
            },
            {
                "rcsb_id": "6XQE",
                "tax_node": {
                    "ncbi_tax_id": 300852,
                    "rank": "strain",
                    "scientific_name": "Thermus thermophilus HB8",
                },
            },
        ],
    ),
]

composite_bsite = {
    "uL4": [
        {
            "source": 4,
            "mapping": "",
        },
        {
            "source": 5120,
            "mapping": "",
        },
    ]
}


CHEM_ID = NewType("CHEM_ID", str)
RCSB_ID = NewType("RCSB_ID", str)
target_rcsb_id = "7K00"
RADIUS = 10


def process_single_structure(
    structure_pair: tuple[RCSB_ID, CHEM_ID],
    target_rcsb_id: str = target_rcsb_id,
    radius: float = RADIUS,
) -> list[tuple[PolymerClass, PredictionTarget]]:
    """Process a single structure and return list of (polymer_class, target) pairs"""
    source_rcsb_id, chem_id = structure_pair
    print(f"Mapping {source_rcsb_id}/{chem_id}(source) into {target_rcsb_id} (target)")

    bsite_source = bsite_ligand(chem_id, source_rcsb_id, radius)
    bsite_target = bsite_transpose(source_rcsb_id, target_rcsb_id, bsite_source)

    results = []
    for chain in bsite_target.constituent_chains:
        results.append((chain.polymer_class, chain.target))
    return results


def prepare_mapping_sources(
    source_structures: list[tuple[RCSB_ID, CHEM_ID]], n_processes: int = 10
) -> dict[PolymerClass, list[PredictionTarget]]:

    per_class_registry: dict[PolymerClass, list[PredictionTarget]] = {}
    with Pool(processes=n_processes) as pool:
        all_results = pool.map(process_single_structure, source_structures)

    for structure_results in all_results:
        for polymer_class, target in structure_results:
            if polymer_class not in per_class_registry:
                per_class_registry[polymer_class] = [target]
            else:
                per_class_registry[polymer_class].append(target)

    return per_class_registry


sources = []
for record in tetracycline_structs:
    chem = record[0]
    structs = record[1]

    chemid = chem["chemicalId"]
    rcsb_id = structs[0]["rcsb_id"]
    sources.append((rcsb_id, chemid))

# Now there are definitely a fuckton of parameters and knobs one can tweak here
# - based on phylogeny(distance)
# - based on compounds and their function
# - based on the conservation of the binding site itself
# we are just going to roll with normalizing on the number of structures for now
# (i.e. a single residue hit in target is weighted via m/n where m is the number of source binding sites this residue figures in and n is the total number of binding sites)

registry = prepare_mapping_sources(sources[:3])


def compact_class(
    # target_chain: BiopythonChain,
    # target_rcsb_id: str,
    # polymer_class: PolymerClass,
    projections: list[PredictionTarget],
    number_of_sources: int,
) -> dict[int, float]:

    weights = {}

    for mapping in projections:
        residues = mapping.target_bound_residues
        for residue in residues:
            if residue.auth_seq_id not in weights:
                weights.update({residue.auth_seq_id: 1})
            else:
                weights[residue.auth_seq_id] += 1

    [ weights.update({x:y})  for x,y in list(map(lambda item: (item[0], item[1] / number_of_sources), weights.items())) ] 
    return weights


pprint(registry['16SrRNA'])
# target_chain_aaid = RibosomeOps(target_rcsb_id).get_poly_by_polyclass(PolymerClass("16SrRNA")).auth_asym_id
# target_chain = RibosomeOps(target_rcsb_id).biopython_structure()[0][target_chain_aaid]

compacted = compact_class(
    # BiopythonChain(target_chain),
    # target_rcsb_id,
    # PolymerClass("16SrRNA"),
   registry['16SrRNA'],
   3)
pprint(compacted)
# class MotifsMapManyToOne:
#     """Given multiple source sequences and residue ranges within them, project the ranges onto a single target sequence"""

#     target_chain : BiopythonChain
#     target_rcsb_id: RCSB_ID
#     polymer_class: PolymerClass
#     projections  : list[PredictionTarget]


#     def __init__(
#         self,
#       target_chain  : BiopythonChain,
#       target_rcsb_id: RCSB_ID,
#       motif_mappings: list[PredictionTarget],
#       polymer_class : PolymerClass
#     ) -> None       :
#         self.target_chain   = target_chain
#         self.target_rcsb_id = target_rcsb_id
#         self.polymer_class  = PolymerClass(polymer_class)
#         self.projections    = motif_mappings

#     def get_weights(self, number_of_instances) -> dict[int, list[RCSB_ID]]:

#         return {}

# return MotifsMapManyToOne(
#     target_chain,
#     target_rcsb_id,
#     polymer_class,


# )
# target chain

# ...
