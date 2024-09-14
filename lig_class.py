from rdkit import Chem
from rdkit.Chem import AllChem
import requests
import json

from neo4j_ribosome.db_lib_reader import Neo4jReader


def get_compound_class(smiles):
    # Generate InChI from SMILES using RDKit
    mol = Chem.MolFromSmiles(smiles)
    inchi = Chem.MolToInchi(mol)

    # Use ClassyFire API to get classification
    url = "http://classyfire.wishartlab.com/queries"
    payload = {"query_input": inchi, "query_type": "INCHI"}
    headers = {"Content-Type": "application/json"}

    response = requests.post(url, data=json.dumps(payload), headers=headers)

    data = response.json()

    # Get the query ID and wait for results
    query_id = data["id"]
    result_url = f"http://classyfire.wishartlab.com/queries/{query_id}.json"

    while True:
        response = requests.get(result_url)
        data = response.json()
        if data["classification_status"] == "Done":
            break

    classification = data["entities"][0]["direct_parent"]
    return classification

# Example usage
smiles = "CC1[C@H]([C@@H]([C@H]([C@@H]([C@H]1N)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)N)O[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O[C@@H]4[C@@H]([C@H]([C@@H]([C@@H](O4)CN)O)O)N)O)O)N"  # SMILES for Paromomycin
compound_class = get_compound_class(smiles)
print(f"The class of the compound is: {compound_class}")


def collect_all_ligands_with_smiles():
    reader = Neo4jReader()
    # ligs = reader.list_ligands()



# curl  -H "Accept: application/json" -H "Content-type: application/json" -d '{"label":"test", "query_input":"MOL1\\tCCCOCC\\nMOL2\\tCOCC=CCCC", "query_type":"STRUCTURE"}' -X POST http://classyfire.wishartlab.com/queries

# curl  -H "Accept: application/json" -H "Content-type: application/json" -d '{"label":"test", "query_input":"PAR\tCC1[C@H]([C@@H]([C@H]([C@@H]([C@H]1N)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)N)O[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O[C@@H]4[C@@H]([C@H]([C@@H]([C@@H](O4)CN)O)O)N)O)O)N\nKSG\tC[C@@H]1[C@H](C[C@@H]([C@H](O1)OC2[C@@H]([C@H](C([C@@H]([C@@H]2O)O)O)O)O)N)N=C(C(=O)O)N", "query_type":"STRUCTURE"}' -X POST http://classyfire.wishartlab.com/queries


# {
#     "id": 11703588,
#     "label": "test",
#     "classification_status": "Done",
#     "number_of_elements": 2,
#     "number_of_pages": 1,
#     "invalid_entities": [],
#     "entities": [
#         {
#             "identifier": "PAR",
#             "smiles": "CC1[C@@H](N)[C@H](O)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@H]3O[C@@H](CN)[C@@H](O)[C@H](O)[C@H]3N)[C@H]2O)[C@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2N)[C@H]1N",
#             "inchikey": "InChIKey=YTYIXKOKVUNKNL-OQYBOFSWSA-N",
#             "kingdom": {
#                 "name": "Organic compounds",
#                 "description": "Compounds that contain at least one carbon atom, excluding isocyanide/cyanide and their non-hydrocarbyl derivatives, thiophosgene, carbon diselenide, carbon monosulfide, carbon disulfide, carbon subsulfide, carbon monoxide, carbon dioxide, Carbon suboxide, and dicarbon monoxide.",
#                 "chemont_id": "CHEMONTID:0000000",
#                 "url": "http://classyfire.wishartlab.com/tax_nodes/C0000000",
#             },
#             "superclass": {
#                 "name": "Lipids and lipid-like molecules",
#                 "description": "Fatty acids and their derivatives, and substances related biosynthetically or functionally to these compounds.",
#                 "chemont_id": "CHEMONTID:0000012",
#                 "url": "http://classyfire.wishartlab.com/tax_nodes/C0000012",
#             },
#             "class": {
#                 "name": "Fatty Acyls",
#                 "description": "Organic molecules synthesized by chain elongation of an acetyl-CoA primer with malonyl-CoA (or methylmalonyl-CoA) groups that might contain a cyclic functionality and/or are substituted with heteroatoms.",
#                 "chemont_id": "CHEMONTID:0003909",
#                 "url": "http://classyfire.wishartlab.com/tax_nodes/C0003909",
#             },
#             "subclass": {
#                 "name": "Fatty acyl glycosides",
#                 "description": "Compounds containing fatty acyl chain linked to a carbohydrate moiety through a glycosidic bond. Fatty acyl glycosides are composed of a glycosyl moiety (one or several units) linked to one hydroxyl group of a fatty alcohol or of a phosphorylated alcohol (phosphoprenols), a hydroxy fatty acid or to one carboxyl group of a fatty acid (ester linkage) or to an amino alcohol.",
#                 "chemont_id": "CHEMONTID:0001766",
#                 "url": "http://classyfire.wishartlab.com/tax_nodes/C0001766",
#             },
#             "intermediate_nodes": [],
#             "direct_parent": {
#                 "name": "Fatty acyl glycosides of mono- and disaccharides",
#                 "description": "Compounds composed of a mono- or disaccharide moiety linked to one hydroxyl group of a fatty alcohol or of a phosphorylated alcohol (phosphoprenols), a hydroxy fatty acid or to one carboxyl group of a fatty acid (ester linkage) or to an amino alcohol.",
#                 "chemont_id": "CHEMONTID:0003861",
#                 "url": "http://classyfire.wishartlab.com/tax_nodes/C0003861",
#             },
#             "alternative_parents": [
#                 {
#                     "name": "Alkyl glycosides",
#                     "description": "Lipids containing a glycosyl moiety (one or several units) linked to the hydroxyl group of a fatty alcohol.",
#                     "chemont_id": "CHEMONTID:0002125",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0002125",
#                 },
#                 {
#                     "name": "Cyclohexanols",
#                     "description": "Compounds containing an alcohol group attached to a cyclohexane ring.",
#                     "chemont_id": "CHEMONTID:0002647",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0002647",
#                 },
#                 {
#                     "name": "Oxanes",
#                     "description": "Compounds containing an oxane (tetrahydropyran) ring, which is a six-member saturated aliphatic heterocycle with one oxygen atom and five carbon atoms.",
#                     "chemont_id": "CHEMONTID:0002012",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0002012",
#                 },
#                 {
#                     "name": "Monosaccharides",
#                     "description": "Compounds containing one carbohydrate unit not glycosidically linked to another such unit, and no set of two or more glycosidically linked carbohydrate units. Monosaccharides have the general formula CnH2nOn.",
#                     "chemont_id": "CHEMONTID:0001540",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0001540",
#                 },
#                 {
#                     "name": "Cyclitols and derivatives",
#                     "description": "Compounds containing a cycloalkane moiety with one hydroxyl group on each of three or more ring atoms. These of also include derivatives where the hydrogen atom of one or more of the hydroxyl groups is replaced with another atom.",
#                     "chemont_id": "CHEMONTID:0002509",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0002509",
#                 },
#                 {
#                     "name": "Tetrahydrofurans",
#                     "description": "Heterocyclic compounds containing a saturated, aliphatic, five-membered ring where a carbon is replaced by an oxygen.",
#                     "chemont_id": "CHEMONTID:0002648",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0002648",
#                 },
#                 {
#                     "name": "Oxacyclic compounds",
#                     "description": "Organic compounds containing an heterocycle with at least one oxygen atom linked to a ring carbon.",
#                     "chemont_id": "CHEMONTID:0004140",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0004140",
#                 },
#                 {
#                     "name": "Organopnictogen compounds",
#                     "description": "Compounds containing a bond between carbon a pnictogen atom. Pnictogens are p-block element atoms that are in the group 15 of the periodic table.",
#                     "chemont_id": "CHEMONTID:0004557",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0004557",
#                 },
#                 {
#                     "name": "Monoalkylamines",
#                     "description": "Organic compounds containing an primary aliphatic amine group.",
#                     "chemont_id": "CHEMONTID:0000469",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0000469",
#                 },
#                 {
#                     "name": "Hydrocarbon derivatives",
#                     "description": "Derivatives of hydrocarbons obtained by substituting one or more carbon atoms by an heteroatom. They contain at least one carbon atom and heteroatom.",
#                     "chemont_id": "CHEMONTID:0004150",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0004150",
#                 },
#             ],
#             "molecular_framework": "Aliphatic heteromonocyclic compounds",
#             "substituents": [
#                 "Fatty acyl glycoside of mono- or disaccharide",
#                 "Alkyl glycoside",
#                 "Cyclohexanol",
#                 "Oxane",
#                 "Monosaccharide",
#                 "Cyclitol or derivatives",
#                 "Tetrahydrofuran",
#                 "Cyclic alcohol",
#                 "Secondary alcohol",
#                 "Oxacycle",
#                 "Organoheterocyclic compound",
#                 "Organic nitrogen compound",
#                 "Organic oxygen compound",
#                 "Organopnictogen compound",
#                 "Hydrocarbon derivative",
#                 "Organooxygen compound",
#                 "Organonitrogen compound",
#                 "Primary aliphatic amine",
#                 "Alcohol",
#                 "Aliphatic heteromonocyclic compound",
#             ],
#             "description": "This compound belongs to the class of organic compounds known as fatty acyl glycosides of mono- and disaccharides. These are compounds composed of a mono- or disaccharide moiety linked to one hydroxyl group of a fatty alcohol or of a phosphorylated alcohol (phosphoprenols), a hydroxy fatty acid or to one carboxyl group of a fatty acid (ester linkage) or to an amino alcohol.",
#             "external_descriptors": [],
#             "ancestors": [
#                 "Alcohols and polyols",
#                 "Alkyl glycosides",
#                 "Amines",
#                 "Carbohydrates and carbohydrate conjugates",
#                 "Chemical entities",
#                 "Cyclic alcohols and derivatives",
#                 "Cyclitols and derivatives",
#                 "Cyclohexanols",
#                 "Fatty Acyls",
#                 "Fatty acyl glycosides",
#                 "Fatty acyl glycosides of mono- and disaccharides",
#                 "Hydrocarbon derivatives",
#                 "Lipids and lipid-like molecules",
#                 "Monoalkylamines",
#                 "Monosaccharides",
#                 "Organic compounds",
#                 "Organic nitrogen compounds",
#                 "Organic oxygen compounds",
#                 "Organoheterocyclic compounds",
#                 "Organonitrogen compounds",
#                 "Organooxygen compounds",
#                 "Organopnictogen compounds",
#                 "Oxacyclic compounds",
#                 "Oxanes",
#                 "Primary amines",
#                 "Secondary alcohols",
#                 "Tetrahydrofurans",
#             ],
#             "predicted_chebi_terms": [],
#             "predicted_lipidmaps_terms": [],
#             "classification_version": "2.1",
#         },
#         {
#             "identifier": "KSG",
#             "smiles": "C[C@H]1O[C@H](OC2[C@@H](O)[C@@H](O)C(O)[C@H](O)[C@H]2O)[C@@H](N)C[C@@H]1\\N=C(/N)C(O)=O",
#             "inchikey": "InChIKey=PVTHJAPFENJVNC-UQTMRZPGSA-N",
#             "kingdom": {
#                 "name": "Organic compounds",
#                 "description": "Compounds that contain at least one carbon atom, excluding isocyanide/cyanide and their non-hydrocarbyl derivatives, thiophosgene, carbon diselenide, carbon monosulfide, carbon disulfide, carbon subsulfide, carbon monoxide, carbon dioxide, Carbon suboxide, and dicarbon monoxide.",
#                 "chemont_id": "CHEMONTID:0000000",
#                 "url": "http://classyfire.wishartlab.com/tax_nodes/C0000000",
#             },
#             "superclass": {
#                 "name": "Organic oxygen compounds",
#                 "description": "Organic compounds that contain one or more oxygen atoms.",
#                 "chemont_id": "CHEMONTID:0004603",
#                 "url": "http://classyfire.wishartlab.com/tax_nodes/C0004603",
#             },
#             "class": {
#                 "name": "Organooxygen compounds",
#                 "description": "Organic compounds containing a bond between a carbon atom and an oxygen atom.",
#                 "chemont_id": "CHEMONTID:0000323",
#                 "url": "http://classyfire.wishartlab.com/tax_nodes/C0000323",
#             },
#             "subclass": {
#                 "name": "Carbohydrates and carbohydrate conjugates",
#                 "description": "Monosaccharides, disaccharides, oligosaccharides, polysaccharides, and their derivatives.",
#                 "chemont_id": "CHEMONTID:0000011",
#                 "url": "http://classyfire.wishartlab.com/tax_nodes/C0000011",
#             },
#             "intermediate_nodes": [
#                 {
#                     "name": "Aminosaccharides",
#                     "description": "Saccharides containing a sugar unit that bears an amino group.",
#                     "chemont_id": "CHEMONTID:0003305",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0003305",
#                 },
#                 {
#                     "name": "Aminoglycosides",
#                     "description": "Molecules or a portion of a molecule composed of amino-modified sugars.",
#                     "chemont_id": "CHEMONTID:0000282",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0000282",
#                 },
#             ],
#             "direct_parent": {
#                 "name": "Aminocyclitol glycosides",
#                 "description": "Organic compounds containing an amicocyclitol moiety glycosidically linked to a carbohydrate moiety. There are two major classes of aminoglycosides containing a 2-streptamine core. They are called 4,5- and 4,6-disubstituted 2-deoxystreptamines.",
#                 "chemont_id": "CHEMONTID:0001675",
#                 "url": "http://classyfire.wishartlab.com/tax_nodes/C0001675",
#             },
#             "alternative_parents": [
#                 {
#                     "name": "Alpha amino acids and derivatives",
#                     "description": "Amino acids in which the amino group is attached to the carbon atom immediately adjacent to the carboxylate group (alpha carbon), or a derivative thereof.",
#                     "chemont_id": "CHEMONTID:0000060",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0000060",
#                 },
#                 {
#                     "name": "Cyclohexanols",
#                     "description": "Compounds containing an alcohol group attached to a cyclohexane ring.",
#                     "chemont_id": "CHEMONTID:0002647",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0002647",
#                 },
#                 {
#                     "name": "Oxanes",
#                     "description": "Compounds containing an oxane (tetrahydropyran) ring, which is a six-member saturated aliphatic heterocycle with one oxygen atom and five carbon atoms.",
#                     "chemont_id": "CHEMONTID:0002012",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0002012",
#                 },
#                 {
#                     "name": "Monosaccharides",
#                     "description": "Compounds containing one carbohydrate unit not glycosidically linked to another such unit, and no set of two or more glycosidically linked carbohydrate units. Monosaccharides have the general formula CnH2nOn.",
#                     "chemont_id": "CHEMONTID:0001540",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0001540",
#                 },
#                 {
#                     "name": "Cyclitols and derivatives",
#                     "description": "Compounds containing a cycloalkane moiety with one hydroxyl group on each of three or more ring atoms. These of also include derivatives where the hydrogen atom of one or more of the hydroxyl groups is replaced with another atom.",
#                     "chemont_id": "CHEMONTID:0002509",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0002509",
#                 },
#                 {
#                     "name": "Amino acids",
#                     "description": "Organic compounds that contain at least one carboxyl group and one amino group.",
#                     "chemont_id": "CHEMONTID:0004176",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0004176",
#                 },
#                 {
#                     "name": "Propargyl-type 1,3-dipolar organic compounds",
#                     "description": "Organic 1,3-dipolar compounds with the general structure  X#N+-Z- \u003c-\u003e X-=N+=Z \u003c-\u003e X-=N-Z+ \u003c-\u003e X-N=Z (X = C or O, Z = C, N, or O).",
#                     "chemont_id": "CHEMONTID:0003633",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0003633",
#                 },
#                 {
#                     "name": "Polyols",
#                     "description": "Organic compounds containing more than one hydroxyl groups.",
#                     "chemont_id": "CHEMONTID:0002286",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0002286",
#                 },
#                 {
#                     "name": "Oxacyclic compounds",
#                     "description": "Organic compounds containing an heterocycle with at least one oxygen atom linked to a ring carbon.",
#                     "chemont_id": "CHEMONTID:0004140",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0004140",
#                 },
#                 {
#                     "name": "Acetals",
#                     "description": "Compounds having the structure R2C(OR')2 ( R' not Hydrogen) and thus diethers of geminal diols. Originally, the term was confined to derivatives of aldehydes (one R = H), but it now applies equally to derivatives of ketones (neither R = H ). Mixed acetals have different R' groups.",
#                     "chemont_id": "CHEMONTID:0001656",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0001656",
#                 },
#                 {
#                     "name": "Monocarboxylic acids and derivatives",
#                     "description": "Carboxylic acids containing exactly one carboxyl groups.",
#                     "chemont_id": "CHEMONTID:0001137",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0001137",
#                 },
#                 {
#                     "name": "Carboxamidines",
#                     "description": "Carboxylic acid derivatives containing the amidine group.",
#                     "chemont_id": "CHEMONTID:0001045",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0001045",
#                 },
#                 {
#                     "name": "Carboxylic acids",
#                     "description": "Compounds containing a carboxylic acid group with the formula -C(=O)OH.",
#                     "chemont_id": "CHEMONTID:0001205",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0001205",
#                 },
#                 {
#                     "name": "Carboximidamides",
#                     "description": "Organonitrogen compounds with the general formula RN=CN(R')(R''), (R,R',R'' = H or organyl group).",
#                     "chemont_id": "CHEMONTID:0003152",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0003152",
#                 },
#                 {
#                     "name": "Hydrocarbon derivatives",
#                     "description": "Derivatives of hydrocarbons obtained by substituting one or more carbon atoms by an heteroatom. They contain at least one carbon atom and heteroatom.",
#                     "chemont_id": "CHEMONTID:0004150",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0004150",
#                 },
#                 {
#                     "name": "Carbonyl compounds",
#                     "description": "Organic compounds containing a carbonyl group, with the general structure RC(=O)R', where R=organyl, R'=H, N, O, organyl group or halide group.",
#                     "chemont_id": "CHEMONTID:0001831",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0001831",
#                 },
#                 {
#                     "name": "Monoalkylamines",
#                     "description": "Organic compounds containing an primary aliphatic amine group.",
#                     "chemont_id": "CHEMONTID:0000469",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0000469",
#                 },
#                 {
#                     "name": "Organic oxides",
#                     "description": "Organic compounds containing an oxide group.",
#                     "chemont_id": "CHEMONTID:0003940",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0003940",
#                 },
#                 {
#                     "name": "Organopnictogen compounds",
#                     "description": "Compounds containing a bond between carbon a pnictogen atom. Pnictogens are p-block element atoms that are in the group 15 of the periodic table.",
#                     "chemont_id": "CHEMONTID:0004557",
#                     "url": "http://classyfire.wishartlab.com/tax_nodes/C0004557",
#                 },
#             ],
#             "molecular_framework": "Aliphatic heteromonocyclic compounds",
#             "substituents": [
#                 "Amino cyclitol glycoside",
#                 "Alpha-amino acid or derivatives",
#                 "Cyclohexanol",
#                 "Cyclitol or derivatives",
#                 "Monosaccharide",
#                 "Oxane",
#                 "Cyclic alcohol",
#                 "Amino acid or derivatives",
#                 "Amino acid",
#                 "Secondary alcohol",
#                 "Acetal",
#                 "Amidine",
#                 "Carboxylic acid amidine",
#                 "Carboxylic acid derivative",
#                 "Carboxylic acid",
#                 "Oxacycle",
#                 "Organoheterocyclic compound",
#                 "Monocarboxylic acid or derivatives",
#                 "Organic 1,3-dipolar compound",
#                 "Propargyl-type 1,3-dipolar organic compound",
#                 "Polyol",
#                 "Carboximidamide",
#                 "Alcohol",
#                 "Amine",
#                 "Organopnictogen compound",
#                 "Primary aliphatic amine",
#                 "Organic oxide",
#                 "Hydrocarbon derivative",
#                 "Organonitrogen compound",
#                 "Carbonyl group",
#                 "Primary amine",
#                 "Organic nitrogen compound",
#                 "Aliphatic heteromonocyclic compound",
#             ],
#             "description": "This compound belongs to the class of organic compounds known as aminocyclitol glycosides. These are organic compounds containing an amicocyclitol moiety glycosidically linked to a carbohydrate moiety. There are two major classes of aminoglycosides containing a 2-streptamine core. They are called 4,5- and 4,6-disubstituted 2-deoxystreptamines.",
#             "external_descriptors": [],
#             "ancestors": [
#                 "Acetals",
#                 "Alcohols and polyols",
#                 "Alpha amino acids and derivatives",
#                 "Amidines",
#                 "Amines",
#                 "Amino acids",
#                 "Amino acids and derivatives",
#                 "Amino acids, peptides, and analogues",
#                 "Aminocyclitol glycosides",
#                 "Aminoglycosides",
#                 "Aminosaccharides",
#                 "Carbohydrates and carbohydrate conjugates",
#                 "Carbonyl compounds",
#                 "Carboxamidines",
#                 "Carboximidamides",
#                 "Carboxylic acids",
#                 "Carboxylic acids and derivatives",
#                 "Chemical entities",
#                 "Cyclic alcohols and derivatives",
#                 "Cyclitols and derivatives",
#                 "Cyclohexanols",
#                 "Ethers",
#                 "Hydrocarbon derivatives",
#                 "Monoalkylamines",
#                 "Monocarboxylic acids and derivatives",
#                 "Monosaccharides",
#                 "Organic 1,3-dipolar compounds",
#                 "Organic acids and derivatives",
#                 "Organic compounds",
#                 "Organic nitrogen compounds",
#                 "Organic oxides",
#                 "Organic oxygen compounds",
#                 "Organoheterocyclic compounds",
#                 "Organonitrogen compounds",
#                 "Organooxygen compounds",
#                 "Organopnictogen compounds",
#                 "Oxacyclic compounds",
#                 "Oxanes",
#                 "Polyols",
#                 "Primary amines",
#                 "Propargyl-type 1,3-dipolar organic compounds",
#                 "Secondary alcohols",
#             ],
#             "predicted_chebi_terms": [
#                 "organonitrogen compound (CHEBI:35352)",
#                 "organooxygen compound (CHEBI:36963)",
#                 "cyclohexanols (CHEBI:23480)",
#                 "oxanes (CHEBI:46942)",
#                 "monosaccharide (CHEBI:35381)",
#                 "amino acid (CHEBI:33709)",
#                 "dipolar compound (CHEBI:51151)",
#                 "polyol (CHEBI:26191)",
#                 "oxacycle (CHEBI:38104)",
#                 "acetal (CHEBI:59769)",
#                 "carbonyl compound (CHEBI:36586)",
#                 "carboxamidine (CHEBI:35359)",
#                 "carboxylic acid (CHEBI:33575)",
#                 "carboxylic acid anion (CHEBI:29067)",
#                 "organic molecule (CHEBI:72695)",
#                 "alkylamine (CHEBI:13759)",
#                 "organic oxide (CHEBI:25701)",
#                 "pnictogen molecular entity (CHEBI:33302)",
#                 "organic molecular entity (CHEBI:50860)",
#                 "glycoside (CHEBI:24400)",
#                 "amino cyclitol (CHEBI:61689)",
#                 "chemical entity (CHEBI:24431)",
#                 "peptide (CHEBI:16670)",
#                 "oxygen molecular entity (CHEBI:25806)",
#                 "organic hydroxy compound (CHEBI:33822)",
#                 "alcohol (CHEBI:30879)",
#                 "secondary alcohol (CHEBI:35681)",
#                 "organic heterocyclic compound (CHEBI:24532)",
#                 "carbohydrates and carbohydrate derivatives (CHEBI:78616)",
#                 "ether (CHEBI:25698)",
#                 "nitrogen molecular entity (CHEBI:51143)",
#                 "amidine (CHEBI:2634)",
#                 "amine (CHEBI:32952)",
#                 "primary amine (CHEBI:32877)",
#                 "amino sugar (CHEBI:28963)",
#                 "aminoglycoside (CHEBI:47779)",
#             ],
#             "predicted_lipidmaps_terms": [],
#             "classification_version": "2.1",
#         },
#     ],
# }
