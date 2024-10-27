from pprint import pprint
from typing import Callable, Literal
from neo4j import Driver, ManagedTransaction, Record, Transaction
from neo4j.graph import Node, Relationship
from neo4j import ManagedTransaction, Transaction
from ribctl.lib.schema.types_ribosome import NonpolymericLigand, RibosomeStructure, RibosomeStructureMetadata

def node__ligand(
    _ligand: NonpolymericLigand,
) -> Callable[[Transaction | ManagedTransaction], Node]:
    def _(tx: Transaction | ManagedTransaction):

            # Prepare the properties dictionary
            properties = {
                "chemicalId"          : _ligand.chemicalId,
                "chemicalName"        : _ligand.chemicalName,
                "formula_weight"      : _ligand.formula_weight,
                "pdbx_description"    : _ligand.pdbx_description,
                "number_of_instances" : _ligand.number_of_instances,
                "drugbank_id"         : None,
                "drugbank_description": None,
                "SMILES"              : _ligand.SMILES,
                "SMILES_stereo"       : _ligand.SMILES_stereo,
                "InChI"               : _ligand.InChI,
                "InChIKey"            : _ligand.InChIKey,
            }

            # this is just to get the props from multiple nullable nestings
            if _ligand.nonpolymer_comp:
                    if _ligand.nonpolymer_comp.drugbank:
                        if _ligand.nonpolymer_comp.drugbank.drugbank_container_identifiers:
                            properties["drugbank_id"] = _ligand.nonpolymer_comp.drugbank.drugbank_container_identifiers.drugbank_id
                        if _ligand.nonpolymer_comp.drugbank.drugbank_info:
                            properties["drugbank_description"] = _ligand.nonpolymer_comp.drugbank.drugbank_info.description

            

            properties = {k: v for k, v in properties.items() if v is not None}
            query = """
            MERGE (ligand:Ligand {chemicalId: $chemicalId})
            SET ligand += $properties
            RETURN ligand
            """
            result = tx.run(query, {"chemicalId": _ligand.chemicalId, "properties": properties})
            return result.single(strict=True)["ligand"]

    return _

# Transaction
def link__ligand_to_struct(
    prot: Node, parent_rcsb_id: str
) -> Callable[[Transaction | ManagedTransaction], list[list[Node | Relationship]]]:
    parent_rcsb_id = parent_rcsb_id.upper()

    def _(tx: Transaction | ManagedTransaction):
        return tx.run(
            """//
  match (ligand:Ligand) where ELEMENTID(ligand)=$ELEM_ID
  match (struct:RibosomeStructure {rcsb_id: $PARENT})
  merge (ligand)<-[contains:contains]-(struct)
  return struct, ligand, contains
""",
            {"ELEM_ID": prot.element_id, "PARENT": parent_rcsb_id},
        ).values("struct", "ligand", "contains")

    return _