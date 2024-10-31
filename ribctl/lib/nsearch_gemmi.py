import gemmi
import sys

def get_ligand_neighbors(structure_path, ligand_id, radius=10.0):
    # 1. Open the mmCIF structure
    structure = gemmi.read_structure(structure_path)

    # 2. Get all residues of the specified ligand
    ligand_residues = []
    for model in structure:
        print("got model")
        for chain in model:
            for residue in chain:
                if residue.name == ligand_id:
                    ligand_residues.append(residue)

    print(ligand_residues)
    if not ligand_residues:
        print(f"No residues found for ligand {ligand_id}")
        return []

    print(structure[0])
    ns = gemmi.NeighborSearch(structure[0], structure.cell, 5).populate()
    print(ns)

    neighboring_residues = set()
    for ligand_residue in ligand_residues:
        # pprint(dir(ligand_residue))
        # print(ligand_residue.seqid)
        # print(ligand_residue.label_seq)
        for atom in ligand_residue:
            # print(dir(atom))
            neighbors = ns.find_neighbors(atom, radius, radius)
            for neighbor in neighbors:
                neighbor_res = neighbor.to_cra(structure[0]).residue
                if neighbor_res.name != ligand_id:
                    neighboring_residues.add((neighbor_res.name, neighbor_res.seqid.num, neighbor_res.chain.name))

    return list(neighboring_residues)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <path_to_mmcif> <ligand_id>")
        sys.exit(1)

    structure_path = sys.argv[1]
    ligand_id = sys.argv[2]

    neighbors = get_ligand_neighbors(structure_path, ligand_id)

    print(f"Neighbors of {ligand_id} within 10 Angstroms:")
    for residue in neighbors:
        print(f"Residue: {residue[0]}, Seq ID: {residue[1]}, Chain: {residue[2]}")
    print(f"Total neighboring residues: {len(neighbors)}")