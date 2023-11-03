import numpy as np
import os
cmd.load("/Volumes/LaCie/project2/euk/structure/4U4R.cif")
cmd.select('ul4', 'chain L4')
cmd.select('ul22', 'chain M7')
cmd.select('el39', 'chain O9')
cmd.create('ul4_group', 'ul4')
cmd.show('spheres', 'ul4_group')
cmd.color('red', 'ul4_group')
cmd.create('ul22_group', 'ul22')
cmd.show('spheres', 'ul22_group')
cmd.color('yellow', 'ul22_group')
cmd.create('el39_group', 'el39')
cmd.show('spheres', 'el39_group')
cmd.color('blue', 'el39_group')
cmd.select('residue_l39', 'resi 37 and chain O9')
cmd.create('exit_conserved', 'residue_l39')
coords = cmd.get_coords('residue_l39')
file_path = '/Volumes/LaCie/project2/euk/exit_residue/l39_coords_4U4R.txt'
if os.path.exists(file_path):
   os.remove(file_path)
np.savetxt('/Volumes/LaCie/project2/euk/exit_residue/l39_coords_4U4R.txt', coords[0])
cutoff = 15.0
cmd.select("close_atoms", f"ul4_group within {cutoff} of ul22_group")
cmd.create('closest', 'close_atoms')
atom_list = cmd.get_model("closest").atom
atoms = cmd.get_model("closest").atom
last_atom = "closest and id " + str(atom_list[-1].id)
cmd.select("last_atom", last_atom)
cmd.create('last_a', 'last_atom')
cmd.select("edited", "br. 4U4R within 80 of last_a")
cmd.create("4U4R_edited", "edited")
cmd.show('spheres', '4U4R_edited')
cmd.save("/Volumes/LaCie/project2/euk/structure/edited/pdb/4U4R_edited.pdb", "4U4R_edited")
