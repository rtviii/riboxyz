import numpy as np
import os
cmd.load("/Volumes/LaCie/project2/structure/3CF5.cif")
cmd.select('ul4', 'chain C')
cmd.select('ul22', 'chain P')
cmd.select('ul23', 'chain Q')
cmd.create('ul4_group', 'ul4')
cmd.show('spheres', 'ul4_group')
cmd.color('red', 'ul4_group')
cmd.create('ul22_group', 'ul22')
cmd.show('spheres', 'ul22_group')
cmd.color('yellow', 'ul22_group')
cmd.create('ul23_group', 'ul23')
cmd.show('spheres', 'ul23_group')
cmd.color('blue', 'ul23_group')
cmd.select('residue_l23', 'resi 14 and chain Q')
cmd.create('exit_conserved', 'residue_l23')
coords = cmd.get_coords('residue_l23')
file_path = '/Volumes/LaCie/project2/exit_residue/l23_coords_3CF5.txt'
if os.path.exists(file_path):
   os.remove(file_path)
np.savetxt('/Volumes/LaCie/project2/exit_residue/l23_coords_3CF5.txt', coords[0])
cutoff = 15.0
cmd.select("close_atoms", f"ul4_group within {cutoff} of ul22_group")
cmd.create('closest', 'close_atoms')
atom_list = cmd.get_model("closest").atom
atoms = cmd.get_model("closest").atom
last_atom = "closest and id " + str(atom_list[-1].id)
cmd.select("last_atom", last_atom)
cmd.create('last_a', 'last_atom')
cmd.select("edited", "br. 3CF5 within 80 of last_a")
cmd.create("3CF5_edited", "edited")
cmd.show('spheres', '3CF5_edited')
cmd.save("/Volumes/LaCie/project2/structure/edited/cif/3CF5_edited.cif", "3CF5_edited")
