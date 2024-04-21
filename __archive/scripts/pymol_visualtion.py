import json
import os
from pprint import pprint
import random
import sys
from typing import List
import chempy
import pandas as pd
from pymol import cmd
from ribctl.etl.ribosome_assets import RibosomeAssets
from mesh_generation.paths import TUNNEL_PATH, poisson_recon_ascii_path
sys.path.append('/home/rtviii/dev/riboxyz')       #! hack until ribctl is a separate pypi project (after that just pip install ribctl)
from ribctl.lib.tunnel import ptc_residues_calculate_midpoint, ptc_resdiues_get
from ribctl import EXIT_TUNNEL_WORK, RIBETL_DATA
from masif_plugin import load_ply

colormap__RNA = {
    "23SrRNA"  : "gray70",    # plants
    "25SrRNA"  : "gray70",    # plants
    "28SrRNA"  : 'gray70',    # eukaryotic
    "5.8SrRNA" : 'lightblue', # eukaryotic
    "5SrRNA"   : "palegreen",
    "16SrRNA"  : "palegreen",
    "18SrRNA"  : "palegreen", # eukaryotic
    "mt12SrRNA": "palegreen", # mitochondrial
    "mt16SrRNA": 'gray70',    # mitochondrial,
    "tRNA"     : "lightblue",

          }
colormap__LSU_Proteins = {
                          "uL1"  : "lawrencium",
                          "uL2"  : "lead",
                          "uL3"  : "lithium",
                          "uL4"  : "lutetium",
                          "uL5"  : "magnesium",
                          "uL6"  : "manganese",
                          "eL6"  : "meitnerium",
                          "eL8"  : "mendelevium",
                          "bL9"  : "mercury",
                          "uL10" : "molybdenum",
                          "uL11" : "neodymium",
                          "bL12" : "neon",
                          "uL13" : "neptunium",
                          "eL13" : "nickel",
                          "uL14" : "niobium",
                          "eL14" : "nitrogen",
                          "uL15" : "nobelium",
                          "eL15" : "osmium",
                          "uL16" : "oxygen",
                          "bL17" : "palladium",
                          "uL18" : "phosphorus",
                          "eL18" : "platinum",
                          "bL19" : "plutonium",
                          "eL19" : "polonium",
                          "bL20" : "potassium",
                          "eL20" : "praseodymium",
                          "bL21" : "promethium",
                          "eL21" : "protactinium",
                          "uL22" : "radium",
                          "eL22" : "radon",
                          "uL23" : "rhenium",
                          "uL24" : "rhodium",
                          "eL24" : "rubidium",
                          "bL25" : "ruthenium",
                          "bL27" : "rutherfordium",
                          "eL27" : "samarium",
                          "bL28" : "scandium",
                          "eL28" : "seaborgium",
                          "uL29" : "selenium",
                          "eL29" : "silicon",
                          "uL30" : "silver",
                          "eL30" : "sodium",
                          "bL31" : "strontium",
                          "eL31" : "sulfur",
                          "bL32" : "tantalum",
                          "eL32" : "technetium",
                          "bL33" : "tellurium",
                          "eL33" : "terbium",
                          "bL34" : "thallium",
                          "eL34" : "thorium",
                          "bL35" : "thulium",
                          "bL36" : "tin",
                          "eL36" : "titanium",
                          "eL37" : "tungsten",
                          "eL38" : "uranium",
                          "eL39" : "vanadium",
                          "eL40" : "xenon",
                          "eL41" : "ytterbium",
                          "eL42" : "yttrium",
                          "eL43" : "zinc",
                          "P1/P2": "zirconium"
                          
                          }
colormap__SSU_Proteins = {
                     "bS1" : "actinium",
                     "eS1" : "aluminum",
                     "uS2" : "americium",
                     "uS3" : "antimony",
                     "uS4" : "argon",
                     "eS4" : "arsenic",
                     "uS5" : "astatine",
                     "bS6" : "barium",
                     "eS6" : "berkelium",
                     "uS7" : "beryllium",
                     "eS7" : "bismuth",
                     "uS8" : "bohrium",
                     "eS8" : "boron",
                     "uS9" : "bromine",
                     "uS10": "cadmium",
                     "eS10": "calcium",
                     "uS11": "californium",
                     "uS12": "carbon",
                     "eS12": "cerium",
                     "uS13": "cesium",
                     "uS14": "chlorine",
                     "uS15": "chromium",
                     "bS16": "cobalt",
                     "uS17": "copper",
                     "eS17": "curium",
                     "bS18": "deuterium",
                     "uS19": "dubnium",
                     "eS19": "dysprosium",
                     "bS20": "einsteinium",
                     "bS21": "erbium",
                     "bTHX": "europium",
                     "eS21": "fermium",
                     "eS24": "fluorine",
                     "eS25": "francium",
                     "eS26": "gadolinium",
                     "eS27": "gallium",
                     "eS28": "germanium",
                     "eS30": "gold",
                     "eS31": "hafnium",
                     }

def build_selection_string( chain_name: str, res_ids: List[int]):
    """pymol command being issued"""
    selection_name    = '23-28S_{}'.format(chain_name)
    selection_string  = f"c. {chain_name} and "
    selection_string += " OR ".join([*map(lambda x: f"resi {x}", res_ids)])
    return selection_name, selection_string

def list_bacteria():
    pprint(os.listdir(f"{RIBETL_DATA}/PTC_COORDINATES"))

def get_markerspath(struct: str):

    struct = struct.upper()
    _path  = os.path.join(RIBETL_DATA, "PTC_MARKERS_RAW", f"{struct}_PTC_MARKERS_RAW.json")

    return _path

def create_marker_at_atom(selection_name:str, posn:List[float], color_:str="red", repr="spheres", label=''):

    cmd.pseudoatom(selection_name, pos=posn, vdw=1, color=color_,  label=label)
    cmd.show(repr, selection_name)

def pseudoatom_ptc(rcsb_id: str):
    assets   = RibosomeAssets(rcsb_id)
    profile  = assets.profile()
    ptc_path = os.path.join(assets._dir_path(), "{}_PTC_COORDINATES.json".format(rcsb_id))


    if not os.path.isfile(ptc_path):
        reslist,auth_asym_id = ptc_resdiues_get(assets.biopython_structure(), profile.rnas, 0)
        midpoint = ptc_residues_calculate_midpoint(reslist,auth_asym_id)
    else:
        with open(ptc_path, 'r') as infile:
            ptc_coords = json.load(infile)
            midpoint = ptc_coords['midpoint_coordinates']


    create_marker_at_atom("centroid",midpoint, color_="red")

def by_chain(pdbid: str):
    pdbid          = pdbid.upper()
    RIBETL_DATA    = str(os.environ.get('RIBETL_DATA'))
    profilepath    = os.path.join(RIBETL_DATA, pdbid, f"{pdbid}.json")

    with open(profilepath, 'rb') as infile:
        profile = json.load(infile)

    if profile['rnas'] != None:
        for rna in profile['rnas']:

            if len( rna['nomenclature'] )>0:
                nomclass = rna['nomenclature'][0]
            else:
                nomclass = random.choice(list(colormap__RNA.keys()))
            prot_tmp = f"{nomclass}.{rna['auth_asym_id']}"
            cmd.create(prot_tmp, f"chain {rna[ 'auth_asym_id' ]}")
            cmd.remove(f"chain {rna[ 'auth_asym_id' ]} and m. {pdbid}")

            #! black surface with white lines inside
            # cmd.hide ( 'everything'           , f"chain {rna['auth_asym_id']}")
            # cmd.show ( 'surface'              , f"chain {rna['auth_asym_id']}")
            # cmd.show ( 'cartoon'              , f"chain {rna['auth_asym_id']}")
            # cmd.color( 'gray70', f"chain {rna['auth_asym_id']}")
            # cmd.set( 'transparency'  , 0.6    , f"chain {rna['auth_asym_id']}")
            # cmd.set( 'cartoon_color' , 'gray80', f"chain {rna['auth_asym_id']}")

            # #! transparent mesh with colored lines inside
            cmd.hide('everything', f"chain {rna['auth_asym_id']}")
            cmd.show('surface', f"chain {rna['auth_asym_id']}")
            cmd.show('cartoon', f"chain {rna['auth_asym_id']}")
            cmd.color(colormap__RNA[nomclass], f"chain {rna['auth_asym_id']}")
            cmd.set('transparency', 0.5, f"chain {rna['auth_asym_id']}")

    for protein in profile['proteins']:
        if len( protein['nomenclature'] )>0:
            nomclass = protein['nomenclature'][0]
        else:
            nomclass = random.choice([ *list(colormap__LSU_Proteins.keys()), *list(colormap__SSU_Proteins.keys()) ])

        if nomclass in colormap__LSU_Proteins:
            CLR = colormap__LSU_Proteins[nomclass]

        elif nomclass in colormap__SSU_Proteins:
            CLR = colormap__SSU_Proteins[nomclass]

        prot_tmp = f"{nomclass}.{protein['auth_asym_id']}"

        cmd.create(prot_tmp, f"chain {protein[ 'auth_asym_id' ]}")
        cmd.remove(f"chain {protein[ 'auth_asym_id' ]} and m. {pdbid}")

        #! black surface with white lines inside

        # cmd.hide ( 'everything'           , f"chain {protein['auth_asym_id']}")
        # cmd.show ( 'surface'              , f"chain {protein['auth_asym_id']}")
        # cmd.show ( 'cartoon'              , f"chain {protein['auth_asym_id']}")
        # cmd.color('gray70', f"chain {protein['auth_asym_id']}")
        # cmd.set( 'transparency'  , 0.5    , f"chain {protein['auth_asym_id']}")
        # cmd.set( 'cartoon_color' , 'black', f"chain {protein['auth_asym_id']}")
        #! transparent mesh with colored lines inside
        cmd.hide ('everything'           ,      f"chain {protein['auth_asym_id']}")
        cmd.show ('surface'                 ,      f"chain {protein['auth_asym_id']}")
        cmd.show ('cartoon'              ,      f"chain {protein['auth_asym_id']}")
        cmd.color(CLR                    ,      f"chain {protein['auth_asym_id']}")
        cmd.set  ('transparency' , 0.5, f"chain {protein['auth_asym_id']}")

    cmd.reset()

def sload(pdbid: str):
    pdbid       = pdbid.upper()
    RIBETL_DATA = str(os.environ.get('RIBETL_DATA'))
    path        = os.path.join(RIBETL_DATA, pdbid, f"{pdbid}.cif")
    cmd.delete('all')
    cmd.load(path)

def ray_picture(pdbid:str):
    cmd.reset()
    cmd.set('ray_trace_mode', 1)
    cmd.png(f"/home/rtviii/dev/riboxyz/rays/{pdbid}.png", ray=1, width=400, height=400, dpi=300)

def test_():
    print("Extended scripts loaded correctly")

def extract_chains(list_auth_asym_ids:str) :
    #! the argument is whitespace separated because pymol doesn't support lists apparently intuitively? 
    #! I.e. extract_chains LC L5

    chains_names   = list_auth_asym_ids.split(" ")
    selection_name = "Chains_{}".format(",".join(chains_names))

    cmd.extract(selection_name,  "c. {}".format("+".join(chains_names)) )
    return selection_name

def visualize_tunnel_mesh(rcsb_id: str):
    load_ply(poisson_recon_ascii_path(rcsb_id))

def centerline_expand(rcsb_id:str):
    print(centerline_expand)
    
def addAtom(model, name, vdw, x, y, z, partialCharge = 0.0):
  a = chempy.Atom()
  a.name = name
  a.vdw = vdw
  a.coord = [x, y, z]
  a.partial_charge = partialCharge
  model.atom.append(a)

def create_centerline(rcsb_id:str):

    rcsb_id = rcsb_id.upper()
    
    def open_tunnel_csv(rcsb_id: str) -> list[list]:

        TUNNEL_PATH = os.path.join( EXIT_TUNNEL_WORK, "mole_tunnels", "tunnel_{}.csv".format(rcsb_id) )
        df          = pd.read_csv(TUNNEL_PATH)
        data        = []

        for index, row in df.iterrows():
            radius       = row["Radius"]
            x_coordinate = row["X"]
            y_coordinate = row["Y"]
            z_coordinate = row["Z"]
            data.append([radius, x_coordinate, y_coordinate, z_coordinate])
        return data

    import chempy
    object_name = "centerline_{}".format(rcsb_id)

    data = open_tunnel_csv(rcsb_id)
    model = chempy.models.Indexed() 
    for i, ( r,x,y,z ) in enumerate(data):
          addAtom(model,str(i),r,x,y,z)
    for a in range(len(model.atom)-1):
        b = chempy.Bond()
        b.index = [a,a+1]
        model.bond.append(b)
    cmd.set('surface_mode', 1)
    cmd.set('sphere_mode', 5)
    cmd.set('mesh_mode', 1)
    cmd.load_model(model, object_name)
    cmd.hide('everything', object_name)
    cmd.set('sphere_color', 'red', object_name)
    cmd.show('spheres', object_name)

    cmd.select("_surf25", "select {} within 25 of {}".format(rcsb_id,object_name))
    cmd.create("surf25", "_surf25")
    cmd.select("surf15", "select {} within 15 of {}".format(rcsb_id,object_name))
    cmd.create("surf15", "_surf15")
    cmd.deselect()





cmd.extend("extract_chains"    , extract_chains        )
cmd.extend("sload"             , sload                 )
cmd.extend("by_chain"          , by_chain              )
cmd.extend("ptc"               , pseudoatom_ptc        )
cmd.extend("list_bacteria"     , list_bacteria         )
cmd.extend("test_"             , test_                 )
cmd.extend("ray_picture"       , ray_picture           )
cmd.extend("tunnel_mesh"       , visualize_tunnel_mesh )
cmd.extend("create_centerline" , create_centerline     )


# cmd.extend("tun_obstructions" , visualize_obstructions )
# cmd.extend("exitport"         , pseudoatom_exitport    )
    

