import json
import os
from pprint import pprint
from typing import List
from pymol import cmd
from api.scratch_tunnel_workflow import exit_port_posn, ptc_residues_calculate_midpoint, ptc_resdiues_get, tunnel_obstructions

RIBETL_DATA = os.environ.get("RIBETL_DATA")

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

def pseudoatom_ptc(struct: str):
    reslist,auth_asym_id = ptc_resdiues_get(struct, 0)
    midpoint = ptc_residues_calculate_midpoint(reslist,auth_asym_id)
    create_marker_at_atom("centroid",midpoint, color_="red")

def visualize_obstructions(rcsb_id):
    ptcres, auth_asym_id = ptc_resdiues_get(rcsb_id, 0)
    midpoint = ptc_residues_calculate_midpoint(ptcres, auth_asym_id)
    polys,nonpolys = tunnel_obstructions(rcsb_id, midpoint)

    cmd.color("white", "all")

    pseudoatom_ptc(rcsb_id)

    for poly in polys:
        cname = "chain_{}".format(poly.auth_asym_id, )
        csele = "c. {} and m. {}".format(poly.auth_asym_id, rcsb_id.upper()) 
        cmd.create(cname, csele)
        cmd.remove(csele)
        cmd.color("red", cname)

    for res in nonpolys:
        resname = "res_{}_{}".format(res.parent_auth_asym_id,res.resname)
        ressele = "m. {} and c. {} and resn {}".format(rcsb_id.upper(),res.parent_auth_asym_id,res.resname)
        cmd.create(resname,  ressele)
        cmd.remove(ressele)
        cmd.color("blue", resname)

def struct_paint_chains(pdbid: str):
    pdbid          = pdbid.upper()
    RIBETL_DATA    = str(os.environ.get('RIBETL_DATA'))
    nomenclaturev2 = os.path.join("/home/rxz/dev/docker_ribxz/api/ribctl/assets/nomenclaturev2", f"{pdbid.upper()}.json")
    profilepath    = os.path.join(RIBETL_DATA, pdbid, f"{pdbid}.json")

    with open(nomenclaturev2, 'rb') as ninfile:
        nomenclaturev2 = json.load(ninfile)

    with open(profilepath, 'rb') as infile:
        profile = json.load(infile)

    if profile['rnas'] != None:
        for rna in profile['rnas']:
            # cmd.color('white', f"chain {rna['auth_asym_id']}")
            cmd.hide('everything', f"chain {rna['auth_asym_id']}")
            # cmd.show('surface', f"chain {rna['auth_asym_id']}")
            cmd.show('cartoon', f"chain {rna['auth_asym_id']}")
            cmd.color('white', f"chain {rna['auth_asym_id']}")
            # cmd.set('transparency', 0.3, f"chain {rna['auth_asym_id']}")
            cmd.set('cartoon_transparency', 0.2, f"chain {rna['auth_asym_id']}")

    for protein in profile['proteins']:
        nomclass = nomenclaturev2[protein['auth_asym_id']]

        if nomclass in colormap__LSU_Proteins:
            CLR = colormap__LSU_Proteins[nomclass]

        elif nomclass in colormap__SSU_Proteins:
            CLR = colormap__SSU_Proteins[nomclass]

        prot_tmp = f"{nomclass}.{protein['auth_asym_id']}"

        cmd.create(prot_tmp, f"chain {protein[ 'auth_asym_id' ]}")

        cmd.remove(f"chain {protein[ 'auth_asym_id' ]} and m. {pdbid}")

        cmd.hide('everything', prot_tmp)
        cmd.show('surface', prot_tmp)
        cmd.show('sticks', prot_tmp)
        cmd.show('cartoon', prot_tmp)
        cmd.show('lines', prot_tmp)
        cmd.set('transparency', 0.75, prot_tmp)
        cmd.color(CLR , prot_tmp)

def sload(pdbid: str):
    pdbid       = pdbid.upper()
    RIBETL_DATA = str(os.environ.get('RIBETL_DATA'))
    path        = os.path.join(RIBETL_DATA, pdbid, f"{pdbid}.cif")
    cmd.delete('all')
    cmd.load(path)

def pseudoatom_exitport(rcsb_id:str):
    posn = exit_port_posn(rcsb_id)
    print("Got positions :", posn)
    create_marker_at_atom("Exitport", posn)


cmd.extend("sload"            , sload                  )
cmd.extend("by_chain"         , struct_paint_chains    )
cmd.extend("ptc"              , pseudoatom_ptc         )
cmd.extend("list_bacteria"    , list_bacteria          )
cmd.extend("tun_obstructions" , visualize_obstructions )
cmd.extend("exitport"         , pseudoatom_exitport    )
    

