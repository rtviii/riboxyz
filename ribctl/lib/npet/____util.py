import numpy as np



#TODO: replace with request.get('https://api.ribosome.xyz/{rcsb_id}/landmark/{ptc,constriction_site}')
# Leaving static data here for now for the structs in the paper. 
# This data will evenutally be served from the dedicated endpoint 


# def landmark_ptc(rcsb_id: str) -> np.ndarray:
#     _ = {
#         "4UG0": {
#             "site_9_residues": [
#                 ["G", " "],
#                 ["A", " "],
#                 ["G", " "],
#                 ["C", " "],
#                 ["U", " "],
#                 ["G", " "],
#                 ["G", " "],
#                 ["G", " "],
#                 ["U", " "],
#                 ["U", " "],
#                 ["U", " "],
#                 ["A", " "],
#             ],
#             "LSU_rRNA_auth_asym_id": "L5",
#             "midpoint_coordinates": [
#                 139.55250549316406,
#                 155.26100158691406,
#                 153.34950256347656,
#             ],
#             "nomenclature_table": {
#                 "S6": {
#                     "nomenclature": ["tRNA"],
#                     "entity_poly_strand_id": "S6",
#                     "rcsb_pdbx_description": "HUMAN INITIATOR MET-TRNA-I",
#                 },
#                 "LB": {
#                     "nomenclature": ["uL3"],
#                     "entity_poly_strand_id": "LB",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L3",
#                 },
#                 "LD": {
#                     "nomenclature": ["uL18"],
#                     "entity_poly_strand_id": "LD",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L5",
#                 },
#                 "LU": {
#                     "nomenclature": ["eL22"],
#                     "entity_poly_strand_id": "LU",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L22",
#                 },
#                 "LW": {
#                     "nomenclature": ["eL24"],
#                     "entity_poly_strand_id": "LW",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L24",
#                 },
#                 "LZ": {
#                     "nomenclature": ["eL27"],
#                     "entity_poly_strand_id": "LZ",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L27",
#                 },
#                 "Lc": {
#                     "nomenclature": ["eL30"],
#                     "entity_poly_strand_id": "Lc",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L30",
#                 },
#                 "Le": {
#                     "nomenclature": ["eL32"],
#                     "entity_poly_strand_id": "Le",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L32",
#                 },
#                 "Lf": {
#                     "nomenclature": ["eL33"],
#                     "entity_poly_strand_id": "Lf",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L35A",
#                 },
#                 "Lk": {
#                     "nomenclature": ["eL38"],
#                     "entity_poly_strand_id": "Lk",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L38",
#                 },
#                 "Lr": {
#                     "nomenclature": ["eL28"],
#                     "entity_poly_strand_id": "Lr",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L28",
#                 },
#                 "SE": {
#                     "nomenclature": ["eS4"],
#                     "entity_poly_strand_id": "SE",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S4, X ISOFORM",
#                 },
#                 "SI": {
#                     "nomenclature": ["eS8"],
#                     "entity_poly_strand_id": "SI",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S8",
#                 },
#                 "SP": {
#                     "nomenclature": ["uS19"],
#                     "entity_poly_strand_id": "SP",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S15",
#                 },
#                 "SR": {
#                     "nomenclature": ["eS17"],
#                     "entity_poly_strand_id": "SR",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S17-LIKE",
#                 },
#                 "SS": {
#                     "nomenclature": ["uS13"],
#                     "entity_poly_strand_id": "SS",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S18",
#                 },
#                 "SV": {
#                     "nomenclature": ["eS21"],
#                     "entity_poly_strand_id": "SV",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S21",
#                 },
#                 "SX": {
#                     "nomenclature": ["uS12"],
#                     "entity_poly_strand_id": "SX",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S23",
#                 },
#                 "Sc": {
#                     "nomenclature": ["eS28"],
#                     "entity_poly_strand_id": "Sc",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S28",
#                 },
#                 "SC": {
#                     "nomenclature": ["uS5"],
#                     "entity_poly_strand_id": "SC",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S2",
#                 },
#                 "LG": {
#                     "nomenclature": ["eL8"],
#                     "entity_poly_strand_id": "LG",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L7A",
#                 },
#                 "LH": {
#                     "nomenclature": ["uL6"],
#                     "entity_poly_strand_id": "LH",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L9",
#                 },
#                 "LI": {
#                     "nomenclature": ["uL16"],
#                     "entity_poly_strand_id": "LI",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L10-LIKE",
#                 },
#                 "LJ": {
#                     "nomenclature": ["uL5"],
#                     "entity_poly_strand_id": "LJ",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L11",
#                 },
#                 "LL": {
#                     "nomenclature": ["eL13"],
#                     "entity_poly_strand_id": "LL",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L13",
#                 },
#                 "LM": {
#                     "nomenclature": ["eL14"],
#                     "entity_poly_strand_id": "LM",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L14",
#                 },
#                 "LN": {
#                     "nomenclature": ["eL15"],
#                     "entity_poly_strand_id": "LN",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L15",
#                 },
#                 "LO": {
#                     "nomenclature": ["uL13"],
#                     "entity_poly_strand_id": "LO",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L13A",
#                 },
#                 "LP": {
#                     "nomenclature": ["uL22"],
#                     "entity_poly_strand_id": "LP",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L17",
#                 },
#                 "LQ": {
#                     "nomenclature": ["eL18"],
#                     "entity_poly_strand_id": "LQ",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L18",
#                 },
#                 "LR": {
#                     "nomenclature": ["eL19"],
#                     "entity_poly_strand_id": "LR",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L19",
#                 },
#                 "LS": {
#                     "nomenclature": ["eL20"],
#                     "entity_poly_strand_id": "LS",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L18A",
#                 },
#                 "LT": {
#                     "nomenclature": ["eL21"],
#                     "entity_poly_strand_id": "LT",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L21",
#                 },
#                 "LV": {
#                     "nomenclature": ["uL14"],
#                     "entity_poly_strand_id": "LV",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L23",
#                 },
#                 "LX": {
#                     "nomenclature": ["uL23"],
#                     "entity_poly_strand_id": "LX",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L23A",
#                 },
#                 "LY": {
#                     "nomenclature": ["uL24"],
#                     "entity_poly_strand_id": "LY",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L26",
#                 },
#                 "La": {
#                     "nomenclature": ["uL15"],
#                     "entity_poly_strand_id": "La",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L27A",
#                 },
#                 "Lb": {
#                     "nomenclature": ["eL29"],
#                     "entity_poly_strand_id": "Lb",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L29",
#                 },
#                 "Ld": {
#                     "nomenclature": ["eL31"],
#                     "entity_poly_strand_id": "Ld",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L31",
#                 },
#                 "Lg": {
#                     "nomenclature": ["eL34"],
#                     "entity_poly_strand_id": "Lg",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L34",
#                 },
#                 "Lh": {
#                     "nomenclature": ["uL29"],
#                     "entity_poly_strand_id": "Lh",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L35",
#                 },
#                 "Li": {
#                     "nomenclature": ["eL36"],
#                     "entity_poly_strand_id": "Li",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L36",
#                 },
#                 "Lj": {
#                     "nomenclature": ["eL37"],
#                     "entity_poly_strand_id": "Lj",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L37",
#                 },
#                 "LA": {
#                     "nomenclature": ["uL2"],
#                     "entity_poly_strand_id": "LA",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L8",
#                 },
#                 "Ll": {
#                     "nomenclature": ["eL39"],
#                     "entity_poly_strand_id": "Ll",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L39",
#                 },
#                 "Lm": {
#                     "nomenclature": ["eL40"],
#                     "entity_poly_strand_id": "Lm",
#                     "rcsb_pdbx_description": "UBIQUITIN-60S RIBOSOMAL PROTEIN L40",
#                 },
#                 "Ln": {
#                     "nomenclature": ["eL41"],
#                     "entity_poly_strand_id": "Ln",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L41",
#                 },
#                 "Lo": {
#                     "nomenclature": ["eL42"],
#                     "entity_poly_strand_id": "Lo",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L36A",
#                 },
#                 "Lp": {
#                     "nomenclature": ["eL43"],
#                     "entity_poly_strand_id": "Lp",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L37A",
#                 },
#                 "Lz": {
#                     "nomenclature": ["uL1"],
#                     "entity_poly_strand_id": "Lz",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L10A",
#                 },
#                 "SA": {
#                     "nomenclature": ["uS2"],
#                     "entity_poly_strand_id": "SA",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN SA",
#                 },
#                 "SB": {
#                     "nomenclature": ["eS1"],
#                     "entity_poly_strand_id": "SB",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S3A",
#                 },
#                 "SD": {
#                     "nomenclature": ["uS3"],
#                     "entity_poly_strand_id": "SD",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S3",
#                 },
#                 "SF": {
#                     "nomenclature": ["uS7"],
#                     "entity_poly_strand_id": "SF",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S5",
#                 },
#                 "SH": {
#                     "nomenclature": ["eS7"],
#                     "entity_poly_strand_id": "SH",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S7",
#                 },
#                 "SK": {
#                     "nomenclature": ["eS10"],
#                     "entity_poly_strand_id": "SK",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S10",
#                 },
#                 "SL": {
#                     "nomenclature": ["uS17"],
#                     "entity_poly_strand_id": "SL",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S11",
#                 },
#                 "SQ": {
#                     "nomenclature": ["uS9"],
#                     "entity_poly_strand_id": "SQ",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S16",
#                 },
#                 "LC": {
#                     "nomenclature": ["uL4"],
#                     "entity_poly_strand_id": "LC",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L4",
#                 },
#                 "ST": {
#                     "nomenclature": ["eS19"],
#                     "entity_poly_strand_id": "ST",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S19",
#                 },
#                 "SU": {
#                     "nomenclature": ["uS10"],
#                     "entity_poly_strand_id": "SU",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S20",
#                 },
#                 "Sa": {
#                     "nomenclature": ["eS26"],
#                     "entity_poly_strand_id": "Sa",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S26",
#                 },
#                 "Sd": {
#                     "nomenclature": ["uS14"],
#                     "entity_poly_strand_id": "Sd",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S29",
#                 },
#                 "Sf": {
#                     "nomenclature": ["eS31"],
#                     "entity_poly_strand_id": "Sf",
#                     "rcsb_pdbx_description": "UBIQUITIN-40S RIBOSOMAL PROTEIN S27A",
#                 },
#                 "Sg": {
#                     "nomenclature": ["RACK1"],
#                     "entity_poly_strand_id": "Sg",
#                     "rcsb_pdbx_description": "GUANINE NUCLEOTIDE-BINDING PROTEIN SUBUNIT BETA-2-LIKE 1",
#                 },
#                 "SG": {
#                     "nomenclature": ["eS6"],
#                     "entity_poly_strand_id": "SG",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S6",
#                 },
#                 "SJ": {
#                     "nomenclature": ["uS4"],
#                     "entity_poly_strand_id": "SJ",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S9",
#                 },
#                 "SM": {
#                     "nomenclature": ["eS12"],
#                     "entity_poly_strand_id": "SM",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN",
#                 },
#                 "SN": {
#                     "nomenclature": ["uS15"],
#                     "entity_poly_strand_id": "SN",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S13",
#                 },
#                 "SO": {
#                     "nomenclature": ["uS11"],
#                     "entity_poly_strand_id": "SO",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S14",
#                 },
#                 "SW": {
#                     "nomenclature": ["uS8"],
#                     "entity_poly_strand_id": "SW",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S15A",
#                 },
#                 "SY": {
#                     "nomenclature": ["eS24"],
#                     "entity_poly_strand_id": "SY",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S24",
#                 },
#                 "SZ": {
#                     "nomenclature": ["eS25"],
#                     "entity_poly_strand_id": "SZ",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S25",
#                 },
#                 "LE": {
#                     "nomenclature": ["eL6"],
#                     "entity_poly_strand_id": "LE",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L6",
#                 },
#                 "Sb": {
#                     "nomenclature": ["eS27"],
#                     "entity_poly_strand_id": "Sb",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S27",
#                 },
#                 "Se": {
#                     "nomenclature": ["eS30"],
#                     "entity_poly_strand_id": "Se",
#                     "rcsb_pdbx_description": "40S RIBOSOMAL PROTEIN S30",
#                 },
#                 "LF": {
#                     "nomenclature": ["uL30"],
#                     "entity_poly_strand_id": "LF",
#                     "rcsb_pdbx_description": "60S RIBOSOMAL PROTEIN L7",
#                 },
#                 "L5": {
#                     "nomenclature": ["28SrRNA"],
#                     "entity_poly_strand_id": "L5",
#                     "rcsb_pdbx_description": "28S ribosomal RNA",
#                 },
#                 "L7": {
#                     "nomenclature": ["5SrRNA"],
#                     "entity_poly_strand_id": "L7",
#                     "rcsb_pdbx_description": "5S ribosomal RNA",
#                 },
#                 "L8": {
#                     "nomenclature": ["5.8SrRNA"],
#                     "entity_poly_strand_id": "L8",
#                     "rcsb_pdbx_description": "5.8S ribosomal RNA",
#                 },
#                 "S2": {
#                     "nomenclature": [],
#                     "entity_poly_strand_id": "S2",
#                     "rcsb_pdbx_description": "18S ribosomal RNA",
#                 },
#             },
#         },
#         "6WD4": {
#             "site_9_residues": [
#                 ["G", " ", ["6WD4", 0, "1", [" ", 2576, " "]]],
#                 ["A", " ", ["6WD4", 0, "1", [" ", 2577, " "]]],
#                 ["G", " ", ["6WD4", 0, "1", [" ", 2578, " "]]],
#                 ["C", " ", ["6WD4", 0, "1", [" ", 2579, " "]]],
#                 ["U", " ", ["6WD4", 0, "1", [" ", 2580, " "]]],
#                 ["G", " ", ["6WD4", 0, "1", [" ", 2581, " "]]],
#                 ["G", " ", ["6WD4", 0, "1", [" ", 2582, " "]]],
#                 ["G", " ", ["6WD4", 0, "1", [" ", 2583, " "]]],
#                 ["U", " ", ["6WD4", 0, "1", [" ", 2584, " "]]],
#                 ["U", " ", ["6WD4", 0, "1", [" ", 2585, " "]]],
#                 ["U", " ", ["6WD4", 0, "1", [" ", 2586, " "]]],
#                 ["A", " ", ["6WD4", 0, "1", [" ", 2587, " "]]],
#             ],
#             "LSU_rRNA_auth_asym_id": "1",
#             "midpoint_coordinates": [
#                 213.7469940185547,
#                 187.90499877929688,
#                 177.95999908447266,
#             ],
#             "nomenclature_table": {
#                 "b": {
#                     "nomenclature": ["uL2"],
#                     "entity_poly_strand_id": "b",
#                     "rcsb_pdbx_description": "50S ribosomal protein L2",
#                 },
#                 "k": {
#                     "nomenclature": ["uL14"],
#                     "entity_poly_strand_id": "k",
#                     "rcsb_pdbx_description": "50S ribosomal protein L14",
#                 },
#                 "l": {
#                     "nomenclature": ["uL15"],
#                     "entity_poly_strand_id": "l",
#                     "rcsb_pdbx_description": "50S ribosomal protein L15",
#                 },
#                 "m": {
#                     "nomenclature": ["uL16"],
#                     "entity_poly_strand_id": "m",
#                     "rcsb_pdbx_description": "50S ribosomal protein L16",
#                 },
#                 "n": {
#                     "nomenclature": ["bL17"],
#                     "entity_poly_strand_id": "n",
#                     "rcsb_pdbx_description": "50S ribosomal protein L17",
#                 },
#                 "o": {
#                     "nomenclature": ["uL18"],
#                     "entity_poly_strand_id": "o",
#                     "rcsb_pdbx_description": "50S ribosomal protein L18",
#                 },
#                 "p": {
#                     "nomenclature": ["bL19"],
#                     "entity_poly_strand_id": "p",
#                     "rcsb_pdbx_description": "50S ribosomal protein L19",
#                 },
#                 "q": {
#                     "nomenclature": ["bL20"],
#                     "entity_poly_strand_id": "q",
#                     "rcsb_pdbx_description": "50S ribosomal protein L20",
#                 },
#                 "r": {
#                     "nomenclature": ["bL21"],
#                     "entity_poly_strand_id": "r",
#                     "rcsb_pdbx_description": "50S ribosomal protein L21",
#                 },
#                 "s": {
#                     "nomenclature": ["uL22"],
#                     "entity_poly_strand_id": "s",
#                     "rcsb_pdbx_description": "50S ribosomal protein L22",
#                 },
#                 "t": {
#                     "nomenclature": ["uL23"],
#                     "entity_poly_strand_id": "t",
#                     "rcsb_pdbx_description": "50S ribosomal protein L23",
#                 },
#                 "c": {
#                     "nomenclature": ["uL3"],
#                     "entity_poly_strand_id": "c",
#                     "rcsb_pdbx_description": "50S ribosomal protein L3",
#                 },
#                 "u": {
#                     "nomenclature": ["uL24"],
#                     "entity_poly_strand_id": "u",
#                     "rcsb_pdbx_description": "50S ribosomal protein L24",
#                 },
#                 "v": {
#                     "nomenclature": ["bL25"],
#                     "entity_poly_strand_id": "v",
#                     "rcsb_pdbx_description": "50S ribosomal protein L25",
#                 },
#                 "w": {
#                     "nomenclature": ["bL27"],
#                     "entity_poly_strand_id": "w",
#                     "rcsb_pdbx_description": "50S ribosomal protein L27",
#                 },
#                 "x": {
#                     "nomenclature": ["bL28"],
#                     "entity_poly_strand_id": "x",
#                     "rcsb_pdbx_description": "50S ribosomal protein L28",
#                 },
#                 "y": {
#                     "nomenclature": ["uL29"],
#                     "entity_poly_strand_id": "y",
#                     "rcsb_pdbx_description": "50S ribosomal protein L29",
#                 },
#                 "z": {
#                     "nomenclature": ["uL30"],
#                     "entity_poly_strand_id": "z",
#                     "rcsb_pdbx_description": "50S ribosomal protein L30",
#                 },
#                 "B": {
#                     "nomenclature": ["bL32"],
#                     "entity_poly_strand_id": "B",
#                     "rcsb_pdbx_description": "50S ribosomal protein L32",
#                 },
#                 "C": {
#                     "nomenclature": ["bL33"],
#                     "entity_poly_strand_id": "C",
#                     "rcsb_pdbx_description": "50S ribosomal protein L33",
#                 },
#                 "D": {
#                     "nomenclature": ["bL34"],
#                     "entity_poly_strand_id": "D",
#                     "rcsb_pdbx_description": "50S ribosomal protein L34",
#                 },
#                 "E": {
#                     "nomenclature": ["bL35"],
#                     "entity_poly_strand_id": "E",
#                     "rcsb_pdbx_description": "50S ribosomal protein L35",
#                 },
#                 "d": {
#                     "nomenclature": ["uL4"],
#                     "entity_poly_strand_id": "d",
#                     "rcsb_pdbx_description": "50S ribosomal protein L4",
#                 },
#                 "F": {
#                     "nomenclature": ["bL36"],
#                     "entity_poly_strand_id": "F",
#                     "rcsb_pdbx_description": "50S ribosomal protein L36",
#                 },
#                 "G": {
#                     "nomenclature": ["uS2"],
#                     "entity_poly_strand_id": "G",
#                     "rcsb_pdbx_description": "30S ribosomal protein S2",
#                 },
#                 "H": {
#                     "nomenclature": ["uS3"],
#                     "entity_poly_strand_id": "H",
#                     "rcsb_pdbx_description": "30S ribosomal protein S3",
#                 },
#                 "I": {
#                     "nomenclature": ["uS4"],
#                     "entity_poly_strand_id": "I",
#                     "rcsb_pdbx_description": "30S ribosomal protein S4",
#                 },
#                 "J": {
#                     "nomenclature": ["uS5"],
#                     "entity_poly_strand_id": "J",
#                     "rcsb_pdbx_description": "30S ribosomal protein S5",
#                 },
#                 "K": {
#                     "nomenclature": ["bS6"],
#                     "entity_poly_strand_id": "K",
#                     "rcsb_pdbx_description": "30S ribosomal protein S6",
#                 },
#                 "L": {
#                     "nomenclature": ["uS7"],
#                     "entity_poly_strand_id": "L",
#                     "rcsb_pdbx_description": "30S ribosomal protein S7",
#                 },
#                 "M": {
#                     "nomenclature": ["uS8"],
#                     "entity_poly_strand_id": "M",
#                     "rcsb_pdbx_description": "30S ribosomal protein S8",
#                 },
#                 "N": {
#                     "nomenclature": ["uS9"],
#                     "entity_poly_strand_id": "N",
#                     "rcsb_pdbx_description": "30S ribosomal protein S9",
#                 },
#                 "O": {
#                     "nomenclature": ["uS10"],
#                     "entity_poly_strand_id": "O",
#                     "rcsb_pdbx_description": "30S ribosomal protein S10",
#                 },
#                 "e": {
#                     "nomenclature": ["uL5"],
#                     "entity_poly_strand_id": "e",
#                     "rcsb_pdbx_description": "50S ribosomal protein L5",
#                 },
#                 "P": {
#                     "nomenclature": ["uS11"],
#                     "entity_poly_strand_id": "P",
#                     "rcsb_pdbx_description": "30S ribosomal protein S11",
#                 },
#                 "Q": {
#                     "nomenclature": ["uS12"],
#                     "entity_poly_strand_id": "Q",
#                     "rcsb_pdbx_description": "30S ribosomal protein S12",
#                 },
#                 "R": {
#                     "nomenclature": ["uS13"],
#                     "entity_poly_strand_id": "R",
#                     "rcsb_pdbx_description": "30S ribosomal protein S13",
#                 },
#                 "S": {
#                     "nomenclature": ["uS14"],
#                     "entity_poly_strand_id": "S",
#                     "rcsb_pdbx_description": "30S ribosomal protein S14",
#                 },
#                 "T": {
#                     "nomenclature": ["uS15"],
#                     "entity_poly_strand_id": "T",
#                     "rcsb_pdbx_description": "30S ribosomal protein S15",
#                 },
#                 "U": {
#                     "nomenclature": ["bS16"],
#                     "entity_poly_strand_id": "U",
#                     "rcsb_pdbx_description": "30S ribosomal protein S16",
#                 },
#                 "V": {
#                     "nomenclature": ["uS17"],
#                     "entity_poly_strand_id": "V",
#                     "rcsb_pdbx_description": "30S ribosomal protein S17",
#                 },
#                 "W": {
#                     "nomenclature": ["bS18"],
#                     "entity_poly_strand_id": "W",
#                     "rcsb_pdbx_description": "30S ribosomal protein S18",
#                 },
#                 "X": {
#                     "nomenclature": ["uS19"],
#                     "entity_poly_strand_id": "X",
#                     "rcsb_pdbx_description": "30S ribosomal protein S19",
#                 },
#                 "Y": {
#                     "nomenclature": ["bS20"],
#                     "entity_poly_strand_id": "Y",
#                     "rcsb_pdbx_description": "30S ribosomal protein S20",
#                 },
#                 "f": {
#                     "nomenclature": ["uL6"],
#                     "entity_poly_strand_id": "f",
#                     "rcsb_pdbx_description": "50S ribosomal protein L6",
#                 },
#                 "Z": {
#                     "nomenclature": ["bS21"],
#                     "entity_poly_strand_id": "Z",
#                     "rcsb_pdbx_description": "30S ribosomal protein S21",
#                 },
#                 "a": {
#                     "nomenclature": ["uL1"],
#                     "entity_poly_strand_id": "a",
#                     "rcsb_pdbx_description": "50S ribosomal protein L1",
#                 },
#                 "8": {
#                     "nomenclature": ["EF_Tu"],
#                     "entity_poly_strand_id": "8",
#                     "rcsb_pdbx_description": "Elongation factor Tu",
#                 },
#                 "g": {
#                     "nomenclature": ["bL9"],
#                     "entity_poly_strand_id": "g",
#                     "rcsb_pdbx_description": "50S ribosomal protein L9",
#                 },
#                 "h": {
#                     "nomenclature": ["uL10"],
#                     "entity_poly_strand_id": "h",
#                     "rcsb_pdbx_description": "50S ribosomal protein L10",
#                 },
#                 "i": {
#                     "nomenclature": ["uL11"],
#                     "entity_poly_strand_id": "i",
#                     "rcsb_pdbx_description": "50S ribosomal protein L11",
#                 },
#                 "j": {
#                     "nomenclature": ["uL13"],
#                     "entity_poly_strand_id": "j",
#                     "rcsb_pdbx_description": "50S ribosomal protein L13",
#                 },
#                 "3": {
#                     "nomenclature": ["rRNA_16S"],
#                     "entity_poly_strand_id": "3",
#                     "rcsb_pdbx_description": "16S ribosomal RNA",
#                 },
#                 "1": {
#                     "nomenclature": ["rRNA_23S"],
#                     "entity_poly_strand_id": "1",
#                     "rcsb_pdbx_description": "23S ribosomal RNA",
#                 },
#                 "2": {
#                     "nomenclature": ["rRNA_5S"],
#                     "entity_poly_strand_id": "2",
#                     "rcsb_pdbx_description": "5S ribosomal RNA",
#                 },
#                 "5": {
#                     "nomenclature": ["tRNA"],
#                     "entity_poly_strand_id": "5,6",
#                     "rcsb_pdbx_description": "tRNAfMet",
#                 },
#                 "6": {
#                     "nomenclature": ["tRNA"],
#                     "entity_poly_strand_id": "5,6",
#                     "rcsb_pdbx_description": "tRNAfMet",
#                 },
#                 "4": {
#                     "nomenclature": [],
#                     "entity_poly_strand_id": "4",
#                     "rcsb_pdbx_description": "mRNA",
#                 },
#                 "7": {
#                     "nomenclature": [],
#                     "entity_poly_strand_id": "7",
#                     "rcsb_pdbx_description": "tRNAPhe",
#                 },
#             },
#         },
#     }
#     return np.array(_[rcsb_id.upper()]["midpoint_coordinates"])


# def landmark_constriction_site(rcsb_id: str) -> np.ndarray:
#     _ = {
#         "4UG0": {
#             "location": [98.67242431640625, 160.83001708984375, 155.1837615966797]
#         },
#         "6WD4": {
#             "location": [234.11978149414062, 178.21743774414062, 185.29786682128906]
#         },
#     }
#     return np.array(_[rcsb_id.upper()]["location"])
