import sys
sys.path.append("/home/rtviii/dev/riboxyz")
from pprint import pprint
from rdkit import Chem
from rdkit.Chem import AllChem
from neo4j_ribosome.db_lib_reader import Neo4jReader
from ribctl.lib.schema.types_ribosome import NonpolymericLigand

profiled = ['SCM',
 'ON0',
 '3TS',
 'AKN',
 '3KD',
 '3K5',
 'A1AE1',
 'SPR',
 '773',
 'PG4',
 'PEG',
 'EDO',
 'PGE',
 '1PE',
 'ACY',
 'GUN',
 'TRS',
 'YXM',
 '6EM',
 '95H',
 'OIY',
 'FSD',
 'WUX',
 'ZIT',
 'NEG',
 'EZG',
 'G6V',
 'SY5',
 'A1AE0',
 'O',
 '5CR',
 '3J2',
 '5GP',
 'A',
 'C',
 'YAT',
 'WC9',
 'EDS',
 'SFG',
 'VIR',
 '3QB',
 'AMP',
 '3LK',
 'CAI',
 '84G',
 'HY0',
 '5I0',
 'MRD',
 'T8B',
 'RD8',
 'W9C',
 'TIX',
 'PHA',
 'MAU',
 'BME',
 'SLD',
 'MUL',
 'AQJ',
 'SJE',
 'CPT',
 'A1H4F',
 'ACA',
 'UAM',
 'U7V',
 '8Q1',
 'WDP',
 'U',
 'OMG',
 'G',
 'PSU',
 'WIN',
 'MT9',
 'A3P',
 '2AE',
 '3V6',
 '6UQ',
 '7MB',
 '7AL',
 '13T',
 'G34',
 '6NO',
 '1F2',
 '1F3',
 'V7A',
 'SPK',
 'K16',
 'PPU',
 '3J6',
 'EVN',
 'N',
 '3K8',
 '4M2',
 'GOL',
 '80P',
 '917',
 '6O1',
 'EUS',
 'ERY',
 'GCP',
 'MVM',
 'SPD',
 'SPM',
 'OHX',
 'HMT',
 'DI0',
 'PAR',
 'PUT',
 'ATP',
 'GTP',
 'SAH',
 'ILE',
 'SRY',
 'NAD',
 'GDP',
 'VAL',
 'UNK',
 'GNP',
 'ADP',
 'HGR',
 'SO1',
 '34G',
 'PHE',
 'GSP',
 'M5Z',
 'CLM',
 'ALA',
 'ANM',
 'FUA',
 'NMY',
 'DOL',
 'HYG',
 'IAS',
 'PM8',
 'LC2',
 'LMA',
 'PNS',
 'SPS',
 'HKO',
 'MYL',
 '7C4',
 'MQ6',
 'LMT',
 'CLY',
 'YMZ',
 'TYK',
 'KSG',
 'FS2',
 'TRP',
 'T1C',
 'BLS',
 'YQM',
 'U6A',
 'MPD',
 'ARG',
 'B3P',
 'KIR',
 'ZLD',
 'YRB',
 'SIS',
 'SER',
 'Z2V',
 'AM2',
 'B6M',
 'PRO',
 'GAL',
 'EDE',
 'ANP',
 'LUJ',
 'PUY',
 'TAC',
 'A1AEZ',
 'SEC',
 'G4P',
 'ASP',
 'OI9',
 'AGS',
 'PCY',
 '84D',
 'GET',
 'TEL',
 'ZC0',
 'LLL',
 'MAN',
 '8UZ',
 'UTP',
 'FAD',
 'VIF',
 'YRW',
 '8AN',
 'LYS',
 'A1D6G',
 '62B',
 'EOH',
 'EPE',
 '1F4']

def collect_all_ligands_with_smiles():
    reader      = Neo4jReader()
    ligs        = reader.list_ligands(nodes_only=True)[0][0]
    ligands     = list(map(NonpolymericLigand.model_validate, ligs))
    ligands = list(filter(lambda x:  x.chemicalId not in profiled, ligands))
    pprint(ligands)
    pprint(len(ligands))

    query_input = ""
    lig:NonpolymericLigand
    for lig in ligands:
        query_input = query_input + f"{lig.chemicalId}\t{lig.SMILES_stereo}\n"
        # query_input = query_input + f"{lig.chemicalId}\t {lig.chemicalName}\n"
    return query_input


def init_classification_class():
  import requests
  url = "http://classyfire.wishartlab.com/queries"

  headers = {
      "Accept": "application/json",
      "Content-Type": "application/json"
  }

  data = {
      "label"      : "All ligands | SMiles stereo",
      # "query_input": "MOL1\\tCCCOCC\\nMOL2\\tCOCC=CCCC",
      "query_input": collect_all_ligands_with_smiles(),
      "query_type" : "STRUCTURE"
  }

  response = requests.post(url, json=data, headers=headers)

  print(response.status_code)
  print(response.json())

init_classification_class()

# 2/3 : {'id': 12124623, 'label': 'All ligands | SMiles stereo', 'finished_at': None, 'created_at': '2024-01-20T03:46:13.000Z', 'updated_at': '2024-01-20T03:46:13.000Z', 'query_errors': None, 'finished_processing_at': None, 'query_type': 'STRUCTURE', 'fstruc_file_name': None, 'fstruc_content_type': None, 'fstruc_file_size': None, 'fstruc_updated_at': None, 'query_input': 'SF4\t[S]12[Fe]3[S]4[Fe]1[S]5[Fe]2[S]3[Fe]45\nFES\tS1[Fe]S[Fe]1\nACE\tNone\nSCM\tC[C@@H]1CC(=O)[C@]2([C@@H](O1)O[C@@H]3[C@H]([C@@H]([C@@H]([C@@H]([C@H]3O2)NC)O)NC)O)O\nON0\tc1ccc(cc1)[C@@H]2OC[C@@H]3[C@@H](O2)[C@@H]([C@H]([C@H](O3)O[C@@H]4[C@H](C[C@H]([C@@H]([C@H]4O[C@H]5[C@@H]([C@@H]([C@H](O5)CO)O[C@@H]6[C@@H]([C@H]([C@@H]([C@@H](O6)CN)O)O)N)O)O)N)N)N)O\n3TS\tc1cc(ccc1CO[C@@H]2[C@H](O[C@@H]([C@@H]([C@H]2O)N)O[C@@H]3[C@H](C[C@H]([C@@H]([C@H]3O[C@H]4[C@@H]([C@@H]([C@H](O4)CO)O[C@@H]5[C@@H]([C@H]([C@@H]([C@@H](O5)CN)O)O)N)O)O)N)N)CO)Cl\nAKN\tC1[C@@H]([C@H]([C@@H]([C@H]([C@@H]1NC(=O)[C@H](CCN)O)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)N)O)O)O[C@@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CN)O)O)O)N\n3KD\tc1c2c(cc3c1OCO3)[C@@H]4[C@@H]([C@H](C=C5[C@H]4N(C2)CC5)O)O\n3K5\tC[C@@H]1CO[C@]2(C[C@@H]1OC(=O)/C=C/c3ccccc3)[C@]4(CO4)[C@@H]5CC[C@@H](C[C@@H]5O2)C(=O)O[C@H]6[C@@H]([C@H]([C@@H]([C@H](O6)C)O)OC(=O)C)O[C@H]7[C@@H]([C@H]([C@@H]([C@H](O7)C)O)OC(=O)C)O\nA1AE1\tCC[C@@H]1[C@@]2([C@@H]([C@H](/C(=N/OC)/[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)OC(=O)NCCCCN3CCN(CC3)c4cc5c(cc4F)C(=O)C(=CN5C6CC6)C(=O)O)C)O[C@H]7[C@@H]([C@H](C[C@H](O7)C)N(C)C)O)(C)OC)C)C)OC(=O)O2)C\nSPR\tC[C@@H]1C\\C=C\\C=C\\[C@@H]([C@@H](C[C@@H]([C@@H]([C@H]([C@@H](CC(=O)O1)O)OC)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)C)O[C@H]3C[C@@]([C@H]([C@@H](O3)C)O)(C)O)N(C)C)O)CC=O)C)O[C@H]4CC[C@@H]([C@H](O4)C)N(C)C\n773\tCC[C@@H]1[C@@]2([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H](C(=O)[C@H](C(=O)O1)C)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)OC\\C=C\\c4cc5ccccc5nc4)C)C)NC(=O)O2)C\nPG4\tC(COCCOCCOCCO)O\nPEG\tC(COCCO)O\nEDO\tC(CO)O\nPGE\tC(COCCOCCO)O\n1PE\tC(COCCOCCOCCOCCO)O\nACY\tCC(=O)O\nGUN\tc1[nH]c2c(n1)C(=O)NC(=N2)N\nTRS\tC(C(CO)(CO)[NH3+])O\nYXM\tc1ccc(cc1)P(CCCCC(=O)N[C@H](CO)[C@@H](c2ccc(cc2)[N+](=O)[O-])O)(c3ccccc3)c4ccccc4\n6EM\tCC[C@@H](C(=O)N)[N+](C)(C)C\n95H\tC[C@H]([C@H]([C@@H]1[C@@H]([C@@H]([C@H]([C@H](O1)SC)O)O)O)NC(=O)c2ccc(cc2)[N+](=O)[O-])O\nOIY\tC1[C@H]([C@@H]2[C@@H](C(=O)N1)N=C(N2)N[C@H]3[C@@H]([C@@H]([C@H]([C@H](O3)CO)OC(=O)N)O)NC(=O)C[C@H](CCCNC(=O)C[C@H](CCCNC(=O)C[C@H](CCCN)N)N)N)O\nFSD\tC[C@@H]1[C@H](CC[C@@H](O1)N2C=CC(=NC2=O)NC(=O)c3ccc(cc3)NC(=O)[C@](C)(CO)N)O[C@@H]4[C@@H]([C@H]([C@@H]([C@H](O4)C)N(C)C)O)O\nWUX\tC[C@H]1[C@@H]2CC[C@]3([C@H]([C@]2(CC[C@H]1O)C)[C@@H](C[C@@H]4[C@@]3(C[C@@H](C4C(CCCC5CCCC5)C(=O)O)OC(=O)C)C)O)C\nZIT\tCC[C@@H]1[C@@]([C@@H]([C@H]([N@@](C[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)C)O)(C)O\nNEG\tC[N@@](CC(=O)O)NC(=O)C[C@@H](C[C@H](CN)O)N\nEZG\tc1cc(ccc1[C@H]([C@@H](CO)NC(=O)[C@H](Cc2c[nH]cn2)N)O)[N+](=O)[O-]\nG6V\tc1cc(c(cc1N2C[C@@H](OC2=O)CNC(=O)C(Cl)Cl)F)N3CCOCC3\nSY5\tCN(CC(=O)O)NC(=O)C[C@@H](C[C@H](CNCCCN)O)N\nA1AE0\tCC[C@@H]1[C@@]2([C@@H]([C@H](/C(=N/OC)/[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)OC(=O)NCCCCCC#Cc3ccc4c(c3)N(C=C(C4=O)C(=O)O)C5CC5)C)O[C@H]6[C@@H]([C@H](C[C@H](O6)C)N(C)C)O)(C)OC)C)C)OC(=O)O2)C\nO\tO\n5CR\tCC(=O)N[C@@H](Cc1ccccc1)C(=O)O\nWO2\t[O][W]1234O[W]567(O[W]89%10(O5[P]5%11O%12[W]%13%14(O6)(O[W]6%15(O1)(O[W]1%16%17(O6[P]6%18O2[W]2(O8)(O3)(O[W]38%19(O6[W](O3)(O1)(O[W]136(O5[W](O8)(O1)(O9)(O[W]15(O%10)(O%118[W](O7)(O%13)(O1)(O[W]8(O3)(O5)(O[W]%12(O%16)(O%14)(O6)[O])[O])[O])[O])[O])[O])(O[W]13(O%17)(O%185[W](O4)(O%15)(O1)(O[W]5(O2)(O%19)(O3)[O])[O])[O])[O])[O])[O])[O])[O])[O])[O])[O]\n3J2\tCC(C)C1=C2[C@H]([C@@H]3[C@H]4[C@]([C@H]([C@@H]5[C@H]([C@@]4(C2=CC(=O)O1)C)O5)O)(C(=O)O3)C)O\n5GP\tc1nc2c(n1[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N=C(NC2=O)N\nA\tc1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N\nC\tC1=CN(C(=O)N=C1N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)O)O)O\nYAT\tCC[C@H]1CN2CCc3cc(c(cc3[C@@H]2C[C@@H]1C[C@@H]4c5cc(c(cc5CCN4)OC)OC)OC)OC\nWC9\tC[C@@H]1/C=C\\CCS[C@@H]2[C@@H]([C@H]([C@H]([C@@H]([C@@H]1NC(=O)[C@@H]3[C@H]4[C@@H](C[C@@H](CCO4)CC(C)C)CN3)O2)O)O)O\nEDS\tC[C@@]1(CO[C@@H]([C@@H]([C@H]1NC)O)O[C@H]2[C@@H](C[C@@H]([C@H]([C@@H]2O)O[C@@H]3[C@@H](CC=C(O3)CNCCO)N)N)NC(=O)[C@H](CCN)O)O\nSFG\tc1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)C[C@H](CC[C@@H](C(=O)O)N)N)O)O)N\nVIR\tC[C@@H]1\\C=C\\C(=O)NC\\C=C\\C(=C\\[C@H](CC(=O)Cc2nc(co2)C(=O)N3CCC=C3C(=O)O[C@@H]1C(C)C)O)\\C\n3QB\tCCC[C@@H]1C[C@H](N(C1)C)C(=O)N[C@@H]([C@@H]2[C@@H]([C@@H]([C@H]([C@H](O2)SC)O)O)O)[C@@H](C)O\nAMP\tc1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N\n3LK\tCC[C@@]1(C[C@H]([C@@]2([C@@H](CC[C@@]3([C@H]2C(=O)CC3)[C@H]([C@@H]1O)C)C)C)OC(=O)CS[C@H]4CCCN(C4)C(=O)[C@H](C(C)C)N)C\nCAI\tC[C@@H]1C[C@@H]([C@@H]([C@H]([C@@H](CC(=O)O[C@@H](C[C@H]2[C@@H](O2)\\C=C\\C1=O)C)OC(=O)C)OC)O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)C)O[C@H]4C[C@@]([C@H]([C@@H](O4)C)OC(=O)CC(C)C)(C)O)N(C)C)O)CC=O\n84G\tC1C[C@H]([C@H](O[C@@H]1CN)O[C@@H]2[C@H](C[C@H]([C@@H]([C@H]2O)O[C@@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)N)O)NC(=O)[C@H](CCN)O)N)N\nHY0\tCN[C@H]1C[C@H]([C@@H]([C@H]([C@@H]1O)O[C@H]2[C@@H]3[C@H]([C@H]([C@H](O2)CO)O)O[C@@]4(O3)[C@@H]([C@H]([C@H]([C@H](O4)[C@H](CO)N)O)O)O)O)N\n5I0\tC[C@H]1[C@@]([C@H]([C@@H](O1)O[C@@H]2[C@H]([C@@H]([C@H]([C@@H]([C@H]2O)O)NC(=[NH2+])N)O)NC(=[NH2+])N)O[C@H]3[C@H]([C@@H]([C@H]([C@@H](O3)CO)O)O)[NH2+]C)(C(O)O)O\nMRD\tC[C@H](CC(C)(C)O)O\nT8B\tCOc1c2c(cc3c1c(c(c(c3)CC(=O)O)C(=O)/C=C(/c4ccccc4O)\\O)OC)cc5c(c2O)C(=O)OC(=C5)C(=O)OC\nRD8\tCC(=O)NC[C@H]1CN(C(=O)O1)c2ccc(c(c2)F)c3ccc(cc3)CNCc4c[nH]nn4\nW9C\tC[C@@H]1CC[C@@H]2[C@@H]1C[C@@]3([C@H]4C[C@H]([C@@]3([C@@]2(C4)C=O)C(=O)[O-])C(C)C)CO[C@@H]5[C@H]([C@@H]([C@@H]([C@H](O5)C)OC)O)O\nTIX\tC[C@@H]1CC(=O)[C@H]([C@]2([C@H]1C[C@@H]3[C@]45[C@@H]2[C@]([C@@H](C(=C)[C@@H]4CC(=O)O3)O)(OC5)O)C)O\nPHA\tc1ccc(cc1)C[C@@H](C=O)N\nMAU\tCC[C@H](C(=O)NC\\C=C\\C=C(/C)\\[C@H]([C@@H](C)[C@H]1[C@H]([C@H]([C@H](O1)\\C=C\\C=C\\C=C(/C)\\C(=O)C2=C(C=CN(C2=O)C)O)O)O)OC)[C@@]3([C@@H]([C@@H](C([C@@H](O3)\\C=C\\C=C/C)(C)C)O)O)O\nBME\tC(CS)O\nSLD\tCC\\1=NC(=O)NC(=O)/C1=C\\CC(=O)NCCC\\C=C\\c2ccc(cc2F)N3C[C@@H](OC3=O)CNC(=O)C\nMUL\tCCN(CC)CCSCC(=O)O[C@@H]1C[C@@]([C@H]([C@@H]([C@@]23CC[C@H]([C@@]1([C@@H]2C(=O)CC3)C)C)C)O)(C)C=C\nAQJ\tCC[C@@H]1[C@@](/C=C/C(=O)[C@@H](C[C@@H]([C@@H]([C@H](C(=O)[C@H](C(=O)O1)C)C)O[C@H]2[C@@H]([C@H](C[C@H](O2)C)N(C)C)O)C)C)(C)O\nSJE\tC[C@@H]1C[C@H]([C@H]([C@@H](O1)O[C@H]2[C@H](C[C@](C(=O)[C@@H]([C@@H]([C@H]([C@H](OC(=O)[C@@H]([C@H]([C@@H]2C)O[C@H]3C[C@@]([C@H]([C@@H](O3)C)O)(C)O)C)[C@@H](C)CO[C@H]4[C@@H]([C@@H]([C@@H]([C@H](O4)C)O)OC)OC)C)OC(=O)CC(C)C)C)(C)OC(=O)NCCNS(=O)(=O)c5ccccc5N(=O)=O)C)O)NOC\nCPT\t[NH3][Pt]([NH3])(Cl)Cl\nA1H4F\tCC(=O)O[C@@H]1[C@H](CN[C@@H]1Cc2ccc(cc2)OCC(=O)N[C@@H](CCCCN)C(=O)O)O\nACA\tC(CCC(=O)O)CCN\nUAM\tCC(C)C[C@@H]([C@@H]1Cc2cccc(c2C(=O)O1)O)NC(=O)[C@H]([C@H]([C@H](CC(=O)N)N)O)O\nU7V\tCn1nc(nn1)c2ccc(cn2)c3ccc(cc3F)N4C[C@@H](OC4=O)CO\n8Q1\tCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@@H](C(C)(C)COP(=O)(O)O)O\nWDP\tC[C@@H]1C[C@@H](O[C@H](CN1C)C)O[C@H]2[C@@H]([C@H]([C@H](C[C@](C(=O)[C@@H]([C@@H]([C@H]([C@H](OC(=O)[C@@H]2C)[C@@H](C)CO[C@H]3[C@@H]([C@@H]([C@@H]([C@H](O3)C)O)OC)OC)C)OC(=O)CC(C)C)C)(C)OC(=O)NC(C)(C)CNS(=O)(=O)c4ccccc4)C)O[C@H]5[C@@H](C(=NOC)C[C@H](O5)C)O)C\nU\tC1=CN(C(=O)NC1=O)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)O)O)O\nOMG\tCO[C@@H]1[C@@H]([C@H](O[C@H]1n2cnc3c2N=C(NC3=O)N)COP(=O)(O)O)O\nG\tc1nc2c(n1[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N=C(NC2=O)N\nPSU\tC1=C(C(=O)NC(=O)N1)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)O)O)O\nWIN\tCC1=C(C(=O)C[C@]2([C@H]1C[C@@H]3[C@]45[C@@H]2[C@H]([C@@H]([C@]([C@@H]4[C@H](C(=O)O3)OC(=O)/C=C(\\C)/C(C)C)(OC5)C(=O)OC)O)O)C)O\nMT9\tCC[C@@H]1[C@@](\\C=C\\C(=O)[C@@H](C[C@@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2[C@@H]([C@H](C[C@H](O2)C)N(C)C)O)C)C)(C)O\nA3P\tc1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)OP(=O)(O)O)O)N\n2AE\tc1ccc(c(c1)C(=O)N)N\nBGC\tNone\n3V6\tC[C@@H]1CC(=C2[C@@H]([C@H]1O)[C@H]([C@@](OC2=O)(C)C(Cl)Cl)NC(=O)[C@H](C)N)O\n6UQ\tCc1c(c(c(c(c1Cl)O)Cl)OC)C(=O)O[C@@H]2[C@H](O[C@H](C[C@H]2O)O[C@@H]3[C@H](O[C@]4(C[C@H]3O)O[C@@H]5[C@H](O[C@H](C[C@]5(O4)C)O[C@@H]6[C@H]([C@@H](O[C@@H]([C@@H]6OC)C)O[C@@H]7[C@H](O[C@H]([C@H]([C@H]7O)OC)O[C@H]8[C@@H]([C@H]9[C@H](CO8)O[C@@]1(O9)[C@H]2[C@H]([C@]([C@H](O1)C)([C@H](C)O)O)OCO2)OC(=O)C(C)C)COC)O)C)C)C\n7MB\tCN1C(=O)N[C@@H]2[C@]1(C[C@@H]3[C@H]2NC(=O)c4n3c(cc4)Br)O\n7AL\tC[C@]12C[C@H](CC([C@@H]1C[C@@H](C(=C)[C@@H]2C[C@@H]([C@H]3CC(=O)NC3=O)O)O)(C)C)Cl\n13T\tC\\C=C/[C@H](C)[C@@H]1[C@@](O1)(C)[C@H]([C@H]2COC(=O)[C@@H]([C@H]([C@@H](C(=O)[C@@H]([C@H](/C(=C/[C@@H](C(=O)CC[C@H](C2=O)C)C)/C)O)C)C)OC)O)O\nG34\tC[C@@H]1CC[C@@]23CCC(=O)[C@H]2[C@@]1([C@@H](C[C@@]([C@H]([C@@H]3C)O)(C)C=C)OC(=O)CSC4C[C@H]5CC[C@@H](C4)N5C)C\n6NO\tCc1c(c(c(c(c1Cl)O)Cl)OC)C(=O)O[C@@H]2[C@H](O[C@H](C[C@H]2O)O[C@@H]3[C@H](O[C@]4(C[C@H]3O)O[C@@H]5[C@H](O[C@H](C[C@]5(O4)C)O[C@@H]6[C@H]([C@@H](O[C@@H]([C@@H]6OC)C)O[C@@H]7[C@H](O[C@H]([C@H]([C@H]7O)OC)O[C@H]8[C@@H]([C@H]9[C@H](CO8)O[C@@]1(O9)[C@H]2[C@H]([C@@]([C@H](O1)C)(C(=O)C)O)OCO2)OC(=O)C(C)C)COC)O)C)C)C\n1F2\tCC[C@@H]1[C@@]2([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)OC(=O)N3CCC[C@@H]3c4cccnc4)C)O[C@H]5[C@@H]([C@H](C[C@H](O5)C)N(C)C)O)(C)OC)C)C)NC(=O)O2)C\n1F3\tCC[C@@H]1[C@@]2([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)OC(=O)N3C=CC[C@@H]3c4ccc(cc4)NC(=O)C)C)O[C@H]5[C@@H]([C@H](C[C@H](O5)C)N(C)C)O)(C)OC)C)C)NC(=O)O2)C\nV7A\tCN(C)[C@H]1[C@@H]2C[C@@H]3Cc4c(ccc(c4C(=O)C3=C([C@@]2(C(=O)C(=C1O)C(=O)N)O)O)O)CN(C)OC\nSPK\tC(CC[NH2+]CCC[NH3+])C[NH2+]CCC[NH3+]\nK16\tCC[C@H]1CN2CCc3cc(c(cc3[C@@H]2C[C@@H]1C[C@@H]4c5cc(c(cc5CCN4)O)OC)OC)OC\nPPU\tCN(C)c1c2c(ncn1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)NC(=O)[C@H](Cc4ccc(cc4)OC)N)O\n3J6\tCC1=C[C@@H]2[C@]([C@@H](C1=O)O)([C@]3(C[C@H]([C@H]([C@@]34CO4)O2)O)C)CO\nEVN\tCc1cc(cc(c1C(=O)O[C@@H]2CO[C@]3([C@H]4[C@H]2OCO4)O[C@H]5CO[C@H]([C@@H]([C@@H]5O3)O)O[C@H]6[C@H]([C@H]([C@@H]([C@H](O6)COC)O[C@H]7[C@@H]([C@H]([C@H]([C@H](O7)C)OC)O[C@H]8[C@H]([C@@]9([C@@H]([C@H](O8)C)O[C@@]1(O9)C[C@H]([C@@H]([C@H](O1)C)O[C@H]1C[C@H]([C@@H]([C@H](O1)C)OC(=O)c1c(c(c(c(c1OC)Cl)O)Cl)C)O[C@@H]1C[C@]([C@H]([C@@H](O1)C)OC)(C)[N+](=O)[O-])O)C)O)O)O)OC)O)O\nN\tC1[C@@H]([C@@H]([C@H](O1)COP(=O)(O)O)O)O\n3K8\tCOc1ccc2c(c1)c3cc(c(cc3c4c2CN5CCCC[C@@H]5C4)OC)OC\n4M2\tC[C@@H]1[C@H]([C@@H]([C@@H]([C@H](O1)OC/C(=C/2\\[C@H]([C@@H]([C@H](O2)Oc3ccc(cc3)/C=C(\\C)/C(=O)N[C@@H]4[C@H](O[C@H]([C@@H]4O)n5cnc6c5ncnc6N(C)C)CO)O)O)/OC)O)OC)OC\nGOL\tC(C(CO)O)O\n80P\tCCN(CC)[C@H]1[C@@H]2C[C@@H]3Cc4c(c(cc(c4C(F)(F)F)[C@@H]5CCCN5)O)C(=O)C3=C([C@@]2(C(=O)C(=C1O)C(=O)N)O)O\n917\tCC(=O)NC[C@H]1CN(C(=O)O1)c2ccc(cc2)c3cncs3\n6O1\tCc1cc(cc(c1C(=O)O[C@@H]2CO[C@]3([C@H]4[C@H]2OCO4)O[C@H]5CO[C@H]([C@@H]([C@@H]5O3)O)O[C@H]6[C@H]([C@H]([C@@H]([C@H](O6)COC)O[C@H]7[C@@H]([C@H]([C@H]([C@H](O7)C)OC)O[C@H]8[C@H]([C@@]9([C@@H]([C@H](O8)C)O[C@@]1(O9)C[C@H]([C@@H]([C@H](O1)C)O[C@H]1C[C@H]([C@@H]([C@H](O1)C)OC(=O)c1c(c(c(c(c1OC)Cl)O)Cl)C)O[C@@H]1C[C@@](C(C(O1)C)OC)(C)[N+](=O)[O-])O)C)O)O)O)OC)O)O\nEUS\tC[C@@]1(CO[C@@H]([C@@H]([C@H]1NC)O)O[C@H]2[C@@H](C[C@@H]([C@H]([C@@H]2O)O[C@@H]3[C@@H](CC=C(O3)CN)N)N)NS(=O)(=O)C)O\nSJH\tCc1ncc(n1CCNC(=O)O[C@]2(C[C@@H]([C@@H]([C@H]([C@@H]([C@H](C(=O)O[C@@H]([C@@H]([C@H]([C@H](C2=O)C)OC(=O)CC(C)C)C)[C@@H](C)CO[C@H]3[C@@H]([C@@H]([C@@H]([C@H](O3)C)O)OC)OC)C)O[C@H]4C[C@@]([C@H]([C@@H](O4)C)O)(C)O)C)O[C@H]5[C@@H](/C(=N/OC)/C[C@H](O5)C)O)C)C)[N+](=O)[O-]\nZIY\tCC[C@H](C)[C@@H]1C(CC(=O)O[C@H](C(=O)[C@@H](C(=O)N[C@H](C(=O)N2CCC[C@H]2C(=O)N([C@H](C(=O)O[C@@H]([C@@H](C(=O)N1)NC(=O)[C@@H](CC(C)C)N(C)C(=O)[C@@H]3CCCN3C(=O)C(=O)C)C)Cc4ccc(cc4)OC)C)CC(C)C)C)C(C)C)O\nM2D\tC[C@@H]1/C=C/C(=O)NC/C=C\\C(=C\\[C@H](C[C@@H](Cc2nc(co2)C(=O)N[C@@H](C(=O)O[C@@H]1C(C)C)C)O)O)\\C\nEM1\tCC[C@@H]1[C@@]2([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H](C(=O)[C@](C(=O)O1)(C)F)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)OC)C)C)N(C(=O)O2)CCCCn4cc(nn4)c5cccc(c5)N)C\nLEU\tCC(C)C[C@@H](C(=O)O)N\nA2G\tCC(=O)N[C@@H]1[C@H]([C@H]([C@H](O[C@@H]1O)CO)O)O\nORN\tC(C[C@@H](C(=O)O)N)CN\nHJO\tCC(=O)N[C@@H]1C[C@@H]([C@@H]([C@H]([C@@H]1O[C@@H]2[C@@H]([C@H](C(CO2)(C)O)NC)O)O)O[C@@H]3[C@@H](CC[C@H](O3)CN)N)N\nPEV\tCCCCCCCCCCCCCCCCCC(=O)O[C@@H](COC(=O)CCCCCCCCCCCCCCC)CO[P@@](=O)(O)OCCN\nPGV\tCCCCCCCCCCCCCCCC(=O)OC[C@H](CO[P@](=O)(O)OC[C@H](CO)O)OC(=O)CCCCCCCCC\\C=C/CCCCCC\nD2C\tCN(C)[C@H]1[C@@H]2C[C@@H]3[C@@H](c4c(ccc(c4C(C3C([C@@H]2C(=O)C(C1=O)C(=O)N)O)O)O)Cl)O\nAB9\tC1[C@@H]([C@H]([C@@H]([C@H]([C@@H]1NC(=O)[C@@H](CCN)O)OCCNCCN)O)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CN)O)O)N)N\n3AB\tc1cc(cc(c1)N)C(=O)N\nTAO\tC[C@@H]1C[C@@H]([C@H]([C@@H](O1)O[C@H]2[C@@H](C[C@@]3(CO3)C(=O)[C@@H]([C@@H]([C@H]([C@H](OC(=O)[C@@H]([C@H]([C@@H]2C)O[C@H]4C[C@@H]([C@H]([C@@H](O4)C)OC(=O)C)OC)C)C)C)OC(=O)C)C)C)OC(=O)C)N(C)C\nRPO\tc1ccc(cc1)CO[C@@H]2[C@H](O[C@@H]([C@@H]([C@H]2O)N)O[C@@H]3[C@H](C[C@H]([C@@H]([C@H]3O[C@H]4[C@@H]([C@@H]([C@H](O4)CO)O[C@@H]5[C@@H]([C@H]([C@@H]([C@@H](O5)CN)O)O)N)O)O)N)N)CO\nY7K\tC[C@@H]1C[C@]([C@]2([C@@H](O1)O[C@@H]3[C@H]([C@@H]([C@@H]([C@@H]([C@H]3O2)NC)O)NC)O)O)(CNCCC4CCOCC4)O\nNO1\tCc1c2c(cccc2[nH]c1C(=O)O)CO\nEZP\tc1cc(ccc1[C@H]([C@@H](CO)NC(=O)[C@@H](Cc2c[nH]cn2)N)O)N(=O)=O\n6IF\tC[C@@H]([C@H]([C@@H]1[C@@H]([C@@H]([C@H]([C@H](O1)SC)O)O)O)NC(=O)[C@@H]2[C@H]3[C@@H](C[C@@H](CCO3)CC(C)C)CN2)Cl\nA1AIX\tCOc1ccc2c(c1OC)c[n+]3c(c2CC(=O)N[C@H](CO)[C@@H](c4ccc(cc4)[N+](=O)[O-])O)-c5cc6c(cc5CC3)OCO6\nG6M\tc1cc(c(cc1N2C[C@@H](OC2=O)CNC(=O)CCl)F)N3CCOCC3\n3KF\tc1c2c(c(c3c1OCO3)O)C(=O)N[C@@H]4C2=C[C@@H]([C@H]([C@H]4O)O)O\nOCW\tCc1c2c(cc3c(c2O)C(=O)[C@@]4(C(=O)C=C([C@H]([C@@]4(C3=O)O)O)OC)OC)cc(c1C(=O)OC)OC\nEZM\tc1cc(ccc1[C@H]([C@@H](CO)NC(=O)[C@H](CCCCN)N)O)[N+](=O)[O-]\nHJR\tC[C@@]1(CO[C@@H](C[C@@H]1NC)O[C@@H]2[C@@H](C[C@H]([C@@H]([C@H]2O)O[C@@H]3CCC[C@H](O3)CN)N)NS(=O)(=O)C(F)(F)F)O\nKAN\tC1[C@H]([C@@H]([C@H]([C@@H]([C@H]1N)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CN)O)O)O)O)O[C@@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)N)O)N\nZBA\tCC1=C[C@@H]2[C@](C[C@@H]1OC(=O)CC(=C)C)([C@]3([C@@H]([C@H]([C@H]([C@@]34CO4)O2)O)OC(=O)C)C)COC(=O)C\n1I7\tCC[C@@H]1C(=O)NC(=C)C(=O)N(CC(=O)N[C@@H](c2nc(cs2)C(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N1)Cc3c[nH]c4c3c(ccc4)OC)Cc5c[nH]c6c5cccc6)C)C\nCTY\tCC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)OC)C)C)O)(C)O\nBTN\tC1[C@H]2[C@@H]([C@@H](S1)CCCCC(=O)O)NC(=O)N2\nTOY\tC1[C@@H]([C@H]([C@@H]([C@H]([C@@H]1N)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)N)O)O)O[C@@H]3[C@@H](C[C@@H]([C@H](O3)CN)O)N)N\nEMK\tCC[C@@H]1[C@@]([C@@H]([C@H]([N@](C[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)[N@](C)CCCc4cn(nn4)[C@H](CF)[C@@H](c5ccc(cc5)S(=O)(=O)C)O)O)(C)O)C)C)C)O)(C)O\nKKL\tC#Cc1ccc(cc1)c2nnc(o2)NC(=O)c3ccc(cc3)N\nHN8\tCO[C@H]1C[C@H]2[C@]3(C=C1)c4cc5c(cc4CN2C[C@@H]3O)OCO5\nSPE\tC(CN)CNCCCNCCCN\nP8F\tCN(C)[C@H]1[C@@H]2C[C@@H]3Cc4cc5ccc(cc5c(c4C(=O)C3=C([C@@]2(C(=O)C(=C1O)C(=O)N)O)O)O)CN6CC(C6)F\nH8T\tC[C@@H]1/C=C/C(=O)NC/C=C/C(=C/[C@H](CC(=O)Cc2nc(co2)C(=O)n3cccc3C(=O)O[C@@H]1C(C)C)O)/C\nCYS\tC([C@@H](C(=O)O)N)S\nGLN\tC(CC(=O)N)[C@@H](C(=O)O)N\nGIR\tc1c([nH]c(n1)N)[C@@H]([C@H](CN)Cl)O\nCPF\tc1c2c(cc(c1F)N3CCNCC3)N(C=C(C2=O)C(=O)O)C4CC4\nUNL\tNone\nFDA\tCc1cc2c(cc1C)N(C3=C(N2)C(=O)NC(=O)N3)C[C@@H]([C@@H]([C@@H](CO[P@@](=O)(O)O[P@@](=O)(O)OC[C@@H]4[C@H]([C@H]([C@@H](O4)n5cnc6c5ncnc6N)O)O)O)O)O\n7NO\tCSCC[C@@H](C(=O)O[C@@H]1[C@H](O[C@H]([C@@H]1O)n2cnc3c2ncnc3N)COP(=O)(O)O)N\nGGM\tCNc1ccccc1C(=O)O[C@@H]2[C@H](O[C@H]([C@@H]2O)n3cnc4c3N=C(NC4=O)N)CO[P@](=O)(O)O[P@](=O)(NP(=O)(O)O)O\nU7Y\tCn1c(nnn1)c2ccc(cn2)c3ccc(cc3F)N4C[C@@H](OC4=O)CO\nJY7\tCc1cnccc1c2c3c4cc(cc(c4nnc3n(n2)C)OC)F\n7YT\tC[C@@H]1CC(=O)N(N=C1c2ccc(cc2)NC(=O)N3Cc4cccnc4C3)C\nA1A1L\tC[C@@H]1/C=C\\C[C@H](CS[C@@H]2[C@@H]([C@H]([C@H]([C@@H]([C@@H]1NC(=O)[C@@H]3[C@H]4[C@@H](C[C@@H](CCO4)CC(C)C)CN3)O2)O)O)O)F\n', 'tag_list': ['All ligands', 'SMiles stereo']}


# 3/3 : {'id': 12124652, 'label': 'All ligands | SMiles stereo', 'finished_at': None, 'created_at': '2024-01-20T03:50:00.000Z', 'updated_at': '2024-01-20T03:50:00.000Z', 'query_errors': None, 'finished_processing_at': None, 'query_type': 'STRUCTURE', 'fstruc_file_name': None, 'fstruc_content_type': None, 'fstruc_file_size': None, 'fstruc_updated_at': None, 'query_input': 'SF4\t[S]12[Fe]3[S]4[Fe]1[S]5[Fe]2[S]3[Fe]45\nFES\tS1[Fe]S[Fe]1\nACE\tNone\nWO2\t[O][W]1234O[W]567(O[W]89%10(O5[P]5%11O%12[W]%13%14(O6)(O[W]6%15(O1)(O[W]1%16%17(O6[P]6%18O2[W]2(O8)(O3)(O[W]38%19(O6[W](O3)(O1)(O[W]136(O5[W](O8)(O1)(O9)(O[W]15(O%10)(O%118[W](O7)(O%13)(O1)(O[W]8(O3)(O5)(O[W]%12(O%16)(O%14)(O6)[O])[O])[O])[O])[O])[O])(O[W]13(O%17)(O%185[W](O4)(O%15)(O1)(O[W]5(O2)(O%19)(O3)[O])[O])[O])[O])[O])[O])[O])[O])[O])[O])[O]\nBGC\tNone\nSJH\tCc1ncc(n1CCNC(=O)O[C@]2(C[C@@H]([C@@H]([C@H]([C@@H]([C@H](C(=O)O[C@@H]([C@@H]([C@H]([C@H](C2=O)C)OC(=O)CC(C)C)C)[C@@H](C)CO[C@H]3[C@@H]([C@@H]([C@@H]([C@H](O3)C)O)OC)OC)C)O[C@H]4C[C@@]([C@H]([C@@H](O4)C)O)(C)O)C)O[C@H]5[C@@H](/C(=N/OC)/C[C@H](O5)C)O)C)C)[N+](=O)[O-]\nZIY\tCC[C@H](C)[C@@H]1C(CC(=O)O[C@H](C(=O)[C@@H](C(=O)N[C@H](C(=O)N2CCC[C@H]2C(=O)N([C@H](C(=O)O[C@@H]([C@@H](C(=O)N1)NC(=O)[C@@H](CC(C)C)N(C)C(=O)[C@@H]3CCCN3C(=O)C(=O)C)C)Cc4ccc(cc4)OC)C)CC(C)C)C)C(C)C)O\nM2D\tC[C@@H]1/C=C/C(=O)NC/C=C\\C(=C\\[C@H](C[C@@H](Cc2nc(co2)C(=O)N[C@@H](C(=O)O[C@@H]1C(C)C)C)O)O)\\C\nEM1\tCC[C@@H]1[C@@]2([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H](C(=O)[C@](C(=O)O1)(C)F)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)OC)C)C)N(C(=O)O2)CCCCn4cc(nn4)c5cccc(c5)N)C\nLEU\tCC(C)C[C@@H](C(=O)O)N\nA2G\tCC(=O)N[C@@H]1[C@H]([C@H]([C@H](O[C@@H]1O)CO)O)O\nORN\tC(C[C@@H](C(=O)O)N)CN\nHJO\tCC(=O)N[C@@H]1C[C@@H]([C@@H]([C@H]([C@@H]1O[C@@H]2[C@@H]([C@H](C(CO2)(C)O)NC)O)O)O[C@@H]3[C@@H](CC[C@H](O3)CN)N)N\nPEV\tCCCCCCCCCCCCCCCCCC(=O)O[C@@H](COC(=O)CCCCCCCCCCCCCCC)CO[P@@](=O)(O)OCCN\nPGV\tCCCCCCCCCCCCCCCC(=O)OC[C@H](CO[P@](=O)(O)OC[C@H](CO)O)OC(=O)CCCCCCCCC\\C=C/CCCCCC\nD2C\tCN(C)[C@H]1[C@@H]2C[C@@H]3[C@@H](c4c(ccc(c4C(C3C([C@@H]2C(=O)C(C1=O)C(=O)N)O)O)O)Cl)O\nAB9\tC1[C@@H]([C@H]([C@@H]([C@H]([C@@H]1NC(=O)[C@@H](CCN)O)OCCNCCN)O)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CN)O)O)N)N\n3AB\tc1cc(cc(c1)N)C(=O)N\nTAO\tC[C@@H]1C[C@@H]([C@H]([C@@H](O1)O[C@H]2[C@@H](C[C@@]3(CO3)C(=O)[C@@H]([C@@H]([C@H]([C@H](OC(=O)[C@@H]([C@H]([C@@H]2C)O[C@H]4C[C@@H]([C@H]([C@@H](O4)C)OC(=O)C)OC)C)C)C)OC(=O)C)C)C)OC(=O)C)N(C)C\nRPO\tc1ccc(cc1)CO[C@@H]2[C@H](O[C@@H]([C@@H]([C@H]2O)N)O[C@@H]3[C@H](C[C@H]([C@@H]([C@H]3O[C@H]4[C@@H]([C@@H]([C@H](O4)CO)O[C@@H]5[C@@H]([C@H]([C@@H]([C@@H](O5)CN)O)O)N)O)O)N)N)CO\nY7K\tC[C@@H]1C[C@]([C@]2([C@@H](O1)O[C@@H]3[C@H]([C@@H]([C@@H]([C@@H]([C@H]3O2)NC)O)NC)O)O)(CNCCC4CCOCC4)O\nNO1\tCc1c2c(cccc2[nH]c1C(=O)O)CO\nEZP\tc1cc(ccc1[C@H]([C@@H](CO)NC(=O)[C@@H](Cc2c[nH]cn2)N)O)N(=O)=O\n6IF\tC[C@@H]([C@H]([C@@H]1[C@@H]([C@@H]([C@H]([C@H](O1)SC)O)O)O)NC(=O)[C@@H]2[C@H]3[C@@H](C[C@@H](CCO3)CC(C)C)CN2)Cl\nA1AIX\tCOc1ccc2c(c1OC)c[n+]3c(c2CC(=O)N[C@H](CO)[C@@H](c4ccc(cc4)[N+](=O)[O-])O)-c5cc6c(cc5CC3)OCO6\nG6M\tc1cc(c(cc1N2C[C@@H](OC2=O)CNC(=O)CCl)F)N3CCOCC3\n3KF\tc1c2c(c(c3c1OCO3)O)C(=O)N[C@@H]4C2=C[C@@H]([C@H]([C@H]4O)O)O\nOCW\tCc1c2c(cc3c(c2O)C(=O)[C@@]4(C(=O)C=C([C@H]([C@@]4(C3=O)O)O)OC)OC)cc(c1C(=O)OC)OC\nEZM\tc1cc(ccc1[C@H]([C@@H](CO)NC(=O)[C@H](CCCCN)N)O)[N+](=O)[O-]\nHJR\tC[C@@]1(CO[C@@H](C[C@@H]1NC)O[C@@H]2[C@@H](C[C@H]([C@@H]([C@H]2O)O[C@@H]3CCC[C@H](O3)CN)N)NS(=O)(=O)C(F)(F)F)O\nKAN\tC1[C@H]([C@@H]([C@H]([C@@H]([C@H]1N)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CN)O)O)O)O)O[C@@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)N)O)N\nZBA\tCC1=C[C@@H]2[C@](C[C@@H]1OC(=O)CC(=C)C)([C@]3([C@@H]([C@H]([C@H]([C@@]34CO4)O2)O)OC(=O)C)C)COC(=O)C\n1I7\tCC[C@@H]1C(=O)NC(=C)C(=O)N(CC(=O)N[C@@H](c2nc(cs2)C(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N1)Cc3c[nH]c4c3c(ccc4)OC)Cc5c[nH]c6c5cccc6)C)C\nCTY\tCC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)OC)C)C)O)(C)O\nBTN\tC1[C@H]2[C@@H]([C@@H](S1)CCCCC(=O)O)NC(=O)N2\nTOY\tC1[C@@H]([C@H]([C@@H]([C@H]([C@@H]1N)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)N)O)O)O[C@@H]3[C@@H](C[C@@H]([C@H](O3)CN)O)N)N\nEMK\tCC[C@@H]1[C@@]([C@@H]([C@H]([N@](C[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)[N@](C)CCCc4cn(nn4)[C@H](CF)[C@@H](c5ccc(cc5)S(=O)(=O)C)O)O)(C)O)C)C)C)O)(C)O\nKKL\tC#Cc1ccc(cc1)c2nnc(o2)NC(=O)c3ccc(cc3)N\nHN8\tCO[C@H]1C[C@H]2[C@]3(C=C1)c4cc5c(cc4CN2C[C@@H]3O)OCO5\nSPE\tC(CN)CNCCCNCCCN\nP8F\tCN(C)[C@H]1[C@@H]2C[C@@H]3Cc4cc5ccc(cc5c(c4C(=O)C3=C([C@@]2(C(=O)C(=C1O)C(=O)N)O)O)O)CN6CC(C6)F\nH8T\tC[C@@H]1/C=C/C(=O)NC/C=C/C(=C/[C@H](CC(=O)Cc2nc(co2)C(=O)n3cccc3C(=O)O[C@@H]1C(C)C)O)/C\nCYS\tC([C@@H](C(=O)O)N)S\nGLN\tC(CC(=O)N)[C@@H](C(=O)O)N\nGIR\tc1c([nH]c(n1)N)[C@@H]([C@H](CN)Cl)O\nCPF\tc1c2c(cc(c1F)N3CCNCC3)N(C=C(C2=O)C(=O)O)C4CC4\nUNL\tNone\nFDA\tCc1cc2c(cc1C)N(C3=C(N2)C(=O)NC(=O)N3)C[C@@H]([C@@H]([C@@H](CO[P@@](=O)(O)O[P@@](=O)(O)OC[C@@H]4[C@H]([C@H]([C@@H](O4)n5cnc6c5ncnc6N)O)O)O)O)O\n7NO\tCSCC[C@@H](C(=O)O[C@@H]1[C@H](O[C@H]([C@@H]1O)n2cnc3c2ncnc3N)COP(=O)(O)O)N\nGGM\tCNc1ccccc1C(=O)O[C@@H]2[C@H](O[C@H]([C@@H]2O)n3cnc4c3N=C(NC4=O)N)CO[P@](=O)(O)O[P@](=O)(NP(=O)(O)O)O\nU7Y\tCn1c(nnn1)c2ccc(cn2)c3ccc(cc3F)N4C[C@@H](OC4=O)CO\nJY7\tCc1cnccc1c2c3c4cc(cc(c4nnc3n(n2)C)OC)F\n7YT\tC[C@@H]1CC(=O)N(N=C1c2ccc(cc2)NC(=O)N3Cc4cccnc4C3)C\nA1A1L\tC[C@@H]1/C=C\\C[C@H](CS[C@@H]2[C@@H]([C@H]([C@H]([C@@H]([C@@H]1NC(=O)[C@@H]3[C@H]4[C@@H](C[C@@H](CCO4)CC(C)C)CN3)O2)O)O)O)F\n', 'tag_list': ['All ligands', 'SMiles stereo']}