import asyncio
from functools import reduce
from io import StringIO
from itertools import tee
import json
import logging
import os
from pprint import pprint
from tempfile import NamedTemporaryFile
from typing import Iterator
from Bio.Align import MultipleSeqAlignment,Seq, SeqRecord
from Bio.Align.Applications import MuscleCommandline
from ribctl import ASSETS, RIBETL_DATA
from ribctl.etl.etl_pipeline import ReannotationPipeline, query_rcsb_api, rcsb_single_structure_graphql
from ribctl.lib.classification import classify_sequence, classify_subchains
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.tunnel import ptc_resdiues_get, ptc_residues_calculate_midpoint
from ribctl import model_species, model_subgenuses
from ete3 import NCBITaxa

logging.getLogger("urllib3.connectionpool").setLevel(logging.CRITICAL)
logging.getLogger('asyncio').setLevel(logging.WARNING)
from ribctl.lib.util_taxonomy import descendants_of_taxid
hmm_cachedir = ASSETS['__hmm_cache']
import sys

logger       = logging.getLogger(__name__)
file_handler = logging.FileHandler('classification.log')
log_format   = logging.Formatter('%(asctime)s [%(levelname)s] [%(name)s] %(message)s')
file_handler.setFormatter(log_format)
logger.addHandler(file_handler)

nomv2dict= '/home/rtviii/dev/riboxyz/nomv2'

if sys.argv[1]    == "process_struct":
   rcsb_id         = sys.argv[2].upper()
   rib             = RibosomeAssets(rcsb_id).profile()
   organism_taxid  = rib.src_organism_ids[0]
   prots           = rib.proteins
   result          = classify_subchains(prots)
   print(result)
elif sys.argv[1] == "process_all":
    for rcsb_id in RibosomeAssets.list_all_structs():
        logger.debug("Processing {}".format(rcsb_id))
        rib            = RibosomeAssets(rcsb_id).profile()
        organism_taxid = rib.src_organism_ids[0]
        prots          = rib.proteins
        results        = classify_subchains(prots)
        with open('/home/rtviii/dev/riboxyz/nomv2/{}.json'.format(rcsb_id), 'w') as f:
            json.dump(results, f)


elif sys.argv[1] =="merge_nomenclature":


    for f in os.listdir(nomv2dict):
        
        print(f)
        with open(os.path.join(nomv2dict, f),'r') as infile:
            nomv2=json.load(infile)
        struct, _ = f.split(".")
        rib_asset = RibosomeAssets(struct)


        filepath = rib_asset._json_profile_filepath()
        profile  = rib_asset.profile()
        
        for chain in profile.proteins:
            if chain.auth_asym_id in nomv2:
                if nomv2[chain.auth_asym_id] != None and nomv2[chain.auth_asym_id] not in chain.nomenclature:
                    chain.nomenclature = [nomv2[chain.auth_asym_id]]
                elif nomv2[chain.auth_asym_id] == None and chain.nomenclature != []:
                    chain.nomenclature = []

        rib_asset.write_own_json_profile(new_profile=json.loads(profile.json()), overwrite=True)


elif sys.argv[1] =="tunnel":
    def list_euk_structs():
        EUK_STRUCTS= []
        with open("eukarya_2023.txt", "r") as data_file:
            for line in data_file:
                structs = line.split(",")
                EUK_STRUCTS = [*EUK_STRUCTS, *structs]
        return EUK_STRUCTS

    if __name__ == "__main__":
        EUK        = list_euk_structs()
   
    print("tunnel")
elif sys.argv[1] == "spec":
    def get_taxonomic_id(organism_name):
        ncbi = NCBITaxa()

        try:
            # Get the taxonomic ID for the given organism name
            taxid = ncbi.get_name_translator([organism_name])[organism_name]
            return taxid
        except KeyError:
            # Handle the case where the organism name is not found
            print(f"Organism '{organism_name}' not found in the NCBI Taxonomy database.")
            return None


    _ ={}
    for organism_name in ['Lactococcus lactis','Mycobacterium smegmatis', 'Mycobacterium tuberculosis', 'Bacillus subtilis'
                          , 'Leishmania donovani', 'Trypanosoma cruzi', 'Trichomonas vaginalis','Giardia duodenalis','Spraguea lophii' ]:
        taxonomic_id = get_taxonomic_id(organism_name)
        
        if taxonomic_id:
            _ = {**_, **{organism_name: taxonomic_id}}
            print(f"Taxonomic ID for {organism_name}: {taxonomic_id}")
    print(_)

elif sys.argv[1] == "test":

    ncbi = NCBITaxa()
    # rp = RibosomeAssets('5MYJ').profile()
    # print(( rp.rcsb_id, rp.src_organism_ids, rp.host_organism_ids ))
    # print(( rp.rcsb_id, rp.src_organism_names, rp.host_organism_names))
    # print(NCBITaxa().get_lineage(rp.src_organism_ids[0]))
    # names=  ncbi.translate_to_names(NCBITaxa().get_lineage(rp.src_organism_ids[0]))
    # print(names)

    # print(descendants_of_taxid( pdbid_taxid_tuples, int(taxid)))

    all_structs = os.listdir(RIBETL_DATA)
    pdbid_taxid_tuples:list = []    

    for struct in all_structs:
        print(struct)
        rp = RibosomeAssets(struct)
        asyncio.run(rp._verify_json_profile())
        rp = rp.profile()
        pdbid_taxid_tuples.append(( rp.rcsb_id, rp.src_organism_ids[0] ))

    for (name, subgenus_taxid) in ( model_subgenuses.items() ):
        print("\n\t###########      Descendants of {} ({})      ##########3".format(subgenus_taxid, ncbi.get_taxid_translator([subgenus_taxid])[subgenus_taxid]))
        rcsb_id_taxid_tuples =  descendants_of_taxid( pdbid_taxid_tuples, subgenus_taxid )

        for i,(rcsb_id, taxid) in enumerate( rcsb_id_taxid_tuples ):
            print("{}.".format(i),rcsb_id, taxid, ncbi.get_taxid_translator([taxid])[taxid])

        ids_in_this_branch = [*list(map(lambda x: x[1], rcsb_id_taxid_tuples)), subgenus_taxid]
        tree =ncbi.get_topology(ids_in_this_branch)
        print(tree.get_ascii(attributes=["taxid"]))


elif sys.argv[1] == "processtax":
    RCSB_ID = '5MYJ'
    ReannotationPipeline(query_rcsb_api(rcsb_single_structure_graphql(RCSB_ID))).process_structure()