import asyncio
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence, SequenceFile, SequenceBlock, TextSequenceBlock
from pyhmmer.plan7 import Pipeline, HMM 
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
from pyhmmer import phmmer
import pyhmmer
from ribctl import ASSETS, RIBETL_DATA
from ribctl.etl.etl_pipeline import ReannotationPipeline, query_rcsb_api, rcsb_single_structure_graphql
from ribctl.etl.obtain import obtain_assets_threadpool
from ribctl.lib.classification import classify_sequence, classify_subchains, hmm_create, hmm_dict_init__candidates_per_organism, hmm_produce
from ribctl.etl.ribosome_assets import Assetlist, RibosomeAssets
from ribctl.lib.ribosome_types.types_ribosome import RNAClassEnum
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
   rnas            = rib.rnas
   result          = classify_subchains(rnas, RNAClassEnum)
   print(result)
elif sys.argv[1] == "process_all":
    for rcsb_id in RibosomeAssets.list_all_structs():
        logger.debug("Processing {}".format(rcsb_id))
        rib            = RibosomeAssets(rcsb_id).profile()
        organism_taxid = rib.src_organism_ids[0]
        rnas          = rib.rnas
        results        = classify_subchains(prots)
        print(results)
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
    all_structs = os.listdir(RIBETL_DATA)
    pdbid_taxid_tuples:list = []    
    __structs = []


    for i,struct in enumerate(all_structs):
        try:
            print(struct, i)
            rp = RibosomeAssets(struct)
            rp = rp.profile()
            pdbid_taxid_tuples.append(( rp.rcsb_id, rp.src_organism_ids[0] ))
        except:
            ...
    for (name, subgenus_taxid) in ( model_subgenuses.items() ):
        print("\n\t###########      Descendants of {} ({})      ##########3".format(subgenus_taxid, ncbi.get_taxid_translator([subgenus_taxid])[subgenus_taxid]))
        rcsb_id_taxid_tuples =  descendants_of_taxid( pdbid_taxid_tuples, subgenus_taxid )

        for i,(rcsb_id, taxid) in enumerate( rcsb_id_taxid_tuples ):
            __structs.append(rcsb_id)
            print("{}.".format(i),rcsb_id, taxid, ncbi.get_taxid_translator([taxid])[taxid])


    print(__structs)
        
        # ids_in_this_branch = [*list(map(lambda x: x[1], rcsb_id_taxid_tuples)), subgenus_taxid]
        # tree =ncbi.get_topology(ids_in_this_branch)
        # print(tree.get_ascii(attributes=["taxid"]))

        # obtain_assets_threadpool([rcsb_id for (rcsb_id, taxid) in rcsb_id_taxid_tuples], Assetlist(ptc_coords=True), overwrite=True)

elif sys.argv[1] == "processtax":
    RCSB_ID = '5MYJ'
    ReannotationPipeline(query_rcsb_api(rcsb_single_structure_graphql(RCSB_ID))).process_structure()
elif sys.argv[1] == 'll':
    import shutil

    structs = ['5MYJ', '6DZI', '5XYM', '5ZEB', '7S0S', '5O61', '5ZET', '7Y41', '5ZEP', '6DZK', '5ZEU', '5O60', '5O5J', '6DZP', '7XAM', '5XYU', '7MT2', '5V7Q', '7MSM', '7MT7', '5V93', '7MSC', '7KGB', '7MSH', '7F0D', '7MSZ', '7MT3', '7SFR', '7AQC', '8BUU', '7SAE', '7S9U', '7QV2', '6PPK', '7O5B', '7QGU', '3J3W', '6HA1', '6HTQ', '6HA8', '7OPE', '3J3V', '3J9W', '7QV1', '7QV3', '7QH4', '6PPF', '7AQD', '6PVK', '6AZ3', '7ANE', '5T2A', '3JCS', '6AZ1', '7AIH', '7AOR', '5OPT', '5T5H', '5XY3', '5XYI', '8BTR', '7PWO', '8BR8', '7PWG', '8BSJ', '8BSI', '7PWF', '8BTD', '8BRM', '7QCA', '8P60', '8P5D']

    ncbi = NCBITaxa()
    all_structs = os.listdir(RIBETL_DATA)
    pdbid_taxid_tuples:list = []    
    __structs = []

    for i,struct in enumerate(all_structs):
        try:
            print(struct, i)
            rp = RibosomeAssets(struct)
            rp = rp.profile()
            pdbid_taxid_tuples.append(( rp.rcsb_id, rp.src_organism_ids[0] ))
        except:
            ...

    for (name, subgenus_taxid) in ( model_subgenuses.items() ):
        print("\n\t###########      Descendants of {} ({})      ##########3".format(subgenus_taxid, ncbi.get_taxid_translator([subgenus_taxid])[subgenus_taxid]))
        rcsb_id_taxid_tuples =  descendants_of_taxid( pdbid_taxid_tuples, subgenus_taxid )

        if not os.path.exists('/home/rtviii/dev/riboxyz/by_spec/{}'.format(subgenus_taxid)):
            os.mkdir('/home/rtviii/dev/riboxyz/by_spec/{}'.format(subgenus_taxid))
        for i,(rcsb_id, taxid) in enumerate( rcsb_id_taxid_tuples ):
                try:
                    shutil.copyfile('/home/rtviii/dev/RIBETL_DATA/{}/{}_PTC_COORDINATES.json'.format(rcsb_id, rcsb_id), '/home/rtviii/dev/riboxyz/by_spec/{}/{}_PTC_COORDINATES.json'.format(subgenus_taxid, rcsb_id))
                except Exception as e:
                    print(e)

        # for i,(rcsb_id, taxid) in enumerate( rcsb_id_taxid_tuples ):
        #     __structs.append(rcsb_id)
        #     print("{}.".format(i),rcsb_id, taxid, ncbi.get_taxid_translator([taxid])[taxid])
           
elif sys.argv[1] == "classify_rna":
    prof = RibosomeAssets('3J7Z').profile()
    # p= prof.rnas
    [ rna23s, rna5s ]=prof.rnas
    seq = rna23s.entity_poly_seq_one_letter_code_can




    hmms = []
    alphabet    = pyhmmer.easel.Alphabet.rna()

    for val in RNAClassEnum:
        hmm_path = "class_{}_taxid_{}.hmm".format(val.value, prof.src_organism_ids[0])
        if os.path.isfile(os.path.join(hmm_cachedir, hmm_path)):
            hmm_path = os.path.join(hmm_cachedir, hmm_path)
            with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
                hmm = hmm_file.read()   
        else:
            hmm = hmm_produce(val, prof.src_organism_ids[0])
        hmms.append(hmm)
        
    
    # k = hmm_dict_init__candidates_per_organism(RNAClassEnum, prof.src_organism_ids[0])

    seq_  = pyhmmer.easel.TextSequence(name=b"template", sequence=seq)
    dsb   = DigitalSequenceBlock(alphabet, [seq_.digitize(alphabet)])
    for hmm in hmms:
        print("\t>>>>>>>>>>>> ", hmm)
        hits =  pyhmmer.plan7.Pipeline(alphabet=alphabet).search_hmm(hmm,dsb)
        for hit in hits:
            print(hit.evalue)
