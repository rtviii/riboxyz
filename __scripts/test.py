import asyncio
import pickle
from pyhmmer.easel import (
    Alphabet,
    DigitalSequenceBlock,
    TextSequence,
    SequenceFile,
    SequenceBlock,
    TextSequenceBlock,
)
from pyhmmer.plan7 import Pipeline, HMM
from functools import reduce
from io import StringIO
from itertools import tee
import json
import os
from pprint import pprint
from tempfile import NamedTemporaryFile
from typing import Iterator
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment, Seq, SeqRecord
from pyhmmer import phmmer
import pyhmmer
from ribctl import ASSETS, LOGS_PATH, NCBI_TAXA_SQLITE, RIBETL_DATA
from ribctl.etl.etl_pipeline import (
    ReannotationPipeline,
    query_rcsb_api,
    rcsb_single_structure_graphql,
)
from ribctl.etl.etl_obtain import obtain_asssets_threadpool
from ribctl.lib.libhmm import (
    HMMClassifier,
    PolymerClassesOrganismScanner,
    hmm_create,
    hmm_produce,
)
from ribctl.etl.etl_assets_ops import Assetlist, RibosomeOps, Structure
from ribctl.lib.schema.types_ribosome import (
    LifecycleFactorClass,
    Polymer,
    PolynucleotideClass,
)
from ribctl import model_subgenuses

from ribctl.lib.taxlib import (
    descendants_of_taxid,
    taxid_superkingdom,
    taxid_is_descendant_of,
)
from ete3 import NCBITaxa

hmm_cachedir = ASSETS["__hmm_cache"]

import sys


if sys.argv[1] == "tunnel":

    def list_euk_structs():
        EUK_STRUCTS = []
        with open("eukarya_2023.txt", "r") as data_file:
            for line in data_file:
                structs = line.split(",")
                EUK_STRUCTS = [*EUK_STRUCTS, *structs]
        return EUK_STRUCTS

    if __name__ == "__main__":
        EUK = list_euk_structs()
    print("tunnel")
elif sys.argv[1] == "spec":

    def get_taxonomic_id(organism_name):
        ncbi = NCBITaxa(dbfile=NCBI_TAXA_SQLITE)
        try:
            # Get the taxonomic ID for the given organism name
            taxid = ncbi.get_name_translator([organism_name])[organism_name]
            return taxid
        except KeyError:
            # Handle the case where the organism name is not found
            print(
                f"Organism '{organism_name}' not found in the NCBI Taxonomy database."
            )
            return None

    _ = {}
    for organism_name in [
        "Lactococcus lactis",
        "Mycobacterium smegmatis",
        "Mycobacterium tuberculosis",
        "Bacillus subtilis",
        "Leishmania donovani",
        "Trypanosoma cruzi",
        "Trichomonas vaginalis",
        "Giardia duodenalis",
        "Spraguea lophii",
    ]:
        taxonomic_id = get_taxonomic_id(organism_name)
        if taxonomic_id:
            _ = {**_, **{organism_name: taxonomic_id}}
            print(f"Taxonomic ID for {organism_name}: {taxonomic_id}")
    print(_)
elif sys.argv[1] == "test":
    ncbi = NCBITaxa()
    all_structs = os.listdir(RIBETL_DATA)
    pdbid_taxid_tuples: list = []
    __structs = []

    for i, struct in enumerate(all_structs):
        try:
            print(struct, i)
            rp = RibosomeOps(struct)
            rp = rp.profile()
            pdbid_taxid_tuples.append((rp.rcsb_id, rp.src_organism_ids[0]))
        except:
            ...
    for name, subgenus_taxid in model_subgenuses.items():
        print(
            "\n\t###########      Descendants of {} ({})      ##########3".format(
                subgenus_taxid,
                ncbi.get_taxid_translator([subgenus_taxid])[subgenus_taxid],
            )
        )
        rcsb_id_taxid_tuples = descendants_of_taxid(pdbid_taxid_tuples, subgenus_taxid)

        for i, (rcsb_id, taxid) in enumerate(rcsb_id_taxid_tuples):
            __structs.append(rcsb_id)
            print(
                "{}.".format(i),
                rcsb_id,
                taxid,
                ncbi.get_taxid_translator([taxid])[taxid],
            )

    print(__structs)

    # ids_in_this_branch = [*list(map(lambda x: x[1], rcsb_id_taxid_tuples)), subgenus_taxid]
    # tree =ncbi.get_topology(ids_in_this_branch)
    # print(tree.get_ascii(attributes=["taxid"]))

    # obtain_assets_threadpool([rcsb_id for (rcsb_id, taxid) in rcsb_id_taxid_tuples], Assetlist(ptc_coords=True), overwrite=True)
elif sys.argv[1] == "processtax":
    RCSB_ID = "5MYJ"
    asyncio.run(
        ReannotationPipeline(
            query_rcsb_api(rcsb_single_structure_graphql(RCSB_ID))
        ).process_structure()
    )
elif sys.argv[1] == "ll":
    import shutil

    structs = [
        "5MYJ",
        "6DZI",
        "5XYM",
        "5ZEB",
        "7S0S",
        "5O61",
        "5ZET",
        "7Y41",
        "5ZEP",
        "6DZK",
        "5ZEU",
        "5O60",
        "5O5J",
        "6DZP",
        "7XAM",
        "5XYU",
        "7MT2",
        "5V7Q",
        "7MSM",
        "7MT7",
        "5V93",
        "7MSC",
        "7KGB",
        "7MSH",
        "7F0D",
        "7MSZ",
        "7MT3",
        "7SFR",
        "7AQC",
        "8BUU",
        "7SAE",
        "7S9U",
        "7QV2",
        "6PPK",
        "7O5B",
        "7QGU",
        "3J3W",
        "6HA1",
        "6HTQ",
        "6HA8",
        "7OPE",
        "3J3V",
        "3J9W",
        "7QV1",
        "7QV3",
        "7QH4",
        "6PPF",
        "7AQD",
        "6PVK",
        "6AZ3",
        "7ANE",
        "5T2A",
        "3JCS",
        "6AZ1",
        "7AIH",
        "7AOR",
        "5OPT",
        "5T5H",
        "5XY3",
        "5XYI",
        "8BTR",
        "7PWO",
        "8BR8",
        "7PWG",
        "8BSJ",
        "8BSI",
        "7PWF",
        "8BTD",
        "8BRM",
        "7QCA",
        "8P60",
        "8P5D",
    ]

    ncbi = NCBITaxa()
    all_structs = os.listdir(RIBETL_DATA)
    pdbid_taxid_tuples: list = []
    __structs = []

    for i, struct in enumerate(all_structs):
        try:
            print(struct, i)
            rp = Structure(struct)
            rp = rp.profile()
            pdbid_taxid_tuples.append((rp.rcsb_id, rp.src_organism_ids[0]))
        except:
            ...

    for name, subgenus_taxid in model_subgenuses.items():
        print(
            "\n\t###########      Descendants of {} ({})      ##########3".format(
                subgenus_taxid,
                ncbi.get_taxid_translator([subgenus_taxid])[subgenus_taxid],
            )
        )
        rcsb_id_taxid_tuples = descendants_of_taxid(pdbid_taxid_tuples, subgenus_taxid)

        if not os.path.exists(
            "/home/rtviii/dev/riboxyz/by_spec/{}".format(subgenus_taxid)
        ):
            os.mkdir("/home/rtviii/dev/riboxyz/by_spec/{}".format(subgenus_taxid))
        for i, (rcsb_id, taxid) in enumerate(rcsb_id_taxid_tuples):
            try:
                shutil.copyfile(
                    "/home/rtviii/dev/RIBETL_DATA/{}/{}_PTC_COORDINATES.json".format(
                        rcsb_id, rcsb_id
                    ),
                    "/home/rtviii/dev/riboxyz/by_spec/{}/{}_PTC_COORDINATES.json".format(
                        subgenus_taxid, rcsb_id
                    ),
                )
            except Exception as e:
                print(e)

        # for i,(rcsb_id, taxid) in enumerate( rcsb_id_taxid_tuples ):
        #     __structs.append(rcsb_id)
        #     print("{}.".format(i),rcsb_id, taxid, ncbi.get_taxid_translator([taxid])[taxid])

        print("\t>>>>>>>>>>>> ", hmm)
        hits = pyhmmer.plan7.Pipeline(alphabet=alphabet).search_hmm(hmm, dsb)
        for hit in hits:
            print(hit.evalue)
            print(hit.score)
            print(hit.bias)
            print(hit.description)
            print(hit.sum_score)
elif sys.argv[1] == "tsv_to_fasta":
    import csv

    mitoproteins_path = (
        "/home/rtviii/dev/riboxyz/ribctl/assets/fasta_proteins_mitochondrial"
    )

    for i in os.listdir(mitoproteins_path):
        tsv_path = os.path.join(mitoproteins_path, i)
        fasta_path = os.path.join(mitoproteins_path, "{}.fasta".format(i.split(".")[0]))
        records = []

        with open(tsv_path, "r", newline="", encoding="utf-8") as tsvfile:
            tsvreader = csv.reader(tsvfile, delimiter="\t")
            next(tsvreader, None)
            for row in tsvreader:
                entry, orgid, sequence = row
                records.append((entry, orgid, sequence))

        seqrecords = []

        for record in records:
            (entry, taxid, sequence) = record
            seqrecords.append(
                SeqRecord(
                    Seq(sequence), id=taxid, description="uniprot_{}".format(entry)
                )
            )

        dest = os.path.join(fasta_path)
        with open(dest, "w") as output_handle:
            SeqIO.write(seqrecords, output_handle, "fasta")
            print("Wrote {} seqs to  to {}".format(len(seqrecords), dest))
elif sys.argv[1] == "struct_factors":
    for struct in Structure.list_all_structs()[:10]:
        print(
            "========================Processing {}=====================".format(struct)
        )
        prof = Structure(struct).profile()
        prots = (
            prof.proteins + prof.polymeric_factors
            if prof.polymeric_factors != None
            else prof.proteins
        )
        candidate_category = LifecycleFactorClass
        hmm_organisms_registry = {}

        for chain in prots:
            chain_organism_taxid = chain.src_organism_ids[0]
            if chain_organism_taxid not in [*hmm_organisms_registry.keys()]:
                hmm_organisms_registry[chain_organism_taxid] = (
                    hmm_dict_init__candidates_per_organism(
                        candidate_category, chain_organism_taxid
                    )
                )

            k = classify_subchain(chain, hmm_organisms_registry[chain_organism_taxid])
            if k[1] != None:
                print(
                    "Classifying chain {}.{}".format(
                        chain.parent_rcsb_id, chain.auth_asym_id
                    )
                )
                pprint(k)
elif sys.argv[1] == "collect_factors":
    factor_structs = [
        "8FL4",
        "8FL6",
        "8FL7",
        "8FL9",
        "8FLA",
        "8FLB",
        "8FLC",
        "8FLD",
        "8FLE",
        "8FLF",
        "8G5Z",
        "8G6J",
        "8HFR",
        "8HKU",
        "8HL2",
        "8HL3",
        "8I9R",
        "8I9T",
        "8I9V",
        "8I9W",
        "8I9X",
        "8I9Y",
        "8I9Z",
        "8IA0",
        "8IDT",
        "8IDY",
        "8IE3",
        "8INE",
        "8INF",
        "8INK",
        "8IPD",
        "8IPX",
        "8IPY",
        "8IR1",
        "8IR3",
        "8OIN",
        "8OIP",
        "8OIQ",
        "8OIR",
        "8OIS",
        "8OIT",
    ]
    factors = {}

    for struct in factor_structs:
        try:
            prof = Structure(struct).profile()
            for p in prof.polymeric_factors:
                if "Factor" in p.nomenclature[0]:
                    factors[struct] = json.loads(p.json())
        except:
            ...

    with open("factors_sample.json", "w") as outfile:
        json.dump(factors, outfile)
elif sys.argv[1] == "hmmt":

    for rcsb_id in Structure.list_all_structs():
        prof = Structure(rcsb_id).profile()

        proteins = [*prof.proteins, *prof.other_polymers, *prof.polymeric_factors]
        rna = [*prof.rnas, *prof.other_polymers]

        pipeline_polypeptides = HMMClassifier(proteins, Alphabet.amino())
        pipeline_polypeptides.___scan_chains()
        prots_report = pipeline_polypeptides.produce_classification()

        pipeline_polynucleotides = HMMClassifier(proteins, Alphabet.rna())
        pipeline_polynucleotides.___scan_chains()
        rna_report = pipeline_polynucleotides.produce_classification()

        report_path = os.path.join(
            LOGS_PATH,
            "classification_reports",
            "{}_classification_report.json".format(rcsb_id),
        )
        pipeline.write_classification_report(report_path)
elif sys.argv[1] == "hmmx":
    p = Structure("3j7z").profile()
    rnas = p.rnas
    alphabet = pyhmmer.easel.Alphabet.rna()
    alphabet = pyhmmer.easel.Alphabet.amino()
    # p.other_polymers
    cca = p.polymeric_factors[1]
    # cca= p.rnas[0]
    # pprint(cca)
    # for rna in [ * p.polymeric_factors]:
    #     pprint(rna)

    hmms = PolymerClassesOrganismScanner(cca["src_organism_ids"][0], [p for p in RNAClass], no_cache=True)
    hmms.scan_seq(alphabet, [SeqRecord(cca["entity_poly_seq_one_letter_code_can"])])
    # pprint(hmms.info())

    # pprint()

    # pipeline_polynucleotides = HMMClassifier( rna, Alphabet.rna())
    # pipeline_polynucleotides.scan_chains()
    # rna_report = pipeline_polynucleotides.produce_classification()

    # report_path = os.path.join(LOGS_PATH,'classification_reports','{}_classification_report.json'.format(rcsb_id))
    # pipeline.write_classification_report(report_path)

# +
#     hmmx        = HMMScanner(
#         prof.src_organism_ids[0],
#         [pc for pc in [*list(ProteinClass), *list(LifecycleFactorClass)]],
#         name         = "max_seed_seq_5",
#         no_cache     = True,
#         max_seed_seq = 5
#     )

#     prots          :list[Polymer] = prof.proteins
#     factor_polymers:list[Polymer] = prof.polymeric_factors
#     other_polymers :list[Polymer] = prof.other_polymers

#     alphabet     = pyhmmer.easel.Alphabet.amino()
#     report_state = {}

#     for chain in  [ *prots, *factor_polymers, *other_polymers ]:

#         report_state[chain.auth_asym_id] = []
#         seq_record = chain.to_SeqRecord()
#         query_seq  = pyhmmer.easel.TextSequence(name=bytes(seq_record.id,'utf-8'), sequence=seq_record.seq)
#         query_seqs = [query_seq.digitize(alphabet)]

#         for scan in list(pyhmmer.hmmscan(query_seqs,[*hmmx.hmms_registry.values()])):

#             for hit in [*scan]:
#                hit_dict = {
#                 "hit.name:"  : hit.name.decode(),
#                 "hit.evalue:": hit.evalue,
#                 "hit.score:" : hit.score,
#                 }

#                report_state[chain.auth_asym_id].append(hit_dict)

#     with open("{}.json".format(report_name), "w") as outfile:
#         json.dump(report_state, outfile)


# # TODO: deliberate on each chain's taxid
# # TODO: pick best hit
# # TODO: save report
