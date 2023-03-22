import argparse
from ctypes import alignment
import json
import re
import os
from typing import List, Union
import sys
import typing
from Bio import pairwise2
import itertools

from pydantic import BaseModel
from ribctl.lib.types.types_ribosome import PolymerClass, RibosomeStructure
from ribctl.lib.mod_extract_bsites import  BindingSite, struct_ligand_ids, struct_liglike_ids
from ribctl.lib.utils import open_structure
import numpy as np

flatten = itertools.chain.from_iterable
n1      = np.array


class SeqMatch(BaseModel):

	def __init__(self,
		sourceseq:str,
		targetseq:str, 
		source_residues:list[int]):
		"""A container for origin and target sequences when matching the resiudes of a ligand binding site
		to another protein's sequence through BioSeq's Align
		 """

		#* Computed indices of the ligand-facing in the source sequence.
		self.src     :str      = sourceseq
		self.src_ids:list[int] = source_residues

		#* Indices of predicted residues in target sequence. To be filled.
		self.tgt     :str      = targetseq
		self.tgt_ids:list[int] = []
		
		_            = pairwise2.align.globalxx(self.src,self.tgt, one_alignment_only=True)
		self.src_aln = _[0].seqA
		self.tgt_aln = _[0].seqB

		self.aligned_ids = []

		for src_resid in self.src_ids:
			self.aligned_ids.append(self.forwards_match(self.src_aln,src_resid))

		self.aligned_ids = list(filter(lambda x: x != None, self.aligned_ids ))

		for aln_resid in self.aligned_ids:
			if self.tgt_aln[aln_resid] == '-':
				continue
			self.tgt_ids.append(self.backwards_match(self.tgt_aln,aln_resid))

	def backwards_match(self, alntgt:str, resid:int):
		"""Returns the target-sequence index of a residue in the (aligned) target sequence"""
		if resid > len(alntgt):
			exit(IndexError(f"Passed residue with invalid index ({resid}) to back-match to target.Seqlen:{len(alntgt)}"))
		counter_proper = 0
		for i,char in enumerate(alntgt):
			if i == resid:
				return counter_proper
			if char =='-':
				continue
			else: 
				counter_proper  +=1

	def forwards_match(self,alnsrc:str, resid:int):
		"""Returns the index of a source-sequence residue in the aligned source sequence."""
		count_proper = 0
		for alignment_indx,char in enumerate( alnsrc ):
			if count_proper == resid:
				return alignment_indx
			if char =='-':
				continue
			else: 
				count_proper  +=1

	@staticmethod
	def hl_subseq(sequence:str, subsequence:str, index:int=None):
		"""Highlight subsequence"""
		CRED = '\033[91m'
		CEND = '\033[0m'
		_ = [ ]
		if index != None:
			return sequence[:index-1] + CRED + sequence[index] + CEND +sequence[index+1:]
		for item in re.split(re.compile(f'({subsequence})'),sequence):
			if item == subsequence:
				_.append(CRED + item + CEND)
			else:
				_.append(item)
		return ''.join(_)

	@staticmethod
	def hl_ixs(sequence:str,  ixs:List[int]):
		"""Highlight indices"""
		CRED = '\033[91m'
		CEND = '\033[0m'
		_ = ''
		for i,v in enumerate(sequence):
			if i in ixs: _ += CRED + v +CEND
			else: 	 	 _ += v
		return _


def struct_bsites(rcsb_id:str):
	"""Returns a list of binding sites from a structure"""
	rcsb_id =rcsb_id.upper()
	struct_profile_handle = RibosomeStructure.parse_obj(open_structure(rcsb_id, 'json'))
	liglike_polys = struct_liglike_ids(struct_profile_handle)
	ligands       = struct_ligand_ids(rcsb_id, struct_profile_handle)
	return ligands, liglike_polys

def open_bsite(path:str)->BindingSite:
	with open(os.path.join(path), 'rb') as infile:
		data = json.load(infile)
	return BindingSite(nbr_chains=data)

class PredictedResiduesPolymer(BaseModel):
	class PredictionSource(BaseModel):
		src: str
		src_ids: list[int]
		auth_asym_id: str
	class PredictionTarget(BaseModel):
		tgt: str
		tgt_ids: list[int]
		auth_asym_id: str
	class PredictionAlignments(BaseModel):
		aln_ids: list[int]
		src_aln: str
		tgt_aln: str

	source   : PredictionSource
	target   : PredictionTarget
	alignment: PredictionAlignments
	
def init_transpose_ligand(
	# source_struct: str,
	# target_struct: str,
	target_profile:RibosomeStructure,
	binding_site : BindingSite
	)->dict[PolymerClass, PredictedResiduesPolymer]: 

	by_class_origin_polymers:dict[PolymerClass, dict] = {}


	origin_polymers = binding_site.nbr_chains
	target_polymers = itertools.chain(target_profile.rnas if target_profile.rnas is not None else [],target_profile.proteins)

	#? For every chain in a bindingsite, if it has nomenclature, append its residues, strand and sequence
	for auth_asym_id, polymer in origin_polymers.items():
		if len(polymer.nomenclature) <1:
			continue
		else:
			by_class_origin_polymers[polymer.nomenclature[0]] = {
				'seq'         : polymer.entity_poly_seq_one_letter_code_can,
				'auth_asym_id': polymer.auth_asym_id,
				'ids'         : [ resid for resid in [*map(lambda x : x.seqid, polymer.residues)] ]
			}

	
	
	by_class_target_polymers:dict[PolymerClass, dict] = {}

	for nomenclature_class, polymer in by_class_origin_polymers.items():
		# find a polymer in target structure with same nomenclature class
		matches =  [*filter(lambda tgt_poly: nomenclature_class in tgt_poly.nomenclature, target_polymers)]
		if len( matches )  < 1:
			continue

		tgt_poly_seq          = matches[0].entity_poly_seq_one_letter_code_can
		tgt_poly_auth_asym_id = matches[0].auth_asym_id

		by_class_target_polymers[nomenclature_class] ={
			'seq'         : tgt_poly_seq,
			'auth_asym_id': tgt_poly_auth_asym_id,
		}


	for nomenclature_class in by_class_origin_polymers:
		if nomenclature_class not in by_class_target_polymers:
			continue

		src_ids = by_class_origin_polymers[nomenclature_class]['ids']
		src     = by_class_origin_polymers[nomenclature_class]['seq']
		tgt     = by_class_target_polymers[nomenclature_class]['seq']

		sq      = SeqMatch(src,tgt,src_ids)

		src_aln = sq.src_aln     # <--- aligned source      sequence (with                        gaps)
		tgt_aln = sq.tgt_aln     # <--- aligned tgt         sequence (with                        gaps)
		aln_ids = sq.aligned_ids # <--- ids     corrected   for                                   gaps
		tgt_ids = sq.tgt_ids     # <--- ids     backtracted to the target polymer (accounting for gaps)



		prediction[nomenclature_class] = PredictedResiduesPolymer.parse_obj({
			"source":{
				"src"         : src,
				"src_ids"     : src_ids,
				"auth_asym_id": by_class_origin_polymers[nomenclature_class]['auth_asym_id']
			},
			"target":{
				"tgt"         : tgt,
				"tgt_ids"     : tgt_ids,
				'auth_asym_id': by_class_target_polymers[nomenclature_class]['auth_asym_id']
			},
			"alignment" :{
				"aln_ids": aln_ids,
				"src_aln": src_aln,
				"tgt_aln": tgt_aln,
			},
		})

	return prediction


if __name__ =="__main__":

	from rbxz_bend.settings import RIBETL_DATA

	def root_self(rootname: str = ''):
		"""Returns the rootpath for the project if it's unique in the current folder tree."""
		root = os.path.abspath(__file__)[:os.path.abspath(
		__file__).find(rootname)+len(rootname)]
		sys.path.append(root)

	root_self('ribetl')

	prs = argparse.ArgumentParser()

	prs.add_argument('-src', '--source_structure', type=str, required=True)
	prs.add_argument('-tgt', '--target_structure', type=str, required=True)
	prs.add_argument('-p', '--path', type=str, required=True)

	args = prs.parse_args()

	SRC_STRUCT = args.source_structure.upper()
	TGT_STRUCT = args.target_structure.upper()
	lig_path = str( args.path )


	target_profile = RibosomeStructure.parse_obj(open_structure(TGT_STRUCT,'json'))

	bsite      = open_bsite(lig_path)
	_type:typing.Literal["LIGAND","POLYMER"] = "LIGAND" if "ligand" in lig_path.lower()  else "POLYMER"

	prediction = init_transpose_ligand(target_profile,bsite)
	fname = f'PREDICTION_{_type}_{SRC_STRUCT}_{TGT_STRUCT}.json'
	with open(os.path.join(RIBETL_DATA,TGT_STRUCT,fname), 'w') as outfile:
		json.dump(prediction,outfile)
		print("\033[093mSucessfully saved prediction {}\033[0m".format(fname))