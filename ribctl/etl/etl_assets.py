from enum import   auto
import enum
import json
import os
import typing
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from loguru import logger
from ribctl import AMINO_ACIDS_3_TO_1_CODE, CLASSIFICATION_REPORTS
from ribctl.lib.libtax import PhylogenyNode, PhylogenyRank, Taxid
from ribctl.lib.tunnel import ptc_resdiues_get, ptc_residues_calculate_midpoint
from ribctl.lib.utils import download_unpack_place, open_structure
from ribctl.lib.schema.types_ribosome import ( RNA, PTCInfo, Polymer, PolymerClass, PolynucleotideClass, PolypeptideClass, RibosomeStructure, )
from ribctl import RIBETL_DATA
from ribctl.logs.loggers import get_etl_logger

