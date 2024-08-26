import sys
import os
from os import path
from zipfile import ZipFile

from cobrakbase.core.kbasefba import FBAModel
# from modelseedpy import AnnotationOntology, MSPackageManager, MSMedia, MSModelUtil, MSBuilder, MSATPCorrection, MSGapfill, MSGrowthPhenotype, MSGrowthPhenotypes, ModelSEEDBiochem
# from modelseedpy.core.annotationontology import convert_to_search_role, split_role
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.core.msgenome import normalize_role
# from modelseedpy.core.msensemble import MSEnsemble
from modelseedpy.community.mscommunity import MSCommunity
from modelseedpy.helpers import get_template

# Add the parent directory to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from baseutil import *
import hashlib
import pandas as pd
from pandas import DataFrame, read_csv, concat, set_option
from cobra.io import write_sbml_model, read_sbml_model



class CliffCommUtil(BaseUtil):
    def __init__(self):
        BaseUtil.__init__(self,"Cliff")
        # self.msseedrecon()
    
    def load_function_data(self,small=True):
        if small:
            return json.load(open('data/annotation_ani_prob_gep_85.json'))
        else:
            return json.load(open('data/annotation_ani_prob_lo_70.json'))
        
    def create_phenotypeset_from_compounds(
        self,
        compounds,
        base_media=None,
        base_uptake=0,
        base_excretion=1000,
        global_atom_limits={},
        type="growth"
    ):
        cpd_hash = {}
        for cpd in compounds:
            cpd_hash[cpd] = 10
        return MSGrowthPhenotypes.from_compound_hash(
            cpd_hash,
            base_media=base_media,
            base_uptake=base_uptake,
            base_excretion=base_excretion,
            global_atom_limits=global_atom_limits,
            type=type
        )

    def translate_protein_to_gene(self,protein):
        back_translation_code = {
            'A': ['GCA', 'GCC', 'GCG', 'GCT'],
            'C': ['TGT', 'TGC'],
            'D': ['GAC', 'GAT'],
            'E': ['GAG', 'GAA'],
            'F': ['TTT', 'TTC'],
            'G': ['GGT', 'GGG', 'GGA', 'GGC'],
            'H': ['CAT', 'CAC'],
            'I': ['ATC', 'ATA', 'ATT'],
            'K': ['AAG', 'AAA'],
            'L': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'],
            'M': ['ATG'],
            'N': ['AAC', 'AAT'],
            'P': ['CCT', 'CCG', 'CCA', 'CCC'],
            'Q': ['CAA', 'CAG'],
            'R': ['AGG', 'AGA', 'CGA', 'CGC', 'CGG', 'CGT'],
            'S': ['AGC', 'AGT', 'TCT', 'TCG', 'TCC', 'TCA'],
            'T': ['ACA', 'ACG', 'ACT', 'ACC'],
            'V': ['GTA', 'GTC', 'GTG', 'GTT'],
            'W': ['TGG'],
            'Y': ['TAT', 'TAC'],
            '*': ['TAA', 'TGA', 'TAG']
        }
        dna = ""
        for aa in protein:
            if aa in back_translation_code:
                dna += back_translation_code[aa][0]
        dna += 'TAA'
        return dna

util = CliffCommUtil() 