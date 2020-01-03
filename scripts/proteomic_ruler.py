import numpy as np
import pandas as pd
from pathlib import Path
import preprocessing
import re
import sys
import argparse
from statistics import mean
from Bio import SeqIO

parent = Path(__file__).resolve().parent
out_dir = parent / 'data'
parser = argparse.ArgumentParser(description='Apply proteomic ruler to MaxQuant output.')
parser.add_argument('txt_folder', type=str, help='MaxQuant txt folder')
parser.add_argument('out_file', type=str, help='Output filename (csv)')
#parser.add_argument('--out_dir', default=out_dir, type=str, help='Output directory')
parser.add_argument('--fasta', default='Reference_proteome/uniprot-proteome-UP000002494.fasta',
                    type=str, help='FASTA for reference proteome')
parser.add_argument('--histone', default='data/rat_histones_091918.csv',
                    type=str, help='Histones UniProt accessions')
parser.add_argument('--email', default='name@domain.com',
                    type=str, help='Email address for UniProt programmatically access.')

args = parser.parse_args()


ref_proteome = SeqIO.parse(args.fasta, "fasta")
proteome = list(ref_proteome)

digested_proteome = list(map(lambda y: preprocessing.digest(y, misscleavage=0),list(map(lambda x: str(x.seq), proteome))))
theoretical_peptide_py = dict(zip(list(map(lambda x: x.name.split("|")[1], proteome)),
    list(map(lambda x: len([pep for pep in x if 7 <= len(pep) <= 30]), digested_proteome))))


txt = Path(args.txt_folder)

#import Histone UniProtAcc for Rat
histone_rat = pd.read_csv(args.histone)
histone_rat = histone_rat['Entry'].tolist()

#import data

proteinGroups = pd.read_csv(txt/"proteinGroups.txt",sep = '\t', low_memory = False)
#peptides = pd.read_csv(txt/"peptides.txt",sep = '\t', low_memory = False)
summary = pd.read_csv(txt/"summary.txt",sep = '\t', low_memory = False)
#evidence = pd.read_csv(txt/"evidence.txt",sep = '\t', low_memory = False)

experiment_name = summary['Experiment'].dropna().unique()
summary = summary[summary['Raw file'] != 'Total']

# Remove 'Only identified by site', 'Reverse', or 'Potential contaminant' from proteinGroups.

proteinGroups_cleaned = preprocessing.cleanup_pg(proteinGroups)
preprocessing.convert_mw(proteinGroups_cleaned)
proteinGroups_cleaned['Number of theoretical peptides (trypsin/P, 7-30)'] = proteinGroups_cleaned['Majority protein IDs'].map(lambda x: mean(list(map(theoretical_peptide_py.get, x.split(";")))))
preprocessing.add_detectability(proteinGroups_cleaned, 'Number of theoretical peptides (trypsin/P, 7-30)')

# Update gene names
proteinGroups_cleaned = preprocessing.update_genes(proteinGroups_cleaned, args.email)
intensity_cols = [col for col in proteinGroups_cleaned.columns if 'Intensity ' in col]
proteinGroups_cleaned = preprocessing.check_histone(proteinGroups_cleaned, histone_rat)
proteinGroups_cleaned_intensity = preprocessing.copy_number(proteinGroups_cleaned, intensity_cols)
preprocessing.export_csv(proteinGroups_cleaned_intensity, Path(args.out_file).name, Path(args.out_file).parent, experiment_name)
