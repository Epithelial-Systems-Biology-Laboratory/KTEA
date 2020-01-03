import numpy as np
import pandas as pd
import urllib
from io import StringIO
from datetime import datetime
from pathlib import Path
from Bio import SeqIO
import re

# Constants
AVOGADRO = 6.022140857e23
BASEPAIR_WEIGHT = 615.8771
RAT_GENOME_SIZE = 2870184193 # Rnor_6.0 https://rgd.mcw.edu/rgdweb/report/genomeInformation/genomeInformation.html
PLOIDY = 2
C_VALUE = RAT_GENOME_SIZE*BASEPAIR_WEIGHT/AVOGADRO
TOTAL_PROTEIN_PER_CELL = 200 # pg


# Functions for preprocessing MaxQuant results

def cleanup_pg(df ,site = True, rev = True, con = True, iRT = True):
    '''
    Remove 'Only identified by site', 'Reverse', or 'Potential contaminant' from proteinGroups.
    '''
    tmp_df = df.copy()

    if site:
        tmp_df = tmp_df[tmp_df['Only identified by site'].isna()]
        tmp_df.drop('Only identified by site', axis = 1, inplace = True)
    if rev:
        tmp_df = tmp_df[tmp_df['Reverse'].isna()]
        tmp_df.drop('Reverse', axis = 1, inplace = True)
    if con:
        tmp_df = tmp_df[tmp_df['Potential contaminant'].isna()]
        tmp_df.drop('Potential contaminant', axis = 1, inplace = True)
    if iRT:
        tmp_df = tmp_df[~tmp_df['Protein IDs'].str.contains('iRT')]
    return tmp_df


### TODO: options for detectability normalizing factor

def convert_mw(df):
    """Convert MW from kDa to Da (*1000)"""
    df['Mol. weight [Da]'] = df['Mol. weight [kDa]'] * 1000

def add_detectability(df, col=None):
    """Add detectability factor by specified a column in df (e.g., number of theoretical peptides). If no column is given or not found in df, assign 'Mol. weight [Da]' instead. (REQUIRED conversion of MW from kDa to Da first!!)"""

    if col == None:
        df['detectability_factor'] = df['Mol. weight [Da]']
    else:
        if col in df.columns:
            df['detectability_factor'] = df[col]
            # No zero value allowed, replace '0' with '1'
            df[df['detectability_factor']==0] = 1
        else:
            print('!!No specified detectability column found!!')
            df['detectability_factor'] = df['Mol. weight [Da]']

def check_histone(df, histone):
    """Add a new column indicating whether that protein match to histone UniProt Acc. [Return as a new DataFrame]"""
    tmp = df.copy()
    tmp['Histone'] = False
    for proteins in tmp['Protein IDs'].index:
        if any(protein in histone for protein in tmp.loc[proteins,'Protein IDs'].split(';')):
            tmp.loc[proteins,'Histone'] = True
        else:
            tmp.loc[proteins,'Histone'] = False
    return tmp


def copy_number(df, expression_cols):
    """Calculate copy number of each protein."""
    histone_index = df[df['Histone'] == True].index
    tmp = df.copy()

    for sample in expression_cols:
        sum_mw_histone = 'sum_mw_histone_' + sample
        factor_histone = 'factor_histone_' + sample
        sum_mw = 'sum_mw_' + sample

        copy_num_histone = 'copy_number_histone_' + sample

        tmp[sum_mw_histone] = (tmp['Mol. weight [Da]'] / tmp['detectability_factor'] * tmp[sample])[histone_index].sum()

        # C_VALUE * PLOIDY is amount of DNA per cell [in gram], which assume to be equal to histone mass
        tmp[factor_histone] = C_VALUE * PLOIDY * AVOGADRO / tmp[sum_mw_histone]
        tmp[copy_num_histone] = tmp[sample] * tmp[factor_histone] / tmp['detectability_factor']
    return tmp


def get_genename(acc, email, cols='id,reviewed,genes(PREFERRED),genes'):

    url = 'https://www.uniprot.org/uploadlists/'
    params = {
    'from':'ACC',
    'to':'ACC',
    'format':'tab',
    'columns':cols,
    'query':acc
    }
    current_time = datetime.now()

    data = urllib.parse.urlencode(params).encode('UTF-8')
    request = urllib.request.Request(url, data)
    contact = email # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urllib.request.urlopen(request)
    page = response.read().decode('UTF-8')

    df = pd.read_csv(StringIO(page), sep = '\t')
    if df['Gene names'].isna().all():
        print('No gene name found!!')
        return None
    else:
        print(f'Gene names downloaded from UniProt at: {current_time}')
        if cols == 'id,reviewed,genes(PREFERRED),genes':
            return df[['Entry','Status', 'Gene names  (primary )', 'Gene names']].fillna('')
        else:
            return df.fillna('')


def id_to_gene(accs):
    try:
        return ';'.join(list(dict.fromkeys(' '.join([gene_dict[acc] for acc in accs.split(';')]).split())))
    except(KeyError):
        return ''


def update_genes(data, email):

    """A function to convert UniProt accession (from data['Majority protein IDs']) to gene name. For each entry in proteinGroup, the "preferred gene name" is chosen from the UniProt entry which is in SwissProt database (first entry will be used if multiple SwissProt entries available). If no "preferred gene name" from SwissProt entry found, the first "preferred gene name" is used instead. If no "preferred gene name" available, 'Gene names' reported from MaxQuant will be used.

    # Parameters
    - data: A dataframe imported from MaxQuant's "proteinGroups.txt" output file.
    - email: User email address for requesting data from UniProt.

    # Return
    - A new dataframe copied from data with additional columns:
        'reviewed'          :  Specify if any accession is in SwissProt (logical)
        'preferred_gene'    :  SwissProt preferred gene name
        'sec_pref_gene'     :  First preferred gene name in each row
        'final_gene'        :  Selected gene name
    """
    uniprotACC = ';'.join(data['Majority protein IDs'])
    tmp_df = get_genename(uniprotACC, email)
    gene_dict = dict(zip(tmp_df['Entry'],tmp_df['Gene names  (primary )']))
    sp = tmp_df.loc[tmp_df['Status'] == 'reviewed', 'Entry'].tolist()
    df = data.copy()
    df['reviewed'] = df['Majority protein IDs'].str.split(';', expand = True).isin(sp).any(axis = 1)
    df['first_sp_protein'] = df['Majority protein IDs'].str.split(';', expand = True).isin(sp).idxmax(axis = 1)
    preferred_gene = df['Majority protein IDs'].str.split(';', expand = True).applymap(gene_dict.get)
    df['preferred_gene'] = preferred_gene.lookup(df.index,df['first_sp_protein'])
    df['sec_pref_gene'] = (preferred_gene.fillna('') + ';').sum(axis = 1).str.strip(';').str.split(';').map(lambda x : x[0])
    df['final_gene'] = df['preferred_gene']
    df.loc[df['preferred_gene'] == "", 'final_gene'] = df.loc[df['preferred_gene'] == "", 'sec_pref_gene']
    df.loc[df['final_gene'] == "", 'final_gene'] =  df.loc[df['final_gene'] == "", 'Gene names']
    return df.drop('first_sp_protein', axis=1)



def export_csv(df, filename, out, experiment_name, additional_cols=None):

    """
    Export processed data to csv.
    ...
    """

    #out = Path('output/')
    try:
        Path.mkdir(out)
        print('Output folder created.')
    except(FileExistsError):
        print('Output folder already exists.')
    out_file = Path.joinpath(out,filename)

    copy_number_cols = [col for col in df.columns if col.startswith('copy')]

    export_cols = ['Majority protein IDs',
                  'Protein names',
                  'Gene names',
                  'final_gene'] + copy_number_cols + (
                  'Razor + unique peptides ' + experiment_name).tolist() + (
                  'Sequence coverage ' + experiment_name + ' [%]').tolist() + [
                  'id',
                  'detectability_factor',
                  'Histone']
    if(additional_cols):
        export_cols = export_cols + additional_cols

    df[export_cols].to_csv(out_file, index=False)
    print(f'Data exported to: {out_file.absolute()}')
    return None


def digest(protein_seq, cut_pattern = "[KR]", misscleavage = 0):

    """
    In silico digestion of protein(s)
    """


    cut_sites = [-1]
    start = 0
    for site in re.finditer(cut_pattern, protein_seq):
        cut_sites.append(site.span()[0])
    cut_sites.append(len(protein_seq)-1)
    digested_peptides = [protein_seq[cut_sites[i]+1:cut_sites[i+1]+1] for i in range(len(cut_sites)-1)]

    if misscleavage == 0:
        return(digested_peptides)
    elif misscleavage == 1:
        digested_peptides_miss1 = list(map("".join, zip(digested_peptides[:-1],digested_peptides[1:])))
        return(digested_peptides + digested_peptides_miss1)
    elif misscleavage == 2:
        digested_peptides_miss1 = list(map("".join, zip(digested_peptides[:-1],digested_peptides[1:])))
        digested_peptides_miss2 = list(map("".join, zip(digested_peptides[:-2],digested_peptides[1:-1],digested_peptides[2:])))
        return(digested_peptides + digested_peptides_miss1 + digested_peptides_miss2)
    else:
        print("Warning!! only 0-2 misscleavage allowed. Return peptides as misscleavage=2")
        digested_peptides_miss1 = list(map("".join, zip(digested_peptides[:-1],digested_peptides[1:])))
        digested_peptides_miss2 = list(map("".join, zip(digested_peptides[:-2],digested_peptides[1:-1],digested_peptides[2:])))
        return(digested_peptides + digested_peptides_miss1 + digested_peptides_miss2)


    
