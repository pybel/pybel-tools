# -*- coding: utf-8 -*-

"""
To run, type :code:`python3 -m pybel_tools.analysis.neurommsig.export` in the command line
"""

import itertools as itt
import logging

import os
import pandas as pd
import re
from functools import partial

from pybel.utils import ensure_quotes, get_bel_resource
from ...document_utils import write_boilerplate
from ...resources import DBSNP_PATTERN, NIFT, HGNC_HUMAN_GENES, NEUROMMSIG, MESHD

log = logging.getLogger(__name__)

HGNCsymbolpattern = re.compile("^[A-Z0-9-]+$|^C[0-9XY]+orf[0-9]+$")

SNPpattern = re.compile("^rs[0-9]+$")
SNPspatternSpace = re.compile("^(rs[0-9]+)\s((rs[0-9]+)\s)*(rs[0-9]+)$")
SNPspatternComma = re.compile("^(rs[0-9]+),((rs[0-9]+),)*(rs[0-9]+)$")
SNPspatternSpaceComma = re.compile("^(rs[0-9]+), ((rs[0-9]+), )*(rs[0-9]+)$")
Checked_by_Anandhi = re.compile("No")

miRNApattern = re.compile("^MIR.*$")
miRNAspattern = re.compile("^(MIR.*),((MIR.*$),)*(MIR.*$)$")


def preprocessing_excel(filepath):
    """Preprocessing of the excel file.

    Parameters
    ----------
    filepath : Filepath of the excel sheet

    Returns
    -------
    dataframe : Preprocessed pandas dataframe

    """

    if not os.path.exists(filepath):
        raise ValueError("Error: %s file not found" % filepath)

    # Import Models from Excel sheet, independent for AD and PD
    df = pd.read_excel(filepath, sheetname=0, header=0)

    # Indexes and column name
    # [log.info(str(x)+': '+str((df.columns.values[x]))) for x in range (0,len(df.columns.values))]

    # Starting from 4: Pathway Name

    # Fill Pathway cells that are merged and are 'NaN' after deleting rows where there is no genes
    df.iloc[:, 4] = pd.Series(df.iloc[:, 4]).fillna(method='ffill')

    # Number of gaps
    # log.info(df.ix[:,6].isnull().sum())

    df = df[df.ix[:, 6].notnull()]
    df = df.reset_index(drop=True)

    # Fill NaN to ceros in PubmedID column
    df.ix[:, 7].fillna(0, inplace=True)

    # Number of gaps in the gene column should be already zero
    if (df.ix[:, 6].isnull().sum()) != 0:
        raise ValueError("Error: Empty cells in the gene column")

    # Check current state
    # df.to_csv('out.csv')

    return df


def munge_cell(cell, line=None, validators=None):
    """

    :param cell:
    :param line:
    :param validators:
    :return:
    """
    if pd.isnull(cell):
        return None

    c = ' '.join(cell.split())

    if validators is not None and all(re.match(validator, c) is None for validator in validators):
        if line:
            log.info("Munge cell error: aprox in line: %s: %s", line, c)
        return None

    return [x.strip() for x in str(c).strip().split(',')]


def preprocessing_br_projection_excel(filepath):
    """Preprocessing of the excel file.

    Parameters
    ----------
    filepath : Filepath of the excel sheet

    Returns
    -------
    dataframe : Preprocessed pandas dataframe

    """

    if not os.path.exists(filepath):
        raise ValueError("Error: %s file not found" % filepath)

    df = pd.read_excel(filepath, sheetname=0, header=0)

    return df


munge_snp = partial(munge_cell, validators=[SNPpattern, SNPspatternSpaceComma])

mesh_alzheimer = "Alzheimer Disease"  # Death to the eponym!
mesh_parkinson = "Parkinson Disease"

pathway_column = 'Pathway Name (Daniel & Apurva)'
columns = [
    'Genes (Daniel & Apurva)',
    'SNPs from Literature (Aybuge)',
    'Genome wide associated SNPs (Mufassra)',
    'Imaging Features (Anandhi)',
    'SNP_Image Feature (Mufassra & Anandhi)',
]


def preprocess(path):
    """

    :param str path:
    :return:
    :rtype: pandas.DataFrame
    """
    df = preprocessing_excel(path)
    df['SNPs from Literature (Aybuge)'] = df['SNPs from Literature (Aybuge)'].map(munge_snp)
    df['Genome wide associated SNPs (Mufassra)'] = df['Genome wide associated SNPs (Mufassra)'].map(munge_snp)
    df['LD block analysis (Mufassra)'] = df['LD block analysis (Mufassra)'].map(munge_snp)
    df['Imaging Features (Anandhi)'] = df['Imaging Features (Anandhi)'].map(munge_cell)
    df['Imaging Features (Anandhi)'] = df['Imaging Features (Anandhi)'].map(
        lambda c: None if c is not None and c[0] == 'No' else c)
    df['SNP_Image Feature (Mufassra & Anandhi)'] = df['SNP_Image Feature (Mufassra & Anandhi)'].map(munge_snp)
    return df


def get_nift_values():
    """

    :return:
    :rtype: dict[str,str]
    """
    r = get_bel_resource(NIFT)
    return {v.lower(): v for v in r['Values']}


def write_neurommsig_bel(file, df, disease, nift_values):
    """Writes the NeuroMMSigDB excel sheet to BEL

    :param file: a file or file-like that can be writen to
    :param pandas.DataFrame df: 
    :param str disease: 
    :param dict nift_values: a dictionary of lowercased to normal names in NIFT
    """
    write_boilerplate(
        document_name='NeuroMMSigDB for {}'.format(disease),
        description='SNP and Clinical Features for Subgraphs in {}'.format(disease),
        authors='Daniel Domingo, Charles Tapley Hoyt, Mufassra Naz, Aybuge Altay, Anandhi Iyappan',
        contact='charles.hoyt@scai.fraunhofer.de',
        namespace_dict={
            'NIFT': NIFT,
            'HGNC': HGNC_HUMAN_GENES,
        },
        namespace_patterns={
            'dbSNP': DBSNP_PATTERN
        },
        annotations_dict={
            'Subgraph': NEUROMMSIG,
            'MeSHDisease': MESHD
        },
        file=file
    )

    print('SET Citation = {"URL", "NeuroMMSigDB", "http://neurommsig.scai.fraunhofer.de/"}', file=file)
    print('SET Evidence = "Serialized from NeuroMMSigDB"', file=file)
    print('SET MeSHDisease = "{}"\n'.format(disease), file=file)

    missing_features = set()
    fixed_caps = set()
    nift_value_originals = set(nift_values.values())

    for pathway, pathway_df in df.groupby(pathway_column):
        print('SET Subgraph = "{}"'.format(pathway), file=file)

        for _, gene, lit_snps, gwas_snps, clinical_features, clinical_snp in pathway_df[columns].itertuples():
            gene = ensure_quotes(gene)

            if lit_snps is None:
                lit_snps = []

            if gwas_snps is None:
                gwas_snps = []

            if clinical_snp is None:
                clinical_snp = []

            for snp in itt.chain(lit_snps, gwas_snps, clinical_snp):
                if not snp.strip():
                    continue
                print('g(HGNC:{}) -- g(dbSNP:{})'.format(gene, snp), file=file)

            for clinical_feature in clinical_features or []:
                if not clinical_feature.strip():
                    continue
                if clinical_feature.lower() not in nift_values:
                    missing_features.add(clinical_feature)
                    continue
                if clinical_feature not in nift_value_originals:
                    fixed_caps.add((clinical_feature, nift_values[clinical_feature.lower()]))
                    clinical_feature = nift_values[clinical_feature.lower()]  # fix capitalization
                print('g(HGNC:{}) -- a(NIFT:{})'.format(gene, ensure_quotes(clinical_feature)), file=file)

        print('UNSET Subgraph\n', file=file)

    print('UNSET MeSHDisease', file=file)
    print('UNSET Evidence', file=file)
    print('UNSET Citation', file=file)

    log.warning('Missing Features in %s', disease)
    for feature in missing_features:
        log.warning(feature)

    log.warning('Fixed capitalization')
    for broken, fixed in fixed_caps:
        log.warning('%s -> %s', broken, fixed)


if __name__ == '__main__':
    bms_base = os.environ['BMS_BASE']
    neurommsig_base = os.environ['NEUROMMSIG_BASE']
    neurommsig_excel_dir = os.path.join(neurommsig_base, 'resources', 'excels')

    nift_values = get_nift_values()

    ad_path = os.path.join(neurommsig_excel_dir, 'AD.xlsx')
    ad_df = preprocess(ad_path)
    with open(os.path.join(bms_base, 'aetionomy', 'alzheimers', 'neurommsigdb_ad.bel'), 'w') as ad_file:
        write_neurommsig_bel(ad_file, ad_df, mesh_alzheimer, nift_values)

    pd_path = os.path.join(neurommsig_excel_dir, 'PD.xlsx')
    pd_df = preprocess(pd_path)
    with open(os.path.join(bms_base, 'aetionomy', 'parkinsons', 'neurommsigdb_pd.bel'), 'w') as pd_file:
        write_neurommsig_bel(pd_file, pd_df, mesh_parkinson, nift_values)
