#!/usr/bin/env python3
# Daniel Stribling
# ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# hkref (hybkit Ref)  http://www.github.com/RenneLab/hkref 
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Settings for the hkref package.
"""

import os

project_id = 'hkref'
project_name = 'Hybkit-Ref'
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '_REF_VERSION.sh'), 'r') as version_file:
    version = version_file.read().split('NEW_DB_VER=')[1].split()[0].strip('"')
description = 'Up-to-date Genomic Sequence Reference Database for Hyb'
project_url = 'https://github.com/RenneLab/hkref'
keywords = 'genetics genomics ribonomics bioinformatics hyb CLASH qCLASH miRNA '
keywords += 'RNA DNA vienna viennad unafold'
name_and_version = project_name + '-' + version

TOTAL_INDENT = 2
BODY_INDENT = 2
PROPERTY_INDENT = 4
FINAL_WIDTH = 80

USE_ALL_NAMES = True
MAX_NAMES = 10

DETAIL_DELIM = '\t'

#Updated 2021-02-16
ensembl_release_num = '103'
ensembl_release_date = '2021-Feb'
ensembl_url = 'http://www.ensembl.org'

genbank_date_str = '202102'

feature_sets = [
    ['ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_transcript_id_version',
     'external_gene_name', 'external_transcript_name',
     'gene_biotype', 'transcript_biotype',],
    ['ensembl_transcript_id', 'refseq_mrna', 'refseq_ncrna',],
    ['ensembl_transcript_id', 'refseq_peptide',],
    ['ensembl_transcript_id', 'hgnc_symbol', 'hgnc_id'],
    ['ensembl_transcript_id', 'entrezgene_accession', 'entrezgene_description',],
    ['ensembl_transcript_id', 'mirbase_id', 'mirbase_accession'],
    ['ensembl_transcript_id',
     '5_utr_start',
     '5_utr_end',
     'cdna_coding_start',
     'cdna_coding_end',
     '3_utr_start',
     '3_utr_end',
     'ensembl_peptide_id',
     'strand',
     'transcript_start',
     'transcript_end',
     'transcription_start_site',
     'transcript_length',
     'cds_start',
     'cds_end',
    ]
]
#unused
 #'genedb',
 #'mirbase_trans_name',

sequence_attributes = [
    'ensembl_gene_id',
    'ensembl_transcript_id',
    'external_gene_name',
    'transcript_biotype',
    ] 

replace_biotypes = ['vtRNA', 'lncRNA', 'misc_RNA', 'Mt_rRNA', 'rRNA', 
                    'snoRNA', 'snRNA', 'Mt_tRNA',
                    'scaRNA', 'scRNA', 'sRNA', 'TEC', 'ribozyme',
                    'IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_V_gene',
                    'IG_C_pseudogene', 'IG_J_pseudogene', 'IG_pseudogene', 
                    'IG_V_pseudogene', 'polymorphic_pseudogene', 'processed_pseudogene', 
                    'rRNA_pseudogene', 'pseudogene', 'TR_J_pseudogene',
                    'TR_V_pseudogene', 'transcribed_processed_pseudogene', 
                    'transcribed_unitary_pseudogene', 'transcribed_unprocessed_pseudogene',
                    'translated_processed_pseudogene', 'translated_unprocessed_pseudogene', 
                    'unitary_pseudogene', 'unprocessed_pseudogene',
                    'TR_C_gene', 'TR_D_gene', 'TR_J_gene', 'TR_V_gene'
]

biotypes = { 
    'protein_coding':'mRNA',
}

replace_biotype_dict = {k: k.replace('_', '-') for k in replace_biotypes}
biotypes.update(replace_biotype_dict)

all_biotypes = [k for k in biotypes.keys()]
other_biotypes = [t for t in all_biotypes if t not in ['protein_coding', 'lncRNA']]

