#!/usr/bin/env python3
# Daniel Stribling
# ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# hkref (hybkit Ref)  http://www.github.com/RenneLab/hkref 
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Settings for the hkref package.
"""

TOTAL_INDENT = 2
BODY_INDENT = 2
PROPERTY_INDENT = 4
FINAL_WIDTH = 80

USE_ALL_NAMES = True

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

biotypes = { 
    'vault_RNA':'vtRNA',
    'lncRNA':'lincRNA',
    'misc_RNA':'miscRNA',
    'protein_coding':'mRNA',
    'Mt_rRNA':'mtrRNA',
    #'':'pr-tr',
    #'':'pseudo',
    'rRNA':'rRNA',
    'snoRNA':'snoRNA',
    'snRNA':'snRNA',
    #'':'microRNA',
    'Mt_tRNA':'tRNA',
}

ig_biotypes = {k:'Ig' for k in ['IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_V_gene']}

pseudo_biotypes = {k:'pseudo' for k in ['IG_C_pseudogene', 'IG_J_pseudogene', 
                                    'IG_pseudogene', 'IG_V_pseudogene',
                                    'polymorphic_pseudogene', 'processed_pseudogene', 
                                    'rRNA_pseudogene', 'pseudogene', 'TR_J_pseudogene',
                                    'TR_V_pseudogene', 'transcribed_processed_pseudogene', 
                                    'transcribed_unitary_pseudogene', 
                                    'transcribed_unprocessed_pseudogene', 
                                    'translated_processed_pseudogene', 
                                    'translated_unprocessed_pseudogene', 
                                    'unitary_pseudogene', 'unprocessed_pseudogene'
                                   ]}

trec_biotypes = {k:'Trec' for k in ['TR_C_gene', 'TR_D_gene', 'TR_J_gene', 'TR_V_gene']}

other_biotypes = {k:k for k in ['scaRNA', 'scRNA', 'sRNA', 'TEC', 'ribozyme']}

for add_biotype in [ig_biotypes, pseudo_biotypes, trec_biotypes, other_biotypes]:
    biotypes.update(add_biotype)

all_biotypes = [k for k in biotypes.keys()]
other_biotypes = [t for t in all_biotypes if t not in ['protein_coding', 'lncRNA']]
