#!/usr/bin/env python3
# Daniel Stribling
# ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# hkref (hybkit Ref)  http://www.github.com/RenneLab/hkref
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Add notes files for non-ensembl databases.
"""

import datetime
import textwrap
import pandas as pd
import sys
import os
import biothings_client
import glob
from Bio import SeqRecord, SeqIO, Seq

# try:
#     import hybkit
# except ModuleNotFoundError:
#     message = 'The "hybkit" module cannot be found. Please ensure this is accessible '
#     message += 'on your $PYTHONPATH'
#     print(message)
#     raise

# ---- Settings ----

from __settings__ import TOTAL_INDENT, BODY_INDENT, PROPERTY_INDENT, FINAL_WIDTH, \
                         USE_ALL_NAMES, DETAIL_DELIM

genbank_fasta_name = glob.glob('hyb_genbank*.fa')[0]
trna_fasta_name = glob.glob('hyb_hg38*.fa')[0]
mirbase_fasta_name = glob.glob('hyb_miRBase*.fa')[0]
ex_ensembl_detail_name = glob.glob('hyb_ensembl*mRNA.csv')[0]

# ---- Write Detail files for fastas ----

# -- Get example header --
with open(ex_ensembl_detail_name, 'r') as ex_ensembl_detail:
    header_line = next(ex_ensembl_detail).rstrip()
header_items = header_line.split(DETAIL_DELIM)

print('\nWriting CSV Detail files for Non-Ensembl Databases.\n')
template_df = pd.DataFrame(columns=header_items)
# -- Query MyGene.info service --
# biothings_scopes = 'symbol,name,alias,accession,hgnc,accession_genomic'
# gene_fields = 'symbol,name,hgnc,alias,entrezgene,ensembl.gene,'
# gene_fields += 'refseq,mirbase,type_of_gene,accession.genomic,other_names'
mg = biothings_client.get_client('gene')
gene_fields_names = {'symbol': 'Symbol MG', 'name': 'Name MG', 'alias': 'Alias MG',
                     'entrezgene': 'NCBI gene MG', 'HGNC': 'HGNC ID num'}
use_field_names = [n for n in gene_fields_names.values()]
gene_fields = ','.join(gene_fields_names.keys())

# -- Write Genbank Details --
genbank_df = template_df.copy()
queries = []
with open(genbank_fasta_name, 'r') as genbank_fasta:
    for i, record in enumerate(SeqIO.parse(genbank_fasta, 'fasta')):
        identifier = record.id
        _, raw_transcript_id, gene_name, hyb_biotype = identifier.split('_')
        transcript_id = raw_transcript_id.replace('NR', 'NR_')
        queries.append(transcript_id)
        series_data = {
            'RefSeq ncRNA ID': transcript_id,
            'Cmb-Symbol': gene_name,
            'identifier': identifier,
            'cleaned_gene_names': gene_name,
            'hyb_biotype': hyb_biotype,
        }
        rec_series = pd.Series(series_data)
        genbank_df = genbank_df.append(rec_series, ignore_index=True)
genbank_df.fillna('', inplace=True)

# Query MyGene.info for additional details
gene_df = mg.querymany(queries, scopes='refseq.rna', fields=gene_fields, species='human',
                       as_dataframe=True, index=True, verbose=False)
gene_df.rename(columns=gene_fields_names, inplace=True)
genbank_df[use_field_names] = gene_df[use_field_names]

# Write details
genbank_detail_name = genbank_fasta_name.replace('.fa', '.csv')
print('Writing Genbank Details to:', genbank_detail_name)
genbank_df.to_csv(genbank_detail_name, sep=DETAIL_DELIM, index=False)

# -- Write miRBase Details --
mirbase_df = template_df.copy()
queries = []
with open(mirbase_fasta_name, 'r') as mirbase_fasta:
    for i, record in enumerate(SeqIO.parse(mirbase_fasta, 'fasta')):
        identifier = record.id
        mirbase_accession, _, mirbase_id, hyb_biotype = identifier.split('_')
        queries.append(mirbase_accession)
        series_data = {
            'miRBase ID': mirbase_id,
            'miRBase accession': mirbase_accession,
            'identifier': identifier,
            'cleaned_gene_names': mirbase_id,
            'hyb_biotype': hyb_biotype,
        }
        rec_series = pd.Series(series_data)
        mirbase_df = mirbase_df.append(rec_series, ignore_index=True)
mirbase_df.fillna('', inplace=True)

# gene_df = mg.querymany(queries, scopes='miRBase', fields=gene_fields, species='human',
#                            as_dataframe=True, index=True, verbose=False)
# print(gene_df)
# gene_df.rename(columns=gene_fields_names, inplace=True)
# mirbase_df[use_field_names] = gene_df[use_field_names]

# Write miRBase details
mirbase_detail_name = mirbase_fasta_name.replace('.fa', '.csv')
print('Writing miRBase Details to:', mirbase_detail_name)
mirbase_df.to_csv(mirbase_detail_name, sep=DETAIL_DELIM, index=False)

# -- Write tRNA Details --
trna_df = template_df.copy()
queries = []
with open(trna_fasta_name, 'r') as trna_fasta:
    for i, record in enumerate(SeqIO.parse(trna_fasta, 'fasta')):
        identifier = record.id
        _, _, trna_id, hyb_biotype = identifier.split('_')
        queries.append(trna_id)
        series_data = {
            'identifier': identifier,
            'cleaned_gene_names': trna_id,
            'hyb_biotype': hyb_biotype,
        }
        rec_series = pd.Series(series_data)
        trna_df = trna_df.append(rec_series, ignore_index=True)
trna_df.fillna('', inplace=True)

# gene_df = mg.querymany(queries, scopes='miRBase', fields=gene_fields, species='human',
#                            as_dataframe=True, index=True, verbose=False)
# print(gene_df)
# gene_df.rename(columns=gene_fields_names, inplace=True)
# trna_df[use_field_names] = gene_df[use_field_names]

# Write tRNA Details
trna_detail_name = trna_fasta_name.replace('.fa', '.csv')
print('Writing tRNA Details to:', trna_detail_name)
trna_df.to_csv(trna_detail_name, sep=DETAIL_DELIM, index=False)

# ---- Write notes files for all detail files ----
for out_file_name in (genbank_detail_name, mirbase_detail_name, trna_detail_name):
    notes_file_name = out_file_name + '.notes.txt'
    print('  - Writing notes to file: %s' % notes_file_name)
    with open(notes_file_name, 'w') as out_notes_file:
        total_prefix = (' ' * TOTAL_INDENT)
        body_prefix = (' ' * BODY_INDENT)
        property_prefix = (' ' * PROPERTY_INDENT)
        total_width = FINAL_WIDTH - (TOTAL_INDENT)
        body_width = FINAL_WIDTH - (TOTAL_INDENT + BODY_INDENT)
        property_width = FINAL_WIDTH - (TOTAL_INDENT + BODY_INDENT + PROPERTY_INDENT)

        header = '-- %s -- ' % out_file_name
        description = 'File: "%s" ' % out_file_name
        description += 'was created by %s ' % (os.path.basename(__file__))
        description += '(part of the hybkit project, release %s)' % 'NA'
        # hybkit.__about__.__version__
        description += ' on %s.\n' % str(datetime.datetime.now().isoformat())
        description = textwrap.fill(description, width=body_width) + '\n'
        properties = ''

        body = textwrap.indent('\n'.join([description, properties]), body_prefix)
        full_text = textwrap.indent('\n'.join([header, body]), total_prefix)

        out_notes_file.write(full_text)

print('Done.\n')
