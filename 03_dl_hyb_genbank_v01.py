#!/usr/bin/env python3
# Daniel Stribling
# ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# hkref (hybkit Ref)  http://www.github.com/RenneLab/hkref
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Download two rRNA Reference sequences from GenBank using BioPython
"""

import datetime
import textwrap
import sys
import os
import natsort
import pandas as pd
from Bio import Entrez, SeqIO

# try:
#     import hybkit
# except ModuleNotFoundError:
#     message = 'The "hybkit" module cannot be found. Please ensure this is accessible '
#     message += 'on your $PYTHONPATH'
#     print(message)
#     raise

from __settings__ import TOTAL_INDENT, BODY_INDENT, PROPERTY_INDENT, FINAL_WIDTH, \
                         genbank_date_str

# ---- Settings ----

seq_accessions = ['NR_003287.4', 'NR_003286.4']

out_name = 'hyb_genbank_{}'.format(genbank_date_str)
out_fasta_name = out_name + '.fa'
out_fasta_notes_name = out_fasta_name + '.notes.txt'

# ---- Download Genbank Sequences ----

print('\nBeginning Genbank Downloads.\n')

if len(sys.argv) < 2:
    Entrez.email = input('Please enter your email address for the Entrez API: ')
else:
    Entrez.email = sys.argv[1]

records = []
for seq_accession in seq_accessions:
    call_params = {
                   'db': 'nuccore',
                   'id': seq_accession,
                   'rettype': 'fasta',
                   'style': 'full',
                   'retmode': 'text',
                  }

    with Entrez.efetch(**call_params) as stream_handle:
        record = SeqIO.read(stream_handle, 'fasta')
        gene_name = record.description.split('(')[1].split(')')[0]
        header = 'Genbank_{}_{}_rRNA'.format(seq_accession.split('.')[0].replace('_', ''),
                                             gene_name)
        record.id = record.name = header
        record.description = ''
        records.append(record)

# ---- Write FASTA output ----

with open(out_fasta_name, 'w') as out_fasta:
    SeqIO.write(records, out_fasta, 'fasta')

# ---- Write notes file ----

print('Writing details to file: %s\n' % out_fasta_notes_name)

with open(out_fasta_notes_name, 'w') as out_notes_file:
    total_prefix = (' ' * TOTAL_INDENT)
    body_prefix = (' ' * BODY_INDENT)
    property_prefix = (' ' * PROPERTY_INDENT)
    total_width = FINAL_WIDTH - (TOTAL_INDENT)
    body_width = FINAL_WIDTH - (TOTAL_INDENT + BODY_INDENT)
    property_width = FINAL_WIDTH - (TOTAL_INDENT + BODY_INDENT + PROPERTY_INDENT)

    header = '-- %s -- ' % out_fasta_name
    description = 'File: "%s" ' % out_fasta_name
    description += 'was created by %s ' % (os.path.basename(__file__))
    description += '(part of the hybkit project, release %s)' % 'NA'
    # hybkit.__about__.__version__
    description += ' on %s.\n' % str(datetime.datetime.now().isoformat())
    description += 'Details for each sequence were downloaded from the '
    description += 'Biopython Entrez interface '
    description += 'to the NCBI Genbank database.\n'
    description = textwrap.fill(description, width=body_width) + '\n'

    body = textwrap.indent(description, body_prefix)
    full_text = textwrap.indent('\n'.join([header, body]), total_prefix)

    out_notes_file.write(full_text)

print('Done.\n')
