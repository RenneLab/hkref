#!/usr/bin/env python3
# Daniel Stribling
# ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# hkref (hybkit Ref)  http://www.github.com/RenneLab/hkref
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Download a reference sequence library for the Hyb program from Ensembl
using the Biomart python module.
"""

import sys
import os
import datetime
from Bio import SeqRecord, SeqIO, Seq

# try:
#     import hybkit
# except ModuleNotFoundError:
#     message = 'The "hybkit" module cannot be found. Please ensure this is accessible '
#     message += 'on your $PYTHONPATH'
#     print(message)
#     raise


# ---- Settings ----
# Import hkref project settings
from __settings__ import TOTAL_INDENT, BODY_INDENT, PROPERTY_INDENT, FINAL_WIDTH, \
                         DETAIL_DELIM, version

analysis_name = os.path.basename(__file__)
in_name = 'hkref_{}_with_coords.fa'.format(version)
out_name = 'hkref_{}.fa'.format(version)

print('\nFormatting Coordinates for consistency with Hyb')

all_ids = set()
dup_ids = set()
with open(in_name, 'r') as in_file:
    for line in in_file:
        if not line.startswith('>'):
            continue
        transcript = line.strip('>').strip().split(':')[0]
        if transcript in all_ids:
            dup_ids.add(transcript)
            #print('Duplicate:', transcript)
        all_ids.add(transcript)

with open(in_name, 'r') as in_file, open(out_name, 'w') as out_file:
    for seq_record in SeqIO.parse(in_file, "fasta"): 
        seq_record.description = ''
        if ':' in seq_record.id:
            transcript, coord = seq_record.id.split(':')
            if transcript in dup_ids:
                new_items = transcript.split('_')
                new_items[1] += '-' + coord
                seq_record.id = '_'.join(new_items)
            else:
                seq_record.id = transcript
            
        SeqIO.write(seq_record, out_file, "fasta")

out_notes_name = out_name + '.notes.txt'
with open(in_name + '.notes.txt', 'r') as notes_file:
    notes_str = notes_file.read()
notes_str += """

  -- {} --
    File: '{}' was created by {}
    (part of the hybkit project, release NA) on {}.
    Transcripts with only one portion had coordinates removed.=:
      "ENSG_ENST_GENE_TYPE:0-100" -> "ENSG_ENST_GENE_TYPE"
    Transcripts with multiple portions had coordinates added to the end 
    the transcript accession:
      "ENSG_ENST_GENE_TYPE:0-100" -> "ENSG_ENST-0-1100_GENE_TYPE"
""".format(out_name, out_name, analysis_name, str(datetime.datetime.now))

with open(out_notes_name, 'w') as out_file:
    out_file.write(notes_str)

print('\nDone\n')

