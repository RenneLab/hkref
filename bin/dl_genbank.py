#!/usr/bin/env python3
# Daniel Stribling
# ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# hkref (hybkit Ref)  http://www.github.com/RenneLab/hkref
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Download sequences from GenBank using BioPython
"""

import argparse
import datetime
import textwrap
import sys
import os
import natsort
import pandas as pd
import yaml
from Bio import Entrez, SeqIO


# ---- Download Genbank Sequences ----
def dl_genbank(db_settings_file, out_settings_file, user_email):
    # ---- Read Database Settings ----
    with open(db_settings_file, "r") as db_settings_file_obj:
        db_settings = yaml.safe_load(db_settings_file_obj)
    
    # ---- Read Output Settings ----
    with open(out_settings_file, "r") as out_settings_file_obj:
        out_settings = yaml.safe_load(out_settings_file_obj)

    out_name = 'hkref_%s_genbank_%s' % (db_settings['title'], db_settings['genbank_db_acc_date'])
    out_fasta_name = out_name + '.fa'
    out_fasta_notes_name = out_fasta_name + '.notes.txt'
    
    print('\nBeginning Genbank Downloads.\n')
    Entrez.email = user_email
    
    records = []
    for seq_type in db_settings['genbank_accessions']:
        for seq_accession in db_settings['genbank_accessions'][seq_type]:
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
                header = 'Genbank_{}_{}_{}'.format(seq_accession.split('.')[0].replace('_', ''),
                                                   gene_name,
                                                   seq_type)
                record.id = record.name = header
                record.description = ''
                records.append(record)
    
    # ---- Write FASTA output ----
    
    with open(out_fasta_name, 'w') as out_fasta:
        SeqIO.write(records, out_fasta, 'fasta')
    
    # ---- Write notes file ----
    
    print('Writing details to file: %s\n' % out_fasta_notes_name)
    
    with open(out_fasta_notes_name, 'w') as out_notes_file:
        total_prefix = (' ' * out_settings['total_indent'])
        body_prefix = (' ' * out_settings['body_indent'])
        property_prefix = (' ' * out_settings['property_indent'])
        total_width = (out_settings['final_width'] - (out_settings['total_indent']))
        body_width = (out_settings['final_width'] 
                      - (out_settings['total_indent'] + out_settings['body_indent'])
                     )
        property_width = (out_settings['final_width'] 
                          - (out_settings['total_indent'] + out_settings['body_indent'] 
                             + out_settings['property_indent']
                            )
                         )
    
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

def parse_args(print_help=False):
    parser = argparse.ArgumentParser(description=__doc__)
    #                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    in_opts = parser.add_argument_group('Input Options')
    parser.add_argument('--db_settings_file', action='store', required=True,
                        help='YAML file containing database settings')
    parser.add_argument('--out_settings_file', action='store', required=True,
                        help='YAML file containing output settings')
    parser.add_argument('--user_email', action='store', required=True,
                        help='Email address of user for Entrez API')

    if print_help:
        parser.print_help()
        return {}

    namespace = parser.parse_args()
    args_info = vars(namespace)
    if not args_info:
        sys.exit()

    return args_info

if __name__ == '__main__':
    args_info = parse_args()
    dl_genbank(**args_info)
