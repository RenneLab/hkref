#!/usr/bin/env python3
# Daniel Stribling
# ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# hkref (hybkit Ref)  http://www.github.com/RenneLab/hkref
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Merge sequences split by seqkit sliding window into a single sequence.
"""

import os
import sys
import yaml
import textwrap
import datetime
import argparse
from Bio import SeqRecord, SeqIO, Seq

SCRIPT_NAME = os.path.basename(__file__)

# ---- Main Function ---
def merge_split_seqs(in_file_name, out_file_name):
    """Merge sequences split by seqkit sliding window into a single sequence."""
    # ---- Read Input File and Identify Duplicate Transcript IDs ----
    all_ids = set()
    dup_ids = set()
    last_id = ''

    print('\nReading input file and identifying duplicate transcript ids...')
    sys.stdout.flush()
    with open(in_file_name, 'r') as in_file_obj:
        for line in in_file_obj:
            if not line.startswith('>'):
                continue
            # Check for more than one occurrence of character ':' in line
            if line.count(':') > 1:
                message = 'ERROR: More than one ":" in transcript id: {}'.format(line.strip())
                raise ValueError(message)
            transcript_id = line.lstrip('>').strip().split(':')[0]
            if transcript_id in all_ids:
                dup_ids.add(transcript_id)
                # Check to ensure duplicates are in sequential order
                if transcript_id != last_id:
                    message = 'ERROR: Duplicate transcript id "{}"'.format(transcript_id)
                    message += ' is not in sequential order.'
                    message += '\nPlease ensure transcripts are sorted in name order before '
                    message += 'running this script.'
                    raise ValueError(message)
            all_ids.add(transcript_id)
            last_id = transcript_id

    # ---- Read Input File, Merge Duplicates, and Output with Corrected Identifiers
    print('\nReading input file, merging duplicates, and outputting with corrected identifiers...')
    print('  Input file: {}'.format(in_file_name))
    print('  Output file: {}'.format(out_file_name))
    sys.stdout.flush()
    with open(in_file_name, 'r') as in_file_obj, open(out_file_name, 'w') as out_file_obj:
        pending_seq = None
        pending_last_coord = 0
        pending_last_id = ''
        for seq_record in SeqIO.parse(in_file_obj, "fasta"):
            seq_record.description = ''
            # Where record ids have coordinates, reformat and/or join.
            if ':' in seq_record.id:
                transcript_id, coord = seq_record.id.split(':')
                first_coord, last_coord = int(coord.split('-')[0]), int(coord.split('-')[1])
                is_first_transcript = coord.startswith('1-')
                # If a previous pending sequence exists with a different ID,
                # write it to the output file.
                if pending_last_id and pending_last_id != transcript_id:
                    SeqIO.write(pending_seq, out_file_obj, "fasta")
                    pending_seq = None
                    pending_last_coord = 0
                    pending_last_id = ''

                # If transcript is split/duplicate, then merge into pending sequence.
                if transcript_id in dup_ids:
                    if pending_seq is None:
                        seq_record.id = transcript_id
                        pending_seq = seq_record
                        pending_last_coord = last_coord
                        pending_last_id = transcript_id
                    else:
                        if not first_coord == pending_last_coord + 1:
                            message = 'ERROR: Coordinates are not sequential for transcript '
                            message += '"{}"'.format(transcript_id)
                            message += '\nPlease ensure transcripts are sorted in name order'
                            message += ' before running this script.'
                            raise ValueError(message)
                        pending_seq.seq += seq_record.seq
                        pending_last_coord = last_coord
                        pending_last_id = transcript_id

                # If transcript is not split/duplicate and is first transcript, remove coordinates.
                elif is_first_transcript:
                    seq_record.id = transcript_id

                # If not split/duplicate and not first transcript, modify coordinates
                #   into transcript portion of identifier
                else:
                    new_items = transcript_id.split('_')
                    new_items[1] += '-' + coord
                    seq_record.id = '_'.join(new_items)

            if pending_seq is None:
                SeqIO.write(seq_record, out_file_obj, "fasta")

        # Write any final pending sequence to the output file.
        if pending_last_id and pending_last_id != transcript_id:
            SeqIO.write(pending_seq, out_file_obj, "fasta")


# ---- Notes Function ----
def write_notes_file(in_notes_name, out_notes_name, out_file_name, out_settings_file_name):
    """Write notes file for this script."""
    # ---- Read Output Settings ----
    with open(out_settings_file_name, "r") as out_settings_file_obj:
        out_settings = yaml.safe_load(out_settings_file_obj)

    with open(in_notes_name, 'r') as in_notes_obj:
        in_notes_str = in_notes_obj.read()

    total_prefix = (' ' * out_settings['total_indent'])
    body_prefix = (' ' * out_settings['body_indent'])
    body_width = (out_settings['final_width'] -
                     (out_settings['total_indent']
                       + out_settings['body_indent']))

    header = '-- {} --\n'.format(out_file_name)
    notes_str = "File: '{}' was created by {}"
    notes_str += ' (part of the hybkit project, release NA) on {}.\n'
    notes_str += 'Transcripts with only one portion had coordinates removed.=\n'
    notes_str += '  "ENSG_ENST_GENE_TYPE:0-100" -> "ENSG_ENST_GENE_TYPE"\n'
    notes_str += 'Split Transcripts were joined into a single sequence and coordinates removed.\n'
    notes_str += 'Transcripts with multiple portions had coordinates added to the end\n'
    notes_str += 'the transcript accession:\n'
    notes_str += '  "ENSG_ENST_GENE_TYPE:0-100" -> "ENSG_ENST-0-1100_GENE_TYPE"\n'
    notes_str = notes_str.format(out_file_name,
                                 SCRIPT_NAME, str(datetime.datetime.now()))
    notes_str = textwrap.fill(notes_str, width=body_width)
    notes_str = textwrap.indent(notes_str, body_prefix)
    add_str = textwrap.indent(header + notes_str, total_prefix)
    write_str = in_notes_str + '\n\n' + add_str
    with open(out_notes_name, 'w') as out_notes_file:
        out_notes_file.write(write_str)





# ---- Argument Parsing Function ----
def parse_args(print_help=False):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    #                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    in_opts = parser.add_argument_group('Input Options')
    parser.add_argument('--in_file', action='store', required=True,
                        help='Input sequence file name')
    parser.add_argument('--in_notes', action='store', required=True,
                        help='Input notes file name')
    parser.add_argument('--out_settings_file', action='store', required=True,
                        help='YAML file containing output settings')

    out_opts = parser.add_argument_group('Output Options')
    parser.add_argument('--out_file', action='store', required=True,
                        help='Output sequence file name')
    parser.add_argument('--out_notes', action='store', required=True,
                        help='Output notes file name')


    if print_help:
        parser.print_help()
        return {}

    namespace = parser.parse_args()
    args_info = vars(namespace)
    if not args_info:
        sys.exit()

    return args_info


# ---- Main Function Call ----
if __name__ == '__main__':
    args_info = parse_args()
    print('\nRunning {} with:\n'.format(SCRIPT_NAME))
    for key, value in args_info.items():
        print('  {:<20}  {}'.format(key, value))

    merge_split_seqs(args_info['in_file'], args_info['out_file'])
    write_notes_file(args_info['in_notes'], args_info['out_notes'],
                     args_info['out_file'], args_info['out_settings_file'])

    print('\nDone\n')

    # Example System call with base path [WORKDIR]:
    # python3 merge_split_seqs.py \
    #       --in_file [WORKDIR]/in_fasta.fa \
    #       --in_notes [WORKDIR]/in_fasta.notes.txt \
    #       --out_file [WORKDIR]/out_fasta.fa \
    #       --out_settings_file [WORKDIR]/out_settings.yaml \
    #       --out_notes [WORKDIR]/out_fasta.notes.txt




