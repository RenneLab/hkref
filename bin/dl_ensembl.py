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

import datetime
import textwrap
import pandas as pd
import sys
import os
import natsort
import pybiomart
import time
import requests
import biothings_client
import yaml
import argparse
import xml
import random
import signal
import functools

from Bio import SeqRecord, SeqIO, Seq

QUERY_TIMEOUT_MIN = 45
PRINT_ERROR = True

# ---- Ensembl Settings ----
feature_sets = [
    ['ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_transcript_id_version',
     'external_gene_name', 'external_transcript_name',
     'gene_biotype', 'transcript_biotype',],
    ['ensembl_transcript_id', 'embl'],
    ['ensembl_transcript_id', 'refseq_mrna'],
    ['ensembl_transcript_id', 'refseq_ncrna',],
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




# ---- DataFrame Helper Functions ----
def flatten_list(items):
    """Flatten any lists (depth 1) into parent list."""
    all_items = []
    for item in items:
        if isinstance(item, list):
            all_items += item
        else:
            all_items.append(item)
    return all_items


def join_unique(items, delim='|'):
    """Join unique entries of passed list of items using 'delim'."""
    return delim.join(natsort.natsorted({i for i in flatten_list(items) if i}))


def join_unique_ordered(items, delim='|'):
    """Join unique entries items in original order using 'delim'."""
    used = set()
    use_items = []
    for item in flatten_list(items):
        if item and item not in used:
            use_items.append(item)
            used.add(item)
    return delim.join(use_items)

def use_first_item(items, delim='|', default='UNKGENE'):
    """Uses first item from entries in order provided, or default."""
    for item in flatten_list(items):
        if item:
            if delim in str(item):
                use_item = item.split(delim)[0]
            else:
                use_item = item
            return use_item
    return default

def make_unique(in_str, delim='|'):
    """Split input string using delim and return corresponding string with unique items."""
    items = in_str.split(delim)
    return join_unique(items, delim)


def split_get_last(in_str, delim=':'):
    """Split by 'delim' and return last item."""
    return in_str.split(delim)[-1]


def remove_gene_version(name):
    """Remove ".##" portion of 'NAME.##' where it exists."""
    if '.' in name and name.split('.')[-1].isnumeric():
        return '.'.join(name.split('.')[:-1])
    else:
        return name


def sanitize_gene_name(all_names, delim='|', max_items=1000):
    """Remove gene version and disallowed characters from gene names."""
    use_names = []
    for name in all_names.split(delim):
        use_name = remove_gene_version(name)
        for c in ['_', '.', ' ', '\t', '\n', '[', ']', ':']:
            use_name = use_name.replace(c, '-')
        use_names.append(use_name)
    return delim.join(use_names[:max_items])

# Timeout methods via:
# https://stackoverflow.com/questions/2281850/timeout-function-if-it-takes-too-long-to-finish
class TimeoutError(Exception):
    pass

def timeout(minutes):
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(minutes * 60)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wrapper

    return decorator

# Example: Timeout after 5 minutes
# @timeout(5)
# def long_running_function2():
#     pass

def retry_query(dataset, *args, tries=40, timeout_min=QUERY_TIMEOUT_MIN, **kwargs):
    """Perform a dataset query, retrying several times as necessary."""
    @timeout(timeout_min)
    def query_with_timeout(dataset, *args, **kwargs):
        return dataset.query(*args, **kwargs)

    time.sleep(random.choice([0,1,2,3,4,5,6,7]))
    pfx = '    -'
    for i in range(1, tries):
        if i > 1:
            print(pfx, datetime.datetime.now(), 'Retrying...')
        else:
            print(pfx, datetime.datetime.now(), 'Querying...')
        sys.stdout.flush()
        try:
            query_df = query_with_timeout(dataset, *args, **kwargs)
            if '<html>' in query_df.columns:
                print(query_df.iloc[0:7])
                raise RuntimeError('Service Unavailable')
            print(pfx, datetime.datetime.now(), 'Success...')
            sys.stdout.flush()
            return query_df
        except requests.exceptions.HTTPError as e:
            err_message = 'HTTP Error'
            error_obj = e
        except requests.exceptions.ChunkedEncodingError as e:
            err_message = 'ChunkedEncodingError'
            error_obj = e
        except requests.exceptions.ConnectionError as e:
            err_message = 'ConnectionError'
            error_obj = e
        except RuntimeError as e:
            err_message = 'Service Unavailable'
            error_obj = e
        except xml.etree.ElementTree.ParseError as e:
            err_message = 'XML Parse Error'
            error_obj = e
        except TimeoutError as e:
            err_message = 'Timeout after %i Minutes.' % timeout_min
            error_obj = e

        next_sleep = ((i * 15) + random.choice([-5,-2,0,2,5]))
        print(pfx, datetime.datetime.now(), err_message + ', sleeping', next_sleep, 'seconds...')
        if PRINT_ERROR:
            print(error_obj)
        sys.stdout.flush()
        time.sleep(next_sleep)

    print(pfx, datetime.datetime.now(), 'Retrying... (Last)')
    sys.stdout.flush()
    return dataset.query(*args, **kwargs)

# ---- Create a text-wrapper for output messages ---
wrapper = textwrap.TextWrapper(width=100,
                               initial_indent='  ',
                               subsequent_indent=(' ' * 10))

# ---- Main Function ---
def dl_ensembl(db_settings_file, out_settings_file, use_query=None, use_cache=True):
    # ---- Read Database Settings ----
    with open(db_settings_file, "r") as db_settings_file_obj:
        db_settings = yaml.safe_load(db_settings_file_obj)

    # ---- Read Output Settings ----
    with open(out_settings_file, "r") as out_settings_file_obj:
        out_settings = yaml.safe_load(out_settings_file_obj)

    all_queries = {}
    biotype_map = {}
    if 'ensembl_biotype_sets' in db_settings:
        for name, set_str in db_settings['ensembl_biotype_sets'].items():
            all_queries[name] = { 'transcript_biotype': set_str.split() }
            for biotype in set_str.split():
                biotype_map[biotype] = biotype.replace('_', '-')

    # Set biotype map custom entries
    biotype_map['miRNA'] = 'pri-microRNA'
    biotype_map['protein_coding'] = 'mRNA'

    if use_query is not None:
        use_queries = {use_query: all_queries[use_query]}
    else:
        use_queries = all_queries

    out_name = 'hkref_%s_ensembl_%s' % (db_settings['title'], db_settings['ensembl_db_ver'])

    # ---- Perform Queries and Save Results ----
    print('\n --- Beginning Ensembl Queries ---')
    print('\nQuerying Ensembl DB at: {}'.format(db_settings['ensembl_db_url']))
    print('Using Cache:', use_cache)
    for query_name, filters in use_queries.items():
        # Setup query-specific paramaters
        detail_df = None
        seq_df = None
        query_strs = ''
        out_fasta_name = '{}_{}.fa'.format(out_name, query_name)
        out_fasta_notes_name = '{}_{}.fa.notes.txt'.format(out_name, query_name)
        out_detail_name = '{}_{}.csv'.format(out_name, query_name)
        out_detail_notes_name = '{}_{}.csv.notes.txt'.format(out_name, query_name)

        # Perform all queries and join results
        print('\nPerforming Query for Group: {}'.format(query_name))
        print('\n'.join(wrapper.wrap('- Using Filters:     ' + str(filters))))
        sys.stdout.flush()

        # Perform a separate query for each feature set, due to Biomart limitations
        dataset = pybiomart.Dataset(name=db_settings['ensembl_biomart_dataset'],
                                    host=db_settings['ensembl_db_url'],
                                    use_cache=use_cache
                                   )
        #print("DEBUG START")
        #if query_name in ['mRNA', 'lncRNA']:
        #    use_sequence_attributes = sequence_attributes + ['cdna']
        #else:
        #    use_sequence_attributes = sequence_attributes + ['transcript_exon_intron']

        #print('\n'.join(wrapper.wrap('- Querying Biomart Sequence: ' + str(use_sequence_attributes))))
        #print("DEBUG")
        #seq_df = retry_query(dataset, tries=3, attributes=use_sequence_attributes, filters=filters)
        #sys.exit()
        for i, feature_set in enumerate(feature_sets):
            print('\n'.join(wrapper.wrap('- Querying Biomart Features: ' + str(feature_set))))
            sys.stdout.flush()
            time.sleep(1)
            query_strs += 'filters={} | features={}\n'.format(str(filters), str(feature_set))

            # Perform biomart query
            query_df = retry_query(dataset, attributes=feature_set, filters=filters)

            # Based on selected features, some queries will return no information.
            #   if this occurs, add empty columns to query results.
            if query_df.empty:
                if i == 0:
                    raise Exception('Failed First Attempt.')
                new_columns = list(query_df.columns)
                for col in new_columns:
                    detail_df[col] = ''

            # If query has results, then add them to the details.
            else:
                query_df.fillna('', inplace=True)
                query_df = query_df.astype(str)
                query_df = query_df.groupby('Transcript stable ID').agg(join_unique)

                if detail_df is None:
                    detail_df = query_df
                else:
                    # assert detail_df.shape[0] == query_df.shape[0]
                    detail_df = detail_df.merge(query_df, how='outer',
                                                left_index=True, right_index=True,
                                                validate='1:1')
                    detail_df.fillna('', inplace=True)

        # If the query is mRNA, then remove all transcripts without a
        #   a corresponding RefSeq peptide ID.
        if query_name == 'mRNA':
            detail_df = detail_df.loc[detail_df['RefSeq peptide ID'] != '']
            keep_mRNA_seqs = {c for c in detail_df.index}

        # -- Perform a query to the MyGene.info service via biothings-client API -- #
        detail_df['HGNC ID num'] = detail_df['HGNC ID'].apply(split_get_last)
        hgnc_ids = list(detail_df['HGNC ID num'].unique())
        rename_columns = {
            'Gene name': 'Gene name Biomart',
            'HGNC symbol': 'HGNC symbol Biomart',
            'NCBI gene (formerly Entrezgene) accession': 'NCBI accession Biomart',
            'NCBI gene (formerly Entrezgene) description': 'NCBI description Biomart',
        }
        detail_df.rename(columns=rename_columns, inplace=True)

        print('  - Querying MyGene.info service for information on {} HGNC IDs'.format(len(hgnc_ids)))
        sys.stdout.flush()

        # biothings_scopes = 'symbol,name,alias,accession,hgnc,accession_genomic'
        # gene_fields = 'symbol,name,hgnc,alias,entrezgene,ensembl.gene,'
        # gene_fields += 'refseq,mirbase,type_of_gene,accession.genomic,other_names'
        mg = biothings_client.get_client('gene')
        # Name format: query_name: nice-name
        gene_fields_names = {'symbol': 'Symbol MG', 'name': 'Name MG', 'alias': 'Alias MG',
                             'entrezgene': 'NCBI gene MG', }
        use_names = list(gene_fields_names.values())
        gene_fields = ','.join(gene_fields_names.keys())
        gene_df = mg.querymany(hgnc_ids, scopes='hgnc', fields=gene_fields, species='human',
                               as_dataframe=True, index=True, verbose=False)

        # Rename to nice names, and merge new query information into details
        gene_df.rename(columns=gene_fields_names, inplace=True)
        detail_df = detail_df.merge(gene_df[use_names], how='outer',
                                    left_on='HGNC ID num', right_index=True,
                                    validate='m:1')
        detail_df.fillna('', inplace=True)

        # Create "first-symbol" field using all available gene identifiers.
        detail_df['First-Symbol'] = detail_df[['Gene name Biomart', 'HGNC symbol Biomart',
                                             'Symbol MG', 'Alias MG', 'European Nucleotide Archive ID']
                                            ].agg(use_first_item, axis='columns')

        # Create combined description field containing one of: NCBI / HGNC description.
        detail_df['Cmb-Description'] = detail_df['NCBI description Biomart']
        blank_description_rows = detail_df['Cmb-Description'] == ''
        detail_df.loc[blank_description_rows,
                      'Cmb-Description'
                      ] = detail_df.loc[blank_description_rows, 'Name MG']

        # Create translated "hyb" biotype
        detail_df['hyb_biotype'] = detail_df['Transcript type'].map(biotype_map)

        # Deprecated...
        # Depending on USE_ALL_NAMES setting, use either combined or original gene name
        #if USE_ALL_NAMES:
        #    def sanitize_with_maxlen(*args):
        #        return sanitize_gene_name(*args, max_items=MAX_NAMES)
        #    detail_df['cleaned_gene_names'] = detail_df['Cmb-Symbol'].apply(sanitize_with_maxlen)
        #else:
        #    detail_df['cleaned_gene_names'] = detail_df['Gene name'].apply(sanitize_gene_name)

        detail_df['cleaned_gene_names'] = detail_df['First-Symbol'].apply(sanitize_gene_name)

        # Create Hyb-style identifier for fasta sequences.
        detail_df['Transcript stable ID Col'] = detail_df.index
        detail_df['identifier'] = detail_df[['Gene stable ID', 'Transcript stable ID Col',
                                             'cleaned_gene_names', 'hyb_biotype']
                                            ].agg('_'.join, axis=1)
        detail_df.sort_values(['Gene stable ID', 'Transcript stable ID Col'], inplace=True)
        detail_df.drop(columns='Transcript stable ID Col', inplace=True)

        # Ensure transcript ID is not duplicated (from index)
        if 'Transcript stable ID' in detail_df.columns:
            detail_df.drop(columns='Transcript stable ID', inplace=True)
        print('  - Writing Detail Output to File:', out_detail_name)

        # Write output to detail file

        detail_df.to_csv(out_detail_name, sep=out_settings['detail_delim_char'])

        # -- Query biomart for sequences --
        if query_name in ['mRNA', 'lncRNA']:
            use_sequence_attributes = sequence_attributes + ['cdna']
        else:
            use_sequence_attributes = sequence_attributes + ['transcript_exon_intron']

        print('\n'.join(wrapper.wrap('- Querying Biomart Sequence: ' + str(use_sequence_attributes))))
        sys.stdout.flush()
        time.sleep(1)
        query_strs += 'filters={} | features={}\n'.format(str(filters), str(use_sequence_attributes))

        seq_df = retry_query(dataset, attributes=use_sequence_attributes, filters=filters)

        # If mRNA query, retain only where a RefSeq Peptide ID is associated
        if query_name == 'mRNA':
            print('  - Culling mRNA to only those with RefSeq Peptide ID')
            sys.stdout.flush()
            seq_df = seq_df.loc[seq_df['Transcript stable ID'].isin(keep_mRNA_seqs)]

        # Ensure same number of entries for sequence and details
        assert seq_df.shape[0] == detail_df.shape[0]
        seq_df = seq_df.join(detail_df['identifier'], on='Transcript stable ID', how='inner')
        assert seq_df.shape[0] == detail_df.shape[0]

        # Sort results
        seq_df.sort_values(['Gene stable ID', 'Transcript stable ID'], inplace=True)

        # Write fasta output
        print('  - Writing Fasta Output to File:', out_fasta_name)
        with open(out_fasta_name, 'w') as out_fasta:
            for i, seq_series in seq_df.iterrows():
                if 'cDNA sequences' in seq_series.index and seq_series['cDNA sequences']:
                    seq = seq_series['cDNA sequences']
                else:
                    seq = seq_series['Unspliced (Transcript)']
                record = SeqRecord.SeqRecord(Seq.Seq(seq),
                                             name=seq_series['identifier'],
                                             id=seq_series['identifier'],
                                             description=''
                                             )
                SeqIO.write(record, out_fasta, 'fasta')

        # Write notes files
        for data_file_name, notes_file_name in ((out_detail_name, out_detail_notes_name),
                                                (out_fasta_name, out_fasta_notes_name)):
            print('  - Writing details to file: %s' % notes_file_name)
            with open(notes_file_name, 'w') as out_notes_file:
                total_prefix = (' ' * out_settings['total_indent'])
                body_prefix = (' ' * out_settings['body_indent'])
                body_width = (out_settings['final_width'] -
                              - (out_settings['total_indent'] + out_settings['body_indent']))

                header = '-- %s -- ' % data_file_name
                description = 'File: "%s" ' % data_file_name
                description += 'was created by %s ' % (os.path.basename(__file__))
                description += '(part of the hybkit project, release %s)' % 'NA'
                # hybkit.__about__.__version__
                description += ' on %s.\n' % str(datetime.datetime.now().isoformat())
                description += 'Details for each sequence were downloaded from the biomart interface '
                description += 'to the Ensembl database, release %s ' % db_settings['ensembl_db_ver']
                description += '(%s), accessed via: "%s"\n' % (db_settings['ensembl_db_release_date'],
                                                               db_settings['ensembl_db_url'])
                description = textwrap.fill(description, width=body_width) + '\n'

                properties = 'The following queries were made:\n' + query_strs
                properties += '\n\n'

                body = textwrap.indent('\n'.join([description, properties]), body_prefix)
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
    parser.add_argument('--use_query', action='store',
                        help='Name of single query to perform if desired')
    parser.add_argument('--no_cache', dest='use_cache', action='store_false',
                        help='Disable use of pybiomart cache (to save space)')

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
    print('Running with Args:', args_info)
    dl_ensembl(**args_info)

