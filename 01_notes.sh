#!/usr/bin/env bash
# Daniel Stribling
# ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# hkref (hybkit Ref)  http://www.github.com/RenneLab/hkref 
# Hybkit Project : http://www.github.com/RenneLab/hybkit

NOTES="""
Download a reference sequence library for the Hyb program from Ensembl 
using the Biomart python module.

Library construction is based on the protocol provided in the supplemental methods of:
Helwak, Aleksandra, et al. 'Mapping the human miRNA interactome by CLASH reveals 
frequent noncanonical binding.' Cell 153.3 (2013): 654-665. 
http://dx.doi.org/10.1016/j.cell.2013.03.043 
( Supplemental methods section found only in PDF-fulltext )

Biomart queries include:
  protein_coding (as cDNA)
  lncRNA (as cDNA)
  All remaining gene_biotypes
      as unspliced transcripts ('transcript_exon_intron')

tRNAs:  genomic tRNA database http://gtrnadb.ucsc.edu/)
rRNAs:  NCBI Genbank Database, rRNA sequences (NR_003287.4, NR_003286.4); 
miRNAs: miRBase release 22.1 (http://www.mirbase.org): mature human miRNAs.

These sequences are then formatted in the required {}_{}_{name}_{biotype} header 
format for Hyb, and all extra '.' and '_' symbols are removed.

Original biotypes from the hOH7 Hyb database are:
Ig, lincRNA, microRNA, miscRNA, mRNA, mtrRNA, pr-tr, pseudo, rRNA, snoRNA, snRNA, Trec, tRNA
In this version, biotypes are passed through as with the ensembl 'transcript_biotype' field.

In order to facillitate unambiguous miRNA alignment, mature iRNA sequences are aligned to the 
reference transcriptome, and any alignemnts within transcripts are masked. This is performed to
ensure both that each given miRNA sequence has only a single reference alignment, as well as 
to allow miRNA precursor transcripts to be identified as hybrid targets.

"""
echo "${NOTES}"
