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

Biomart queries include:
  mRNA (protein_coding; as cDNA) where a RefSeq Protein Identifier Exists
  lncRNA (as cDNA)
  All remaining gene_biotypes, excluding 'miRNA', 
      as unspliced transcripts ('transcript_exon_intron')

tRNAs:  genomic tRNA database http://gtrnadb.ucsc.edu/)
rRNAs:  NCBI Genbank Database, rRNA sequences (NR_003287.4, NR_003286.4); 
miRNAs: miRBase release 22.1 (http://www.mirbase.org): mature human miRNAs.

These sequences are then formatted in the required {}_{}_{name}_{biotype} header 
format for Hyb, and all extra '.' and '_' symbols are removed.

Original biotypes from the hOH7 Hyb database are:
Ig, lincRNA, microRNA, miscRNA, mRNA, mtrRNA, pr-tr, pseudo, rRNA, snoRNA, snRNA, Trec, tRNA
Other types are passed through as with the ensembl 'transcript_biotype' field.

The original protocol deduplicated sequences and removed subsequences from the reference. 
In this protocol, these steps are omitted to include all possible annotations and splicing
variants in the analysis.

"""
echo "${NOTES}"
