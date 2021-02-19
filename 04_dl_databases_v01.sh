#!/usr/bin/env bash
# Daniel Stribling
# ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# hkref (hybkit Ref)  http://www.github.com/RenneLab/hkref 
# Hybkit Project : http://www.github.com/RenneLab/hybkit

set -e -u -o pipefail
NOTES="""
Download current (non-ensembl) databases for use and prepare.
Utilizes seqkit as a dependency.

tRNAs:  genomic tRNA database http://gtrnadb.ucsc.edu/)
miRNAs: miRBase release 22.1 (http://www.mirbase.org): mature human miRNAs.
"""

TRNA_VER="2013-12"
TRNA_SOURCE="Genomic tRNA Database"
TRNA_URL="http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz"

MIRNA_VER="R22.1-2018-Oct"
MIRNA_URL="ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz"
MIRNA_SOURCE="miRBase"

TRNA_FINAL_FILE="hyb_hg38-tRNAs_${TRNA_VER}.fa"
MIRNA_FINAL_FILE="hyb_miRBase_${MIRNA_VER}.fa"

echo "$NOTES"

# Download and prepare tRNA sequences
echo "Downloading tRNAs from ${TRNA_URL}"
wget -vc ${TRNA_URL}
tar -xvf hg38-tRNAs.tar.gz hg38-tRNAs.fa 
mv hg38-tRNAs.fa hg38-tRNAs_2013-12.fa
rm hg38-tRNAs.tar.gz
seqkit replace -p "\s.+" hg38-tRNAs_2013-12.fa | \
    seqkit replace -p "_" | \
    seqkit replace -p "Homosapiens" -r "tRNASCAN_tRNASeq_" | \
    seqkit replace -p "$" -r "_tRNA" | \
    seqkit sort -N -o ${TRNA_FINAL_FILE}
rm hg38-tRNAs_2013-12.fa

echo """
  -- ${TRNA_FINAL_FILE} -- 
    File: '${TRNA_FINAL_FILE}' was created by $0
    (part of the hybkit project, release NA) on $(date).
    tRNA sequences (version ${TRNA_VER}) were downloaded from the ${TRNA_SOURCE}:
    http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz
    and identifiers were renamed to match Hyb convention.
""" > ${TRNA_FINAL_FILE}.notes.txt

# Download and prepare miRNA sequences
echo "Downloading miRNAs from ${MIRNA_URL}"
wget -vc ${MIRNA_URL}
gunzip -f mature.fa.gz
mv mature.fa miRBase_mature.fa
seqkit grep -r -p "^hsa" miRBase_mature.fa | \
  seqkit seq --rna2dna | \
  seqkit replace -p " Homo sapiens " -r "_miRBase_hsa-" | \
  seqkit replace -p ".+MIMAT" -r "MIMAT" | \
  seqkit replace -p "$" -r "_microRNA" | \
  seqkit sort -N -o ${MIRNA_FINAL_FILE}
rm miRBase_mature.fa

echo """
  -- ${MIRNA_FINAL_FILE} -- 
    File: '${MIRNA_FINAL_FILE}' was created by $0
    (part of the hybkit project, release NA) on $(date).
    Mature miRNA sequences (version ${MIRNA_VER}) were downloaded from ${MIRNA_SOURCE}:
    ${MIRNA_URL}
    were converted to DNA alphabet, 
    and identifiers were renamed to match Hyb convention.
""" > ${MIRNA_FINAL_FILE}.notes.txt

echo -e "\nDone\n"
