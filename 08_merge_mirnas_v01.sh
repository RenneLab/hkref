#!/usr/bin/env bash
# Daniel Stribling
# ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# hkref (hybkit Ref)  http://www.github.com/RenneLab/hkref 
# Hybkit Project : http://www.github.com/RenneLab/hybkit

set -e -u -o pipefail
NOTES="""
Combine individual RNA fasta components into combined transcript reference.
"""

SEQKIT_VERSION="seqkit/0.10.2"

# Load NEW_DB_VER and NEW_DB_NAME variables:
if [ -f ./_REF_VERSION.sh ]; then
    source ./_REF_VERSION.sh
elif [ -f ../_REF_VERSION.sh ]; then
    source ../_REF_VERSION.sh
else
    echo "Can't find _REF_VERSION.sh"
    exit 1
fi 
CMB_FA="${NEW_DB_NAME}_with_coords.fa"
CMB_DETAIL="${NEW_DB_NAME}.csv"

# Get names of all output files
#FASTA_NAME_STEMS="hyb_ensembl*.fa  hyb_genbank_*.fa hyb_*tRNAs*.fa  hyb_miRBase*.fa"
FASTA_NAME_STEMS="${NEW_DB_NAME}_*_ext.fa hyb_miRBase*.fa "
DETAIL_NAME_STEMS="${NEW_DB_NAME}*_only.csv hyb_miRBase*.csv "
FASTA_NAMES=""
FASTA_NOTES_NAMES=""
DETAIL_NAMES=""
DETAIL_NOTES_NAMES=""
VERSION_STR=""

for file in $(find ${FASTA_NAME_STEMS}); do
  FASTA_NAMES="${FASTA_NAMES} ${file}"
  FASTA_NOTES_NAMES="${FASTA_NOTES_NAMES} ${file/.fa/.fa.notes.txt}"
done
for file in $(find ${DETAIL_NAME_STEMS}); do
  DETAIL_NAMES="${DETAIL_NAMES} ${file}"
  DETAIL_NOTES_NAMES="${DETAIL_NOTES_NAMES} ${file}.notes.txt"
done

# Combine FASTA files
echo -e "\nCombining Files: ${FASTA_NAMES}"
echo -e "Into File: ${CMB_FA}"
set -v
cat ${FASTA_NAMES} > ${CMB_FA} 
set +v

# Combine FASTA notes files
echo -e "\nCombining Notes Files: ${FASTA_NOTES_NAMES}"
echo -e "Into File: ${CMB_FA}.notes.txt"
set -v
cat ${FASTA_NOTES_NAMES} > ${CMB_FA}.notes.txt 
set +v

echo """

  -- ${CMB_FA} --
    File: '${CMB_FA}' was created by $0
    (part of the hybkit project, release NA) on $(date).
    by concatenating files: ${FASTA_NAMES} 
""" >> ${CMB_FA}.notes.txt

# Combine Detail Files
echo -e "\nCombining Detail Files: ${DETAIL_NAMES}"
echo -e "Into File: ${CMB_DETAIL}"
set -v
awk 'FNR==1 && NR!=1{next;}{print}' ${DETAIL_NAMES} > ${CMB_DETAIL}
set +v

# Combine Detail Notes Files
echo -e "\nCombining Detail Notes Files: ${DETAIL_NOTES_NAMES}"
echo -e "Into File: ${CMB_DETAIL}.notes.txt"
set -v
cat ${DETAIL_NOTES_NAMES} > ${CMB_DETAIL}.notes.txt 
set +v

echo """

  -- ${CMB_DETAIL} --
    File: '${CMB_DETAIL}' was created by $0
    (part of the hybkit project, release NA) on $(date).
    by concatenating files: ${DETAIL_NAMES} 
""" >> ${CMB_DETAIL}.notes.txt

# -- Final sorting is disabled.

# -- Sort output sequences
#if [ ${USE_MODULE} > 0 ]; then
#  module load ${SEQKIT_VERSION}
#  VERSION_STR="${SEQKIT_VERSION}"
#fi
#set -v -H -o history
#seqkit sort --by-name ${CMB_FA} > ${CMB_FA_SORT}
#SEQKIT_COMMAND="!!"
#set +v +H +o history
#cp -v ${CMB_FA}.notes.txt ${CMB_FA_SORT}.notes.txt
#
#echo """
#
#  -- ${CMB_FA_SORT} --
#    File: '${CMB_FA_SORT}' was created by $0
#    (part of the hybkit project, release NA) on $(date).
#    by using seqkit ${VERSION_STR} command:
#    ${SEQKIT_COMMAND}
#""" >> ${CMB_FA_SORT}.notes.txt

# Deduplication is DISABLED
#NODUP_FASTA_NAME_STEMS=""

#MEM_MB_AVAIL=7000
#CPUS_AVAIL=8

#echo -e "\nRemoving Duplicates from Files:\n$(find ${FASTA_NAME_STEMS})"
#for in_file in ${FASTA_NAME_STEMS}; do
#    echo -e "\nRemoving duplicates and subsequences from file: ${in_file}"
#    cd-hit -c 1.0 -t 1 -T ${CPUS_AVAIL} -M ${MEM_MB_AVAIL} -uS 0.0 -i ${in_file} -o ${in_file/.fa/_nodup.fa}
#    # -c  : required identity proportion
#    # -uS : Maximum percent not aligned of smaller read
#    NODUP_FASTA_NAME_STEMS="${NODUP_FASTA_NAME_STEMS} ${in_file/.fa/_nodup.fa}"
#done
#echo -e "\nCombining Files:\n$(find ${NODUP_FASTA_NAME_STEMS})"
#echo -e "Into File: ${CMB_FA}"
#cat ${NODUP_FASTA_NAME_STEMS} > ${CMB_FA} 

#ZIP_FILES="${CMB_FA}  ${CMB_FA/.fa/.csv}"
#echo "Zipping Files: ${ZIP_FILES}"
#gzip --best ${ZIP_FILES}

echo -e "\nDone.\n"
