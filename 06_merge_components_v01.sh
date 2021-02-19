#!/usr/bin/env bash
# Daniel Stribling
# ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# hkref (hybkit Ref)  http://www.github.com/RenneLab/hkref 
# Hybkit Project : http://www.github.com/RenneLab/hybkit

set -e -u -o pipefail
NOTES="""
Combine individual fasta components into combined Hyb reference.
"""

# Load NEW_DB_VER and NEW_DB_NAME variables:
if [ -f ./_REF_VERSION.sh ]; then
    source ./_REF_VERSION.sh
elif [ -f ../_REF_VERSION.sh ]; then
    source ../_REF_VERSION.sh
else
    echo "Can't find _REF_VERSION.sh"
    exit 1
fi 
NEW_DB_FILE="${NEW_DB_NAME}.fa"

# Get names of all output files
FASTA_NAME_STEMS="hyb_ensembl*.fa  hyb_genbank_*.fa hyb_*tRNAs*.fa  hyb_miRBase*.fa"
FASTA_NAMES=""
FASTA_NOTES_NAMES=""
DETAIL_NAMES=""
DETAIL_NOTES_NAMES=""

for file in $(find ${FASTA_NAME_STEMS}); do
  FASTA_NAMES="${FASTA_NAMES} ${file}"
  FASTA_NOTES_NAMES="${FASTA_NOTES_NAMES} ${file/.fa/.fa.notes.txt}"
  DETAIL_NAMES="${DETAIL_NAMES} ${file/.fa/.csv}"
  DETAIL_NOTES_NAMES="${DETAIL_NOTES_NAMES} ${file/.fa/.csv.notes.txt}"
done

# Combine FASTA files
echo -e "\nCombining Files: ${FASTA_NAMES}"
echo -e "Into File: ${NEW_DB_FILE}"
set -v
cat ${FASTA_NAMES} > ${NEW_DB_FILE} 
set +v

# Combine FASTA notes files
echo -e "\nCombining Notes Files: ${FASTA_NOTES_NAMES}"
echo -e "Into File: ${NEW_DB_FILE}.notes.txt"
set -v
cat ${FASTA_NOTES_NAMES} > ${NEW_DB_FILE}.notes.txt 
set +v

echo """

  -- ${NEW_DB_FILE} --
    File: '${NEW_DB_FILE}' was created by $0
    (part of the hybkit project, release NA) on $(date).
    by concatenating files: ${FASTA_NAMES} 
""" >> ${NEW_DB_FILE}.notes.txt

# Combine Detail Files
echo -e "\nCombining Detail Files: ${DETAIL_NAMES}"
echo -e "Into File: ${NEW_DB_FILE/.fa/.csv}"
set -v
awk 'FNR==1 && NR!=1{next;}{print}' ${DETAIL_NAMES} > ${NEW_DB_FILE/.fa/.csv}
set +v

# Combine Detail Notes Files
echo -e "\nCombining Detail Notes Files: ${DETAIL_NOTES_NAMES}"
echo -e "Into File: ${NEW_DB_FILE/.fa/.csv}.notes.txt"
set -v
cat ${DETAIL_NOTES_NAMES} > ${NEW_DB_FILE/.fa/.csv}.notes.txt 
set +v

echo """

  -- ${NEW_DB_FILE/.fa/.csv} --
    File: '${NEW_DB_FILE/.fa/.csv}' was created by $0
    (part of the hybkit project, release NA) on $(date).
    by concatenating files: ${DETAIL_NAMES} 
""" >> ${NEW_DB_FILE/.fa/.csv}.notes.txt

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
#echo -e "Into File: ${NEW_DB_FILE}"
#cat ${NODUP_FASTA_NAME_STEMS} > ${NEW_DB_FILE} 

ZIP_FILES="${NEW_DB_FILE}  ${NEW_DB_FILE/.fa/.csv}"
echo "Zipping Files: ${ZIP_FILES}"
gzip --best ${ZIP_FILES}

echo -e "\nDone.\n"
