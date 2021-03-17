#!/usr/bin/env bash
# Daniel Stribling
# ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# hkref (hybkit Ref)  http://www.github.com/RenneLab/hkref 
# Hybkit Project : http://www.github.com/RenneLab/hybkit

set -e -u -o pipefail 
NOTES="""
Align miRBase to downloaded ensembl sequences.
"""
if [ ! -z ${2+x} ]; then
  USE_MODULES=$2
else
  USE_MODULES=0
fi

MIRNA_FA="$(find hyb_miRBase*.fa)"
SEQ_FA="$(find *_transcripts_only.fa)"
SEQ_FAI="${SEQ_FA}.fai"
SEQ_GENOME="${SEQ_FA/.fa/.genome}"
ALN_DIR="miRBase_align"
ALN_SAM="${ALN_DIR}/$(basename ${MIRNA_FA/.fa/.sam})"
ALN_BED="${ALN_DIR}/$(basename ${MIRNA_FA/.fa/.bed})"
ALN_BED_TRIM="${ALN_DIR}/$(basename ${MIRNA_FA/.fa/_trim.bed})"
COMPLEMENT_BED="${ALN_DIR}/$(basename ${MIRNA_FA/.fa/_complement.bed})"
OUT_TRANS_FA="${SEQ_FA/_transcripts_only.fa/_transcripts_only_ext.fa}"
OUT_TRANS_FA_NOCOORD="${OUT_TRANS_FA/_ext.fa/_ext_nocoord.fa}"
REF_DIR="${ALN_DIR}/ensembl_ref"
ALL_COMMANDS=""
ALL_VERSIONS=""
mkdir -pv ${ALN_DIR}
mkdir -pv ${REF_DIR}

CPUS=1

if [ $USE_MODULES -gt 0 ]; then
  CPUS=8
  module load bowtie2/2.4.2
  ALL_VERSIONS="bowtie2/2.4.2"
fi

echo -e "\nAligning miRNAs: ${MIRNA_FA}"
echo -e "  to Transcripts: ${SEQ_FA}"
echo -e "  to find overlap."

# -- Align miRNA To Transcripts
DB_NAME="${REF_DIR}/$(basename ${SEQ_FA/.fa/})"

set -v -H -o history
bowtie2-build --threads ${CPUS} ${SEQ_FA} ${DB_NAME}
COMMAND="!!" ; ALL_COMMANDS="${ALL_COMMANDS}\n${COMMAND}"
bowtie2 --no-unal -x ${DB_NAME} -f ${MIRNA_FA} --threads ${CPUS} -S ${ALN_SAM} 2>&1
COMMAND="!!" ; ALL_COMMANDS="${ALL_COMMANDS}\n${COMMAND}"
set +v +H +o history

echo -e "\nProcessing alignments into BED file and inverting to represent "
echo -e "non-aligned regions"

# -- Convert Alignments to BED file
if [ $USE_MODULES -gt 0 ]; then
  module load bedops/2.4.30 samtools/1.10 seqkit/0.10.2
  ALL_VERSIONS="${ALL_VERSIONS}, bedops/2.4.30, samtools/1.10, seqkit/0.10.2"
fi
set -v -H -o history
sam2bed < ${ALN_SAM} | cut -f1-6> ${ALN_BED}
COMMAND="!!" ; ALL_COMMANDS="${ALL_COMMANDS}\n${COMMAND}"

# -- Create "genome" file reference for transcripts
samtools faidx ${SEQ_FA}
COMMAND="!!" ; ALL_COMMANDS="${ALL_COMMANDS}\n${COMMAND}"

awk -v OFS='\t' {'print $1,$2'} ${SEQ_FAI} > ${SEQ_GENOME}
COMMAND="!!" ; ALL_COMMANDS="${ALL_COMMANDS}\n${COMMAND}"
set +v +H +o history

# -- Reduce size of BED features to remove middle of transcript
if [ $USE_MODULES -gt 0 ]; then
  module load bedtools/2.30.0 seqkit/0.10.2
  ALL_VERSIONS="${ALL_VERSIONS}, bedtools/2.30.0"
fi
set -v -H -o history
bedtools slop -i ${ALN_BED} -g ${SEQ_GENOME} -b "-6" > ${ALN_BED_TRIM} #/.bed/_unsort.bed}
COMMAND="!!" ; ALL_COMMANDS="${ALL_COMMANDS}\n${COMMAND}"

#bedtools sort -i ${ALN_BED_TRIM/.bed/_unsort.bed} > ${ALN_BED_TRIM}
#COMMAND="!!" ; ALL_COMMANDS="${ALL_COMMANDS}\n${COMMAND}"

# -- Get sequences with middle of identified mature miRNAs extracted
bedtools complement -i ${ALN_BED_TRIM} -g ${SEQ_GENOME} > ${COMPLEMENT_BED}
COMMAND="!!" ; ALL_COMMANDS="${ALL_COMMANDS}\n${COMMAND}"

set -v -H -o history
echo -e "\nExtracting non-miRNA regions of transcripts from FASTA."
# -- Extract sequences from fasta.
bedtools getfasta -fi ${SEQ_FA} -bed ${COMPLEMENT_BED} | \
  seqkit seq --min-len 15 --upper-case --validate-seq > ${OUT_TRANS_FA}
COMMAND="!!" ; ALL_COMMANDS="${ALL_COMMANDS}\n${COMMAND}"
set +v +H +o history

cp -v ${SEQ_FA}.notes.txt ${OUT_TRANS_FA}.notes.txt
echo -e """

  -- ${OUT_TRANS_FA} --
    File: '${OUT_TRANS_FA}' was created by $0
    (part of the hybkit project, release NA) on $(date).
    using ${ALL_VERSIONS}
    by commands: 
    ${ALL_COMMANDS} 
""" >> ${OUT_TRANS_FA}.notes.txt

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

#ZIP_FILES="${NEW_DB_FILE}  ${NEW_DB_FILE/.fa/.csv}"
#echo "Zipping Files: ${ZIP_FILES}"
#gzip --best ${ZIP_FILES}

echo -e "\nDone.\n"
