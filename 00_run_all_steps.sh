#!/usr/bin/env bash
# Daniel Stribling
# ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# hkref (hybkit Ref)  http://www.github.com/RenneLab/hkref 
# Hybkit Project : http://www.github.com/RenneLab/hybkit

set -e -u -o pipefail

NOTES="""
See 01_notes.sh for detailed pipeline information.
"""

# Check if email passed on command line
if [ ! -z ${1+x} ]; then
  USER_EMAIL="${1}"
else
  echo "Please enter your email address:"
  read USER_EMAIL
  #USER_EMAIL=""
fi
if [ ! -z ${2+x} ]; then
  USE_MODULES="${2}"
else
  USE_MODULES=0
fi

# Databse version information

# Load NEW_DB_VER and NEW_DB_NAME variables:
source _REF_VERSION.sh

LOG="${NEW_DB_NAME}_log.txt"
RUN_DIR="work_dir"

# Get names of scripts.
SCRIPT_NAMES=$(sort --version-sort <<< "$(find 0?_*.sh 0?_*.py 1?_*.sh)" )

# Perform work in a separate directory
mkdir -pv ${RUN_DIR}
cd ${RUN_DIR}
echo "Starting on $(date)" > ${LOG} 

# Run each script 
for script in ${SCRIPT_NAMES} ; do
  if [ "${script}" != "$(basename $0)" ] ; then
    ../$script ${USER_EMAIL} ${USE_MODULES} | tee ${LOG}
  fi
done

echo -e "\nDone at $(date)\n" | tee ${LOG}
# Return output to main directory. 
mv -v ${NEW_DB_NAME}.* ../ 
mv -v ${LOG} ../ 

