#!/usr/bin/env bash
# Daniel Stribling
# ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# hkref (hybkit Ref)  http://www.github.com/RenneLab/hkref 
# Hybkit Project : http://www.github.com/RenneLab/hybkit

set -e -u -o pipefail
NOTES="""
Zip final reference file outputs
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
CMB_FA="${NEW_DB_NAME}.fa"

ZIP_FILES="${CMB_FA}  ${CMB_FA/.fa/.csv}"
echo "Zipping Files: ${ZIP_FILES}"
gzip --best ${ZIP_FILES}

echo -e "\nDone.\n"
