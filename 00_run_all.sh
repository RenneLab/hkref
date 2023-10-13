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

# Load local resources if present
if [ -f "./load_resources.source.sh" ]; then
    source ./load_resources.source.sh
fi

bash 01_notes.sh

nextflow run 02_hkref_build.nf -resume

echo -e "\nDone at $(date)\n"


