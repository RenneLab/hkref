
hkref (Hybkit-Reference)
==================================
.. image:: https://img.shields.io/github/v/release/RenneLab/hkref?include_prereleases
   :target: https://github.com/RenneLab/hkref/releases
   :alt: GitHub release (latest by date including pre-releases)


| This repository is a part of the hybkit project: available at
  http://www.github.com/RenneLab/hybkit .
| Full hybkit project documentation is available at
  `hybkit's ReadTheDocs <https://hybkit.readthedocs.io/>`_.

Description:
============
| This repository includes an up-to-date human genomic sequence reference designed to be 
  compaitble with the Hyb program for chimeric (hybrid) read calling for ribonomics experiments,
  data, available at https://github.com/gkudla/hyb . 
| The method for reference library construction is based on the protocol provided in the 
  supplemental methods of:

  | Helwak, Aleksandra, et al. 'Mapping the human miRNA interactome by CLASH reveals
    frequent noncanonical binding.' Cell 153.3 (2013): 654-665.
    http://dx.doi.org/10.1016/j.cell.2013.03.043

| One significant departure from the original method used to create the hOH7 database included
  with Hyb is that duplicate sequences and subsequences are allowed within the reference.
  This was chosen to capture the full potential variation in transcript isoform that may
  be present within a sequence dataset.   
  The reference library is primarily based on sequences downloaded from Ensembl via the
  Biomart API, with the use of miRBase for mature miRNA sequences and a few other sequence
  sources.
  
Biomart queries include:
  * mRNA : transcript_biotype=protein_coding; cDNA; 
    (limited to where a RefSeq Protein Identifier Exists)
  * lncRNA (as cDNA)
  * All remaining Ensembl transcript_biotype values, (excluding miRNA).
  
For a detailed description of the current sequences and queries utilized, see 
"Current Reference Details" below.

Run Reference Creation Pipeline:
================================
The source code utilized to create the current version of the reference is provided here,
and has been utilized on a unix system (bash shell). 
This code depends on SeqKit ( https://bioinf.shenwei.me/seqkit/ )
for sequence maniuplation, Python 3.6+, and well as multiple python packages:

Required Python Packages:
  * pandas 
  * natsort
  * pybiomart
  * biothings-client
  * biopython

The scripts can be run by executing the first script: "00_run_all_steps.sh" with all 
required resources (seqkit, python) available on the system path.

Hyb Reference Specification:
============================
The Hyb program has requirements about the formatting of the fasta file used for the reference.

Currently identified requirements include:
  * No description in header (no whitespace characters)
  * | Sequence identifier be of the form of "{1}_{2}_{seqid}_{biotype}"
    | {1}: Arbitrary Identifier (ENSG... for Ensembl Sequences)
    | {2}: Arbitrary Identifier (ENST... for Ensembl Sequences)
    | {3}: Name of gene/miRNA
    | {4}: Ensembl-style transcript_biotype. 
    |     (Note, "microRNA" must be used in place of "miRNA" for recognition by Hyb)

  * {1}, {2}, {seqid}, and {biotype} should contain only [a-z], [A-Z], [0-9],
     "-", and "|" characters.
  * "_", ".", and "," characters are *specifically* excluded from identifiers.

Examples:
      
.. code-block:: bash  

    >ENSG00000003137_ENST00000001146_CYP26B1_mRNA
    .....
    >MIMAT0000062_MirBase_let-7a_microRNA
    TGAGGTAGTAGGTTGTATAGTT

Thanks to Grzegorz Kudla ( https://github.com/gkudla ) for providing information on
Hyb reference creation. 

Current Reference Version:
==========================

.. include:: ./_REF_VERSION.sh
   :code: bash
   :start-line: 8
   :caption: _REF_VERSION.sh (Current Reference Version)
   :name: version_file
 

Current Reference Details:
==========================

.. include:: ./01_notes.sh
   :literal:
   :start-line: 9
   :end-line: 37
   :caption: 01_notes.sh (Current Reference Details)
   :name: notes_file
 
