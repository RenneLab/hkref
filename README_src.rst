hkref (Hybkit-Reference)
==================================
.. image:: https://img.shields.io/github/v/release/RenneLab/hkref?include_prereleases&logo=github
   :target: https://github.com/RenneLab/hkref/releases
   :alt: GitHub release (latest by date including pre-releases)
.. image:: https://img.shields.io/badge/Made%20with-Python-1f425f.svg?logo=python&logoColor=white
   :target: https://www.python.org/
   :alt: Made with Python
.. image:: https://img.shields.io/badge/Made%20with-Bash-1f425f.svg?logo=gnu-bash
   :target: https://www.gnu.org/software/bash/
   :alt: Made with Bash
.. image:: https://img.shields.io/badge/Bio-python-yellow?logo=data%3Aimage%2Fpng%3Bbase64%2CiVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAMAAAAoLQ9TAAABvFBMVEUAAAD%2F0UEoaaY3cJ%2F%2F0kD2zUWxq2Zvi4VVgJDvyUlIeJfhwk%2FItlyFlnpBdZs%2Bc5w1b6A4cZ%2F%2F%2FwD%2F0UH%2F0UE3cJ83cJ83cJ83cJ%2F%2F0UL%2F0UD%2F0T%2F%2F0UH%2F0kDXvVQ1b6A3cJ83cJ83cJ83cJ%2F%2Fzzr%2F0kT%2F00j%2F0kT%2F0z%2FMuFk7cp02cJ83cJ83cJ83cJ%2F%2F0UH%2F0UH%2F0kX%2F0kb60ELvykiHl3oAT784cZ43cJ83cJ83cJ%2F%2F0UH%2F0UH%2F0UH%2F0T9YgJDQulf%2B0UH51D43cJ83cJ81b6D%2F4TH%2F0UH%2F0UH%2F0UE3cJ8zbqGJmHn%2F0kD%2F0UE3cJ83cJ80bqDSvFX%2F0kD%2F0UE3cJ82b6BghIz%2F00D%2F0UH%2F0UE2cJ80b6CtqWj%2F0z%2F%2F0UE3cJ83cJ83cJ8YYa7%2F0kD%2F0UH%2F0UFxjIRFd5lnh4nfwlA1bp43cJ83cJ83cJ%2F%2F0UH%2F0UH%2F0UH%2F2Tm%2Fsl9Fd5k7cp08dKI7c6E3cJ83cJ%2F%2F0UH%2F0UH%2F0UH%2F0UH5z0NTfpI0b6A6cqE%2FdaM7cqEwa5z%2F0UH%2F0UH%2F0UH%2F0UH%2F0z8vbKM3cJ81b542cJ84caD%2F0UH%2F0UE3cJ83cJ%2F%2F%2F%2F%2BMZpkLAAAAk3RSTlMAAAAAAAAAAAAAAAAAAAAAAAAAAQUKTlYRBiYNBmYkifDyog4EdORih83qcl7scBzb9lBH%2BrIKA5bYH17tdwxW%2Fa0GMOV0DcOtCAu08Og8AprNa%2BtLSOtuzZsDPOjutwwHqsUOcuYyCbr%2BWQtz7WIe1Z4q3PhHTfTdHWXuhrfVsolh5nkECYTo0UQcaQ0pBiYaBQETquGoAAAAl0lEQVQY02NgIAQYGRiFGN8JMzIyvmb4yygBpBkYpRkh4ME%2FZSDJwsDwQgUiwMDCDRb4z8gNETA%2FwwMWZ2axAFLHrCGiO4BanMEssCzjWgYmZgZuEAjcDaYY%2FrMwxEItASlZxMDAlMHDw7MKKMcJUsD0n4EJTE%2Fi4soAM5ghAl8g5nBzBzCABXr%2F%2Fq%2BZxsFR0cFmxEAEAAAfKxn6VT4rZAAAAABJRU5ErkJggg%3D%3D%0A
   :target: https://github.com/biopython/biopython
   :alt: BioPython Project
.. image:: http://biothings.io/static/img/powered-by-mygene.png
   :width: 128 px
   :target: https://mygene.info/  
   :alt: Powered by MyGene.info

| This repository is a part of the `hybkit project <http://www.github.com/RenneLab/hybkit>`_.
| Full hybkit project documentation is available at
  `hybkit's ReadTheDocs <https://hybkit.readthedocs.io/>`_.

Description:
============
| This repository includes an up-to-date human genomic sequence reference designed to be 
  compaitble with the `Hyb <https://github.com/gkudla/hyb>`_ program 
  for chimeric (hybrid) read calling for ribonomics experiments.
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
  * mRNA : transcript_biotype=protein_coding; cdna; 
    (limited to where a RefSeq Protein Identifier Exists)
  * lncRNA : transcript_biotype=lncrna; cdna
  * other : transcript_biotype=[all remaining, excluding miRNA]; transcript_exon_intron
  
For a detailed description of the current sequences and queries utilized, see 
"Current Reference Details" below.

Run Reference Creation Pipeline:
================================
The source code utilized to create the current version of the reference is provided here,
and has been utilized on a unix system (bash shell). 
This code depends on `SeqKit <https://bioinf.shenwei.me/seqkit/>`_ 
for sequence maniuplation, Python 3.6+, and multiple python packages:

Required Python Packages:
  * `pandas <https://pandas.pydata.org/>`_
  * `natsort <https://pypi.org/project/natsort/>`_
  * `pybiomart <https://pypi.org/project/pybiomart/>`_
  * `biothings-client <https://pypi.org/project/biothings-client/>`_
  * `biopython <https://biopython.org/>`_

The scripts can be run by executing the first script: "00_run_all_steps.sh" with all 
required resources (seqkit, python3) available on the system path.

Hyb Reference Specification:
============================
The Hyb program has requirements about the formatting of the FASTA file used for the reference.

Currently identified requirements include:
  * No description in FASTA sequence header (no whitespace characters)
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

Text of: *./_REF_VERSION.sh*

.. include:: ./_REF_VERSION.sh
   :code: bash
   :start-line: 7

Current Reference Details:
==========================

Text of: *./01_notes.sh*

.. include:: ./01_notes.sh
   :code: bash
   :start-line: 8
   :end-line: 37
 
