*******************************************************
*                                                     *
*        I  N  F  O  R  M  A  T  I  V  I  T  Y        *
*                                                     *
*                                                     *
*******************************************************

Developed for Windows and UNIX.

0. Prerequisites
Requires blast+ be installed and PATH variable set.
BLAST executables are available: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST
Refer to NCBI for instructions/troubleshooting install and setting the path.
Please note, this code was tested for ncbi-blast-2.4.0. No guarantee using prior versions of blast.


1. Download/Install
Download version for your OS and build.
Build:
g++ -o informativity informativity-unix.cpp


2. Running Informativity
The code requires 4 (or 5) parameters:
- Nucleotide sequences for coding regions of the genome of interest (X)
- Amino acid sequences for coding regions of the genome of interest (X)
- Amino acid sequences for coding regions of the related genome for thresholding (G)
- Amino acid sequences for coding regions of outgroup
- Contigs (nucleotide sequences) for metagenomic dataset *optional*

Example, execution via UNIX:
./informativity NC_011810.fna NC_011810.faa NC_009015.faa all-non-Pbuna.fasta sample_metagenome2.fasta

The user will be prompted for a name for the run. A file will be created in which all temporary and result files will be written.


3. Sample data provided
Sample data files for identifying informative genes for Pseudomonas phage PB1 (NC_011810) is provided. The outgroup (Burkholderia ambifaria phage BcepF1; NC_009015) files are included. Sequences for phages not classified as Pbunalikeviruses (file name: "all-non-Pbuna.fasta"), a subset of those available through GenBank, are also included in FASTA in a gz file that needs to be unpacked prior to use. To test the identification of PB1 genes within a metagenome data set, we have also included "sample_metagenome2.fasta". This is a small dataset from sequencing of a PB1-like virus. Sample data is available as a compressed folder through github.


Questions/Problems: Email cputonti@luc.edu.