----------------------------------------------------------
*jsa.seq.gff2fasta*: Join multiple sequences into one file 
----------------------------------------------------------

*jsa.seq.gff2fasta* extract the functional sequences (genes, CDS, etc) from
a gff file and a sequence file.

~~~~~~~~
Synopsis
~~~~~~~~

*jsa.seq.gff2fasta*: Extract sequences from a gff annotation

~~~~~
Usage
~~~~~
::

   jsa.seq.gff2fasta [options]

~~~~~~~
Options
~~~~~~~
  --sequence=s    The sequence (whole chromosome)
                  (REQUIRED)
  --gff=s         Annotation file in gff format
                  (REQUIRED)
  --type=s        types of features to be extracted (all, gene, CDS etc)
                  (default='gene')
  --output=s      Name of the output file, - for standard output
                  (REQUIRED)
  --alphabet=s    Alphabet of the input file. Options: DNA (DNA=DNA16), DNA4
                  (ACGT), DNA5(ACGTN), DNA16 and Protein
                  (default='DNA')
  --help          Display this usage and exit
                  (default='false')




