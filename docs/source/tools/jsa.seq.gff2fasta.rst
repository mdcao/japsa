-------------------------------------------
*jsa.seq.gff2fasta*: Extract gene sequences 
-------------------------------------------

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
  --flank=i       Size of flanking regions
                  (default='0')
  --output=s      Name of the output file, - for standard output
                  (REQUIRED)
  --help          Display this usage and exit
                  (default='false')




