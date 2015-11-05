-----------------------------------------------------
*jsa.seq.extract*: Extract subsequences from a genome 
-----------------------------------------------------


*jsa.seq.extract* is included in the 
`Japsa package <http://mdcao.github.io/japsa/>`_. 
Please see check the installation_ page for instructions.  

.. _installation: ../install.html

~~~~~~~~
Synopsis
~~~~~~~~

*jsa.seq.extract*: Extract subsequences

~~~~~
Usage
~~~~~
::

   jsa.seq.extract [options] <chr:start-end> <chr:start-end> ...

~~~~~~~
Options
~~~~~~~
  --input=s       Name of the input file, - for standard input
                  (REQUIRED)
  --output=s      Name of the output file, - for standard output
                  (REQUIRED)
  --alphabet=s    Alphabet of the input file. Options: DNA (DNA=DNA16), DNA4
                  (ACGT), DNA5(ACGTN), DNA16 and Protein
                  (default='DNA')
  --reverse       Reverse complement the subsequence
                  (default='false')
  --format=s      format of the output file (jsa and fasta)
                  (default='fasta')
  --help          Display this usage and exit
                  (default='false')




