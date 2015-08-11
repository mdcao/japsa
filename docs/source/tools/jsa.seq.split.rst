---------------------------------------------
*jsa.seq.split*: Split multiple sequence file
---------------------------------------------

*jsa.seq.split* splits a file containing multiple sequences to each file
containing a sequence. It is included in the 
`Japsa package <http://mdcao.github.io/japsa/>`_. 
Please see check the installation_ page for instructions.  

.. _installation: ../install.html

~~~~~~~~
Synopsis
~~~~~~~~

*jsa.seq.split*:Break a multiple sequence files to each sequence per file

~~~~~
Usage
~~~~~
::

   jsa.seq.split [options]

~~~~~~~
Options
~~~~~~~
  --input=s       Name of the input file, - for standard input
                  (REQUIRED)
  --alphabet=s    Alphabet of the input file. Options: DNA (DNA=DNA16), DNA4
                  (ACGT), DNA5(ACGTN), DNA16 and Protein
                  (default='DNA')
  --output=s      Prefix of the output files
                  (default='out\_')
  --format=s      Format of output files. Options : japsa or fasta
                  (default='fasta')
  --help          Display this usage and exit
                  (default='false')




