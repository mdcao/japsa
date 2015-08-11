--------------------------------------------
*jsa.seq.sort*: Sort the sequences in a file
--------------------------------------------

*jsa.seq.sort* sort the sequences from a file or from a standard input into
some order.

*jsa.seq.sort* is included in the 
`Japsa package <http://mdcao.github.io/japsa/>`_. 
Please see check the installation_ page for instructions.  

.. _installation: ../install.html

~~~~~~~~
Synopsis
~~~~~~~~

*jsa.seq.sort*:Sort sequences based on their lengths

~~~~~
Usage
~~~~~
::

   jsa.seq.sort [options]

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
  --number        Add the order number to the beginning of contig name
                  (default='false')
  --reverse       Reverse sort order
                  (default='false')
  --sortKey=s     Sort key
                  (default='length')
  --help          Display this usage and exit
                  (default='false')




