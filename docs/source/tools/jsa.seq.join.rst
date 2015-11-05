-----------------------------------------------------
*jsa.seq.join*: Join multiple sequences into one file 
-----------------------------------------------------


~~~~~~~~
Synopsis
~~~~~~~~

*jsa.seq.join*: Join multiple sequences into one

~~~~~
Usage
~~~~~
::

   jsa.seq.join [options] file1 file2 ...

~~~~~~~
Options
~~~~~~~
  --alphabet=s    Alphabet of the input file. Options: DNA (DNA=DNA16), DNA4
                  (ACGT), DNA5(ACGTN), DNA16 and Protein
                  (default='DNA')
  --output=s      Name of the output file, - for standard output
                  (default='-')
  --name=s        Name of the new sequence
                  (default='newseq')
  --removeN       Remove wildcards (N)
                  (default='false')
  --help          Display this usage and exit
                  (default='false')




