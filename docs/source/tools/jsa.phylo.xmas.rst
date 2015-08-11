----------------------------------------------------------------------
*XMas*: Robust estimation of genetic distances with information theory
----------------------------------------------------------------------

*XMas* (jsa.phylo.xmas) is a tool to measure genetics distances between
aligned sequences. It reads in a list of sequences from a fasta file format,
and outputs a distance matrix in the format required by the PHYLIP package
to run neighbour joining.
 
*XMas* is included in the `Japsa package <http://mdcao.github.io/japsa/>`_. 
Please see check the installation page for instructions.  


~~~~~~~~
Synopsis
~~~~~~~~

*jsa.phylo.xmas*:Generate a distance matrix from aligned sequences

~~~~~
Usage
~~~~~
::

   jsa.phylo.xmas [options]

~~~~~~~
Options
~~~~~~~
  --input=s       Name of the input file, - for standard input
                  (REQUIRED)
  --output=s      Name of the file for output (distances in phylip format)
                  (default='output')
  --adapt         Use adaptive
                  (default='false')
  --help          Display this usage and exit
                  (default='false')






