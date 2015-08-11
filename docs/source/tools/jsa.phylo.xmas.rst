----------------------------------------------------------------------
*XMas*: Robust estimation of genetic distances with information theory
----------------------------------------------------------------------

*XMas* (jsa.phylo.xmas) is a tool to measure genetics distances between
aligned sequences. It reads in a list of sequences from a fasta file format,
and outputs a distance matrix in the format required by the PHYLIP package
to run neighbour joining.
 
*XMas* is included in the `Japsa package <http://mdcao.github.io/japsa/>`_. 
Please see check the installation_ page for instructions.  

.. _installation: ../install.html

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




-------------
Usage samples
-------------

At the moment, XMas is designed to worked with aligned sequences, with indels 
and wildcards (*e.g.*, N) removed. XMas reads in these aligned sequences from
a fasta file, and output the distances to a file in a format ready to run
neibour-joining with PHYLIP. For examples::
   	
   	jsa.phylo.xmas -i sequences.fas -o infile

And run phylip neighbor-joining from the distances in *infile*::

   phylip neighbor
   
   	 




