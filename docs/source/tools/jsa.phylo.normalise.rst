-------------------------------------------------------------
*jsa.phylo.normalise*: Normalise branch length of a phylogeny 
-------------------------------------------------------------

*jsa.phylo.normalise* scales the branches of a phylogeny so that their sum
equates to a value.
 
*jsa.phylo.normalise* is included in the `Japsa package <http://mdcao.github.io/japsa/>`_. 
Please see check the installation_ page for instructions.  

.. _installation: ../install.html

~~~~~~~~
Synopsis
~~~~~~~~

*jsa.phylo.normalise*: Scale branches of a phylogeny so that the sum of branch lengths is equal to a value

~~~~~
Usage
~~~~~
::

   jsa.phylo.normalise [options]

~~~~~~~
Options
~~~~~~~
  --input=s       Name of the input file, - for standard input
                  (REQUIRED)
  --sum=d         Sum of branches after normalising
                  (default='1.0')
  --scale=d       Scale factor, if set to a positive number will override the sum parameter
                  (default='0.0')
  --output=s      Name of the file for output, - for stdout
                  (default='-')
  --help          Display this usage and exit
                  (default='false')





