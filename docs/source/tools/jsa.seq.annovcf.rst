--------------------------------------
*jsa.seq.annovcf*: Annotate a vcf file 
--------------------------------------

*jsa.seq.annovcf*: reads annotations from a gff file and annotates a vcf file 
(i.e., identify the functional of variation). 

*jsa.seq.annovcf* is included in the 
`Japsa package <http://mdcao.github.io/japsa/>`_. 
Please see check the installation_ page for instructions.  

.. _installation: ../install.html

~~~~~~~~
Synopsis
~~~~~~~~

*jsa.seq.annovcf*: Annotate variation in a vcf file using annotation from gff file

~~~~~
Usage
~~~~~
::

   jsa.seq.annovcf [options]

~~~~~~~
Options
~~~~~~~
  --gffin=s       GFF file
                  (default='-')
  --upstream=i    Add upstream 
                  (default='0')
  --downstream=i  Add downstream 
                  (default='0')
  --output=s      Name of output file, - for standard out
                  (default='-')
  --vcf=s         Name of vcf file
                  (REQUIRED)
  --help          Display this usage and exit
                  (default='false')




