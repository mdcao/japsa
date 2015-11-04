----------------------------------------------------------------------
*jsa.np.rtStrainTyping*: Bacterial strain typing with Oxford Nanopore sequencing
----------------------------------------------------------------------

*jsa.np.rtStrainTyping* strain types a bacterial sample from Oxford Nanopore
sequencing in real-time. It is included in the 
`Japsa package <http://mdcao.github.io/japsa/>`_. 
Please see check the installation_ page for instructions.  

.. _installation: ../install.html

~~~~~~~~
Synopsis
~~~~~~~~

*jsa.np.rtStrainTyping*:Realtime strain typing using Nanopore sequencing data

~~~~~
Usage
~~~~~
::

   jsa.np.rtStrainTyping [options]

~~~~~~~
Options
~~~~~~~
  --geneDB=s       Path to the gene database
                  (REQUIRED)
  --bamFile=s     The bam file
                  (REQUIRED)
  --qual=d        Minimum alignment quality
                  (default='0.0')
  --twodonly      Use only two dimentional reads
                  (default='false')
  --read=i        Minimum number of reads between analyses
                  (default='50')
  --time=i        Minimum number of seconds between analyses
                  (default='30')
  --output=s      Output file
                  (default='output.dat')
  --help          Display this usage and exit
                  (default='false')




