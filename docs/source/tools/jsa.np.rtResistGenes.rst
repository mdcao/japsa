-------------------------------------------------------------------------------------------------------
*jsa.np.rtResistGenes*: Antibiotic resistance gene identification in real-time with Nanopore sequencing 
-------------------------------------------------------------------------------------------------------

*jsa.np.rtResistGenes* identifies antibiotic resistance genes from real-time sequencing
with Nanopore MinION. 

~~~~~~~~
Synopsis
~~~~~~~~

*jsa.np.rtResistGenes*: Realtime identification of antibiotic resistance genes from Nanopore sequencing

~~~~~
Usage
~~~~~
::

   jsa.np.rtResistGenes [options]

~~~~~~~
Options
~~~~~~~
  --output=s      Output file
                  (default='output.dat')
  --bamFile=s     The bam file
                  (default='null')
  --score=d       The alignment score threshold
                  (default='1.0E-4')
  --msa=s         Name of the msa method, support poa, kalign, muscle and clustalo
                  (default='kalign')
  --tmp=s         Temporary folder
                  (default='\_tmpt')
  --resDB=s       Path to resistance database
                  (REQUIRED)
  --qual=d        Minimum alignment quality
                  (default='0.0')
  --twodonly      Use only two dimentional reads
                  (default='false')
  --read=i        Minimum number of reads between analyses
                  (default='50')
  --time=i        Minimum number of seconds between analyses
                  (default='1800')
  --thread=i      Number of threads to run
                  (default='4')
  --help          Display this usage and exit
                  (default='false')


~~~~~~~~
See also
~~~~~~~~

jsa.np.f5reader_, jsa.np.rtSpeciesTyping_, jsa.np.rtStrainTyping_, jsa.util.streamServer_, jsa.util.streamClient_

.. _jsa.np.f5reader: jsa.np.f5reader.html
.. _jsa.np.rtSpeciesTyping: jsa.np.rtSpeciesTyping.html
.. _jsa.np.rtStrainTyping: jsa.np.rtStrainTyping.html
.. _jsa.util.streamServer: jsa.util.streamServer.html
.. _jsa.util.streamClient: jsa.util.streamClient.html



~~~~~~~~~~
Setting up
~~~~~~~~~~

Refer to the documentation at https://github.com/mdcao/npAnalysis/ for more 
details.


