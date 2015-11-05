------------------------------------------------------------------------------------
*jsa.np.rtMLST*: Multi-locus Sequencing Typing in real-time with Nanopore sequencing 
------------------------------------------------------------------------------------

*jsa.np.rtMLST* performs MLST typing from real-time sequencing with Nanopore MinION. 

~~~~~~~~
Synopsis
~~~~~~~~

*jsa.np.rtMLST*: Realtime Multi-Locus Strain Typing using Nanopore Sequencing data

~~~~~
Usage
~~~~~
::

   jsa.np.rtMLST [options]

~~~~~~~
Options
~~~~~~~
  --output=s      Output file
                  (default='output.dat')
  --mlstScheme=s  Path to mlst scheme
                  (REQUIRED)
  --bamFile=s     The bam file
                  (default='null')
  --qual=d        Minimum alignment quality
                  (default='0.0')
  --twodonly      Use only two dimentional reads
                  (default='false')
  --read=i        Minimum number of reads between analyses
                  (default='50')
  --time=i        Minimum number of seconds between analyses
                  (default='30')
  --help          Display this usage and exit
                  (default='false')


~~~~~~~~
See also
~~~~~~~~

jsa.np.f5reader_, jsa.np.rtSpeciesTyping_, jsa.np.rtStrainTyping_, jsa.np.rtResistGenes_, jsa.util.streamServer_, jsa.util.streamClient_

.. _jsa.np.f5reader: jsa.np.f5reader.html
.. _jsa.np.rtSpeciesTyping: jsa.np.rtSpeciesTyping.html
.. _jsa.np.rtStrainTyping: jsa.np.rtStrainTyping.html
.. _jsa.np.rtResistGenes: jsa.np.rtResistGenes.html
.. _jsa.util.streamServer: jsa.util.streamServer.html
.. _jsa.util.streamClient: jsa.util.streamClient.html



~~~~~~~~~~
Setting up
~~~~~~~~~~

Refer to real-time analyais page at https://github.com/mdcao/npAnalysis/

