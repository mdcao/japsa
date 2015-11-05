----------------------------------------------------------------------------------
*jsa.np.rtSpeciesTyping*: Bacterial species typing with Oxford Nanopore sequencing
----------------------------------------------------------------------------------

*jsa.np.rtSpeciesTyping* identify propotions of species from a DNA sample 
using Oxford Nanopore sequencing in real-time. It reads data in SAM/BAM format
of the alignments of sequence reads to a collection of species genomes.

We provide a genome collection of nearly 1500 bacterial species http://genomicsresearch.org/public/researcher/npAnalysis/SpeciesTyping.tar.gz 
(or https://swift.rc.nectar.org.au:8888/v1/AUTH_15574c7fb24c44b3b34069185efba190/npAnalysis/SpeciesTyping.tar.gz).
Obtain them by::
   wget https://swift.rc.nectar.org.au:8888/v1/AUTH_15574c7fb24c44b3b34069185efba190/npAnalysis/SpeciesTyping.tar.gz.tar.gz
   tar zxvf SpeciesTyping.tar.gz

which will generate three folders for the three species.

~~~~~~~~
Synopsis
~~~~~~~~

*jsa.np.rtSpeciesTyping*: Realtime species typing using Nanopore Sequencing data

~~~~~
Usage
~~~~~
::

   jsa.np.rtSpeciesTyping [options]

~~~~~~~
Options
~~~~~~~
  --output=s      Output file
                  (default='output.dat')
  --bamFile=s     The bam file
                  (REQUIRED)
  --indexFile=s   indexFile 
                  (REQUIRED)
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

jsa.np.f5reader_, jsa.np.rtStrainTyping_, jsa.np.rtResistGenes, jsa.util.streamServer_, jsa.util.streamClient_

.. _jsa.np.f5reader: jsa.np.f5reader.html
.. _jsa.np.rtStrainTyping: jsa.np.rtStrainTyping.html
.. _jsa.util.streamServer: jsa.util.streamServer.html
.. _jsa.util.streamClient: jsa.util.streamClient.html



~~~~~~~~~~~~~~
Usage examples
~~~~~~~~~~~~~~

If there is a sam/bam file of aligning the Nanopore sequencing to the genome 
collection, the program can read from this
::
   jsa.np.rtSpeciesTyping -bam alignment.sam -index SpeciesTyping/Bacteria/speciesIndex --read 50 -time 60 -out speciesTypingResults.out
   
   
This program can read data from the output stream of an alignment program to
perform analysis in real-time. For example, one can create such a pipeline
to listen on port 3456
::
  jsa.util.streamServer -port 3456 \
  | bwa mem -t 10 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -Y -K 10000 SpeciesTyping/Bacteria/genomeDB.fasta - 2> /dev/null \
  | jsa.np.rtSpeciesTyping -bam - -index SpeciesTyping/Bacteria/speciesIndex --read 50 -time 60 -out speciesTypingResults.out 2>  speciesTypingResults.log &
  
  
and streams data to this pipeline using npReader:
::
  jsa.np.f5reader -GUI -realtime -folder <DownloadFolder> -fail -output data.fastq -stream serverAddress:3456


