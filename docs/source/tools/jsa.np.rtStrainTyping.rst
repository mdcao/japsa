--------------------------------------------------------------------------------
*jsa.np.rtStrainTyping*: Bacterial strain typing with Oxford Nanopore sequencing
--------------------------------------------------------------------------------

*jsa.np.rtStrainTyping* strain types a bacterial sample from Oxford Nanopore
sequencing in real-time. It reads data in SAM/BAM format from a file or from
a stream and identifies the genes in the samples. Based on the patterns of 
gene presence, it makes an inference of the strain, together with the confidence
interval of 95%.

We provide the gene databases for three bacterial species  K. pneumoniae, 
E. coli and S. aureus on  http://genomicsresearch.org/public/researcher/npAnalysis/StrainTyping.tar.gz 
(or https://swift.rc.nectar.org.au:8888/v1/AUTH_15574c7fb24c44b3b34069185efba190/npAnalysis/StrainTyping.tar.gz).
Obtain them by::
   wget https://swift.rc.nectar.org.au:8888/v1/AUTH_15574c7fb24c44b3b34069185efba190/npAnalysis/StrainTyping.tar.gz
   tar zxvf StrainTyping.tar.gz

which will generate three folders for the three species.

~~~~~~~~
Synopsis
~~~~~~~~

*jsa.np.rtStrainTyping*: Realtime strain typing using Nanopore sequencing data

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


~~~~~~~~
See also
~~~~~~~~

jsa.np.f5reader_, jsa.np.rtSpeciesTyping, jsa.np.rtResistGenes, jsa.util.streamServer_, jsa.util.streamClient_

.. _jsa.np.f5reader: jsa.np.f5reader.html
.. _jsa.util.streamServer: jsa.util.streamServer.html
.. _jsa.util.streamClient: jsa.util.streamClient.html



~~~~~~~~~~~~~~
Usage examples
~~~~~~~~~~~~~~
If there is a sam/bam file of aligning the Nanopore sequencing to the gene 
database (ie, geneFam.fasta in one of the said databases), the program
can read from this file (note, this is not real-time analysis):
::
   jsa.np.rtStrainTyping -geneDB StrainTyping/Escherichia_coli/ -bamFile ecoli.bam -read 100000 -time 1000 -output output.dat
   
This program can read data from the output stream of an alignment program to
perform analysis in real-time. For example, one can create such a pipeline
to listen on port 3457
::
  jsa.util.streamServer -port 3457 \
  | bwa mem -t 2 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -Y -K 10000 -a StrainTyping/Escherichia_coli/geneFam.fasta - 2> /dev/null \
  | jsa.np.rtStrainTyping -bam - -geneDB StrainTyping/Escherichia_coli/ -read 0 -time 20 --out EcStrainTyping.dat 2>  kPStrainTyping.log
  
and streams data to this pipeline using npReader:
::
  jsa.np.f5reader -GUI -realtime -folder <DownloadFolder> -fail -output data.fastq -stream serverAddress:3457

