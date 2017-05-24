----------------------------------------------------------------------------
*capsim*: Simulating the Dynamics of Targeted Capture Sequencing with CapSim
----------------------------------------------------------------------------

*capsim* (jsa.sim.capsim) is a tool to simulate target capture sequencing.

~~~~~~~~
Synopsis
~~~~~~~~

*jsa.sim.capsim*: Simulate capture sequencing

~~~~~
Usage
~~~~~
::

   jsa.sim.capsim [options]

~~~~~~~
Options
~~~~~~~
  --reference=s   Name of genome to be 
                  (REQUIRED)
  --probe=s       File containing probes mapped to the reference in bam format
                  (default='null')
  --logFile=s     Log file
                  (default='-')
  --ID=s          A unique ID for the data set
                  (default='')
  --miseq=s       Name of read file if miseq is simulated
                  (default='null')
  --pacbio=s      Name of read file if pacbio is simulated
                  (default='null')
  --fmedian=i     Median of fragment size at shearing
                  (default='2000')
  --fshape=d      Shape parameter of the fragment size distribution
                  (default='6.0')
  --num=i         Number of fragments 
                  (default='1000000')
  --pblen=i       PacBio: Average (polymerase) read length
                  (default='30000')
  --illen=i       Illumina: read length
                  (default='300')
  --seed=i        Random seed, 0 for a random seed
                  (default='0')
  --help          Display this usage and exit
                  (default='false')




~~~~~~~~~~~~~
Usage samples
~~~~~~~~~~~~~

