------------------------------------------------
*jsa.hts.countReads*: Count reads from bam files 
------------------------------------------------


~~~~~~~~
Synopsis
~~~~~~~~

*jsa.hts.countReads*: Count the number of reads in some regions from a sorted, indexed bam file

~~~~~
Usage
~~~~~
::

   jsa.hts.countReads [options]

~~~~~~~
Options
~~~~~~~
  --bamFile=s     Name of the bam file
                  (REQUIRED)
  --bedFile=s     Name of the regions file in bed format
                  (REQUIRED)
  --output=s      Name of output file, - for from standard out.
                  (default='-')
  --flanking=i    Size of the flanking regions, effectively expand the region by flanking
                  (default='0')
  --qual=i        Minimum quality
                  (default='0')
  --filterBits=i  Filter reads based on flag. Common values:
                   0    no filter
                   256  exclude secondary alignment 
                   1024 exclude PCR/optical duplicates
                   2048 exclude supplementary alignments
                  (default='0')
  --contained     Count reads contained in the region
                  (default='false')
  --overlap       Count number of read overlap with the region
                  (default='false')
  --span          Count reads span the region
                  (default='false')
  --help          Display this usage and exit
                  (default='false')




