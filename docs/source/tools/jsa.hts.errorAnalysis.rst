----------------------------------------------------------
*jsa.hts.errorAnalysis*: Error analysis of sequencing data
----------------------------------------------------------

*jsa.hts.errorAnalysis* assesses the error profile of sequencing data by getting the numbers
of errors (mismatches, indels etc) from a bam file. Obviously, it does not distinguish
sequencing errors from mutations, and hence consider mutations as errors. It is best to use
with the bam file from aligning sequencing reads to a reliable assembly of the sample.

~~~~~~~~
Synopsis
~~~~~~~~

*jsa.hts.errorAnalysis*: Error analysis of sequencing data

~~~~~
Usage
~~~~~
::

   jsa.hts.errorAnalysis [options]

~~~~~~~
Options
~~~~~~~
  --bamFile=s     Name of bam file
                  (REQUIRED)
  --reference=s   Name of reference genome
                  (REQUIRED)
  --pattern=s     Pattern of read name, used for filtering
                  (default='null')
  --qual=i        Minimum quality required
                  (default='0')
  --help          Display this usage and exit
                  (default='false')




