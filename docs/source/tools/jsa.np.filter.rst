---------------------------------------
*jsa.np.filter*: Filter sequencing data 
---------------------------------------

*jsa.np.filter* filters sequencing data based on sequence read type, length and
quality. Examples of its usage can be found on jsa.np.npreader_.

~~~~~~~~
Synopsis
~~~~~~~~

*jsa.np.filter*: Filter nanopore reads data from fastq file

~~~~~
Usage
~~~~~
::

   jsa.np.filter [options]

~~~~~~~
Options
~~~~~~~
  --input=s       Name of the input file, - for standard input
                  (REQUIRED)
  --output=s      Name of the output file, - for standard output
                  (REQUIRED)
  --lenMin=i      Minimum sequence length
                  (default='0')
  --lenMax=i      Minimum sequence length
                  (default='2147483647')
  --qualMin=d     Minimum average quality
                  (default='0.0')
  --qualMax=d     Maximum average quality
                  (default='1000.0')
  --group=s       Group need to be extracted, leave blank for selecting all groups
                  (default='')
  --excl2D        Exclude 2D reads
                  (default='false')
  --exclTemp      Exclude template reads
                  (default='false')
  --exclComp      Exclude complement reads
                  (default='false')
  --format=s      Format of the output file
                  (default='fastq')
  --help          Display this usage and exit
                  (default='false')


~~~~~~~~
See also
~~~~~~~~

jsa.np.npreader_, jsa.util.streamServer_, jsa.util.streamClient_

.. _jsa.np.npreader: jsa.np.npreader.html
.. _jsa.util.streamServer: jsa.util.streamServer.html
.. _jsa.util.streamClient: jsa.util.streamClient.html




