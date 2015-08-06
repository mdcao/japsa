----------------------------------------
*jsa.np.filter*: Filter sequencing reads
----------------------------------------

*jsa.np.filter* is a program that filters sequencing reads based
on the read length, quality and type. The program can read/write from/to a file
or a stream and can be integrated into a streaming pipeline.

*jsa.np.filter* is included in the `Japsa package <http://mdcao.github.io/japsa/>`_.


~~~~~
Usage
~~~~~
::

   jsa.np.filter [options]
   
~~~~~~~~
Options:
~~~~~~~~

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
jsa.util.streamServer, jsa.util.streamClient, jsa.np.f5reader_.

.. _jsa.np.f5reader: jsa.np.f5reader.html


