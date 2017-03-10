---------------------------------------------------------
*Expert Model*: tool for compression of genomic sequences 
---------------------------------------------------------

*jsa.xm.compress* in the implementation of the expert model (XM) algorithm for 
compression of genomics sequences. The source code is included in the 
`Japsa package <http://mdcao.github.io/japsa/>`_. 
Please see check the installation_ page for instructions.  

.. _installation: ../install.html

~~~~~~~~
Synopsis
~~~~~~~~

*jsa.xm.compress*: Compression of DNA/protein sequences

~~~~~
Usage
~~~~~
::

   jsa.xm.compress [options] file1 file2 ...

~~~~~~~
Options
~~~~~~~
  --hashSize=i    Hash size
                  (default='11')
  --context=i     Length of the context
                  (default='15')
  --limit=i       Expert Limit
                  (default='200')
  --threshold=d   Listen threshold
                  (default='0.15')
  --chance=i      Chances
                  (default='20')
  --binaryHash    Use binary hash or not
                  (default='false')
  --offsetType=s  Way of update offset/palindrome expert: possible value count, subs
                  (default='counts')
  --real=s        File name of the real compression
                  (default='null')
  --decode=s      File name of the encoded
                  (default='null')
  --output=s      The output file of decoded file
                  (default='decoded')
  --info=s        File name of the infomation content
                  (default='null')
  --markov=s      File name of the markov infomation content
                  (default='null')
  --optimise      Running in optimise mode, just report the entropy,recommended for long sequence
                  (default='false')
  --checkPoint=i  Frequency of check point
                  (default='1000000')
  --hashType=s    Type of Hash table: hash=hashtable, sft=SuffixTree,sfa = SuffixArray
                  (default='hash')
  --selfRep       Propose experts from the sequence to compressed?
                  (default='true')
  --help          Display this usage and exit
                  (default='false')




~~~~~~~~
Citation
~~~~~~~~

If you find XM useful for your research, please cite

Cao MD, Dix TI, Allison L, and Mears C, 
A simple statistical algorithm for biological sequence compression
Data Compression Conference, 2007 (DCC'07), Snowbird, UT, pp43-52.

