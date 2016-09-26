----------------------------------------------------
*jsa.util.streamClient*: Streams data over a network
----------------------------------------------------

*jsa.util.streamClient* streams data over the network to a listening server
(*jsa.util.streamServer*).
 
~~~~~~~~
Synopsis
~~~~~~~~

*jsa.util.streamClient*: Forward data from a stream input or a file over the network to a jsa.util.streamServer

~~~~~
Usage
~~~~~
::

   jsa.util.streamClient [options]

~~~~~~~
Options
~~~~~~~
  --input=s       Name of the input file, - for standard input
                  (REQUIRED)
  --server=s      Stream output to one or more servers, format IP:port,IP:port
                  (REQUIRED)
  --help          Display this usage and exit
                  (default='false')


~~~~~~~~
See also
~~~~~~~~

jsa.util.streamServer_, jsa.np.filter_, jsa.np.npreader_

.. _jsa.util.streamServer: jsa.util.streamServer.html
.. _jsa.np.filter: jsa.np.filter.html
.. _jsa.np.npreader: jsa.np.npreader.html


 
