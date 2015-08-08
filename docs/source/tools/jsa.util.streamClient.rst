----------------------------------------------------
*jsa.util.streamClient*: Streams data over a network
----------------------------------------------------

*jsa.util.streamClient* streams data over the network to a listening server
(*jsa.util.streamServer*).
 
~~~~~~~~
Synopsis
~~~~~~~~

*jsa.util.streamClient*:Listen for input from the standard input and output to a stream

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

jsa.np.streamServer_, jsa.np.filter_, jsa.np.f5reader_

.. _jsa.np.streamServer: jsa.np.streamServer.html
.. _jsa.np.filter: jsa.np.filter.html
.. _jsa.np.f5reader: jsa.np.f5reader.html


 
