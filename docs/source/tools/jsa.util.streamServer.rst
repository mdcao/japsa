----------------------------------------------------------------
*jsa.util.streamServer*: Receiving streaming data over a network
----------------------------------------------------------------

*jsa.util.streamServer* implements a server that listen at a specified port. 
Upon receiving data from a client, it forwards the stream data to standard 
output. *jsa.util.streamServer* and *jsa.util.streamClient* can be used to
set up streaming applications such as real-time analyses. By default, 
the server listens on port 3456, unless specified otherwise.

~~~~~~~~
Synopsis
~~~~~~~~

*jsa.util.streamServer*:Listen for input from a stream and output to the standard output

~~~~~
Usage
~~~~~
::

   jsa.util.streamServer [options]

~~~~~~~
Options
~~~~~~~
  --port=i        Port to listen to
                  (default='3456')
  --help          Display this usage and exit
                  (default='false')


~~~~~~~~
See also
~~~~~~~~

jsa.util.streamClient_, jsa.np.filter.jsa.np.f5reader

.. _jsa.util.streamClient: jsa.util.streamClient.html





