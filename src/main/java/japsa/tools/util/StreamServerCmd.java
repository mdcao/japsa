/*****************************************************************************
 * Copyright (c) Minh Duc Cao, Monash Uni & UQ, All rights reserved.         *
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  * 
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 * 3. Neither the names of the institutions nor the names of the contributors*
 *    may be used to endorse or promote products derived from this software  *
 *    without specific prior written permission.                             *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS   *
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, *
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR    *
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR         *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 ****************************************************************************/

/**************************     REVISION HISTORY    **************************
 * 5 Mar 2015 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/
package japsa.tools.util;

import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

import java.io.IOException;
import java.net.ServerSocket;
import java.net.Socket;

import com.google.common.io.ByteStreams;


/**
 * @author minhduc
 *
 */
@Deployable(
		scriptName = "jsa.util.streamServer",
		scriptDesc = "Listen for input from a stream and output to the standard output",
		seeAlso = "jsa.util.streamClient, jsa.np.filter.jsa.np.f5reader"
		//scriptDocs = "jsa.util.streamServer implements a server that listen at "
		//		+ "a specified port. Upon receiving data from a client, it forwards the stream "
		//		+ "data to standard output"
		)
public class StreamServerCmd extends CommandLine{
	public static int DEFAULT_PORT = 3456;
	public StreamServerCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addInt("port", DEFAULT_PORT,  "Port to listen to");	
		
		addStdHelp();		
	} 
	
	
/**
 * @param args
 * @throws InterruptedException 
 * @throws Exception 
 * @throws OutOfMemoryError 
 */
	public static void main(String[] args) throws IOException, InterruptedException{		 		
		CommandLine cmdLine = new StreamServerCmd();				
		args = cmdLine.stdParseLine(args);			
		/**********************************************************************/
				
		int port = cmdLine.getIntVal("port");
		
		ServerSocket serverSocket = new ServerSocket(port);
		Logging.info("Listen on port " + port);		
	    Socket clientSocket = serverSocket.accept();
	    Logging.info("Connection establised");
	    ByteStreams.copy(clientSocket.getInputStream(), System.out);
	    serverSocket.close();
	    Logging.info("Connection closed");	    
	}
}

/*RST*
----------------------------------------------------------------
*jsa.util.streamServer*: Receiving streaming data over a network
----------------------------------------------------------------

*jsa.util.streamServer* implements a server that listen at a specified port. 
Upon receiving data from a client, it forwards the stream data to standard 
output. *jsa.util.streamServer* and *jsa.util.streamClient* can be used to
set up streaming applications such as real-time analyses. By default, 
the server listens on port 3456, unless specified otherwise.

<usage>



*RST*/

