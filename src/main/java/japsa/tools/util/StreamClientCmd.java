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
import japsa.util.net.StreamClient;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.Socket;
import java.util.Date;

/**
 * @author minhduc
 *
 */
@Deployable(
	scriptName = "jsa.util.streamClient",
	scriptDesc = "Forward data from a stream input or a file over the network to a jsa.util.streamServer",
	seeAlso = "jsa.util.streamServer, jsa.np.filter, jsa.np.npreader"
	)

public class StreamClientCmd extends CommandLine{	
	public StreamClientCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addStdInputFile();
		addString("server",null, "Stream output to one or more servers, format IP:port,IP:port",true);
		
		addStdHelp();
	} 


	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 * @throws Exception 
	 * @throws OutOfMemoryError 
	 */
	public static void main(String[] args) throws IOException{		 		
		CommandLine cmdLine = new StreamClientCmd();				
		args = cmdLine.stdParseLine(args);						
		/**********************************************************************/
		String input = cmdLine.getStringVal("input");
		StreamClient client = new StreamClient(cmdLine.getStringVal("server"));
		Logging.info("Connection established at " + new Date());

		InputStream ins = input.equals("-")? System.in : new FileInputStream(input);
		byte[] buffer = new byte[8192];

		while (true){
			int ret = ins.read(buffer);
			if (ret < 0)
				break;
			int count = 0;
			
			for (Socket socket:client.getSockets()){
				if (!socket.isClosed()){
					try {
						socket.getOutputStream().write(buffer,0, ret);
						count ++;
					} catch (IOException e) {
						Logging.info("Connection to " + socket.getRemoteSocketAddress() + " closed at "+ new Date());
						socket.close();
					}					
				}				
			}
			if (count <= 0)
				break;
		}
		client.close();
	}
}
/*RST*
----------------------------------------------------
*jsa.util.streamClient*: Streams data over a network
----------------------------------------------------

*jsa.util.streamClient* streams data over the network to a listening server
(*jsa.util.streamServer*).
 
 <usage>
 
*RST*/
