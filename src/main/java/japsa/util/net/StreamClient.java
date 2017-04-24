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
package japsa.util.net;


import japsa.tools.util.StreamServerCmd;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.Closeable;
import java.io.IOException;
import java.net.Socket;
import java.net.UnknownHostException;
import java.util.ArrayList;

/**
 * @author minhduc
 *
 */
public class StreamClient implements Closeable{
    private static final Logger LOG = LoggerFactory.getLogger(StreamClient.class);
	private ArrayList<Socket> sockets;

	public StreamClient(String serverList){
		sockets = new ArrayList<Socket> ();
		String [] servers = serverList.split(",");
		for (String server:servers){
			String [] toks = server.trim().split(":");
			int portNumber = StreamServerCmd.DEFAULT_PORT;
			if (toks.length > 1)
				portNumber = Integer.parseInt(toks[1]);			
			LOG.info("Trying to connect " + toks[0] + ":" + portNumber);
			try {
				Socket socket = new Socket(toks[0], portNumber);
				sockets.add(socket);
				LOG.info("Connection to " + toks[0] + ":" + portNumber + " established");
			} catch (UnknownHostException e) {
				//e.printStackTrace();
				LOG.warn("Could not connect to " + toks[0] + ":" + portNumber);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}


		}		
	}
	public ArrayList<Socket>  getSockets(){
		return sockets;
	}
	/* (non-Javadoc)
	 * @see java.io.Closeable#close()
	 */
	@Override
	public void close() throws IOException {
		for (Socket socket:sockets){
			socket.close();
			LOG.info("Connection to " + socket.getRemoteSocketAddress() + " closed");
		}
	}	

}
