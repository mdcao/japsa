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
 * 21/07/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.seq.nanopore;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.Socket;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import japsa.seq.FastqSequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.nanopore.Fast5NPReader.BaseCalledFastq;
import japsa.util.DoubleArray;
import japsa.util.IntArray;
import japsa.util.JapsaException;
import japsa.util.net.StreamClient;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * Read nanopore data (read sequence, events, alignment, models etc) from a raw
 * (fast5) format. 
 * @author minhduc
 *
 */
public class NanoporeReaderStream{
    private static final Logger LOG = LoggerFactory.getLogger(NanoporeReaderStream.class);


    public String prepareIO(){
		String msg = null;
		try{
			sos = SequenceOutputStream.makeOutputStream(output);
			if (streamServers != null && streamServers.trim().length() > 0){
				@SuppressWarnings("resource")
				StreamClient streamClient = new StreamClient(streamServers);
				ArrayList<Socket>  sockets = streamClient.getSockets();
				networkOS = new ArrayList<SequenceOutputStream>(sockets.size());				
				for (Socket socket:sockets)					
					networkOS.add(new SequenceOutputStream(socket.getOutputStream()));			
			}			
		}catch (Exception e){
			msg = e.getMessage();
		}finally{
			if (msg != null){
				sos = null;
				return msg;
			}
		}
		return msg;
	}


	public void close() throws IOException{
        LOG.info("npReader closing");
		sos.close();
		if (networkOS != null){
			for (SequenceOutputStream out:networkOS)
				out.close();
		}
		if(dmplx!=null)
			dmplx.close();
        LOG.info("npReader closed");
		//done = true;		
	}

	double tempLength = 0, compLength = 0, twoDLength = 0;
	int tempCount = 0, compCount = 0, twoDCount = 0;
	IntArray lengths = new IntArray();
	DoubleArray qual2D = new DoubleArray(), qualComp = new DoubleArray(), qualTemp = new DoubleArray();
	IntArray lengths2D = new IntArray(), lengthsComp = new IntArray(), lengthsTemp = new IntArray();

	SequenceOutputStream sos;
	ArrayList<SequenceOutputStream> networkOS = null;
	public boolean stats, number;
	public String folder = null;
	public int minLength = 1;
	public volatile boolean wait = true;
	public boolean realtime = true;
	public int interval = 1, age = 30000;
	public boolean doFail = false;
	public String output = "-";
	public String streamServers = null;
    public boolean getTimeStamp = false;
    private double rps = 4000.0;

	
	public String format = "fastq";
	public boolean ready = true;
	private static final byte MIN_QUAL = '!';//The minimum quality
	
	public boolean exhaustive = false;
	Demultiplexer dmplx = null;
	private String bcFile = null;
	
	public void updateDemultiplexFile(String file){
		if(file==null)
			return;
		bcFile = file;
		try{
			dmplx = new Demultiplexer(file);
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	public String getBCFileName(){
		return bcFile;
	}

	/**
	 * Compute average quality of a read
	 * @param fq
	 * @return
	 */
	public static double averageQuality(FastqSequence fq){
		if (fq.length() > 0){
			double sumQual  = 0;
			for (int p = 0; p < fq.length(); p++){
				sumQual += (fq.getQualByte(p) - MIN_QUAL);

			}
			return (sumQual/fq.length());
		}	
		else return 0;
	}

	public void print(FastqSequence fq) throws IOException{
		if (format.equals("fasta"))
			fq.writeFasta(sos);
		else
			fq.print(sos);

		if (networkOS != null){
			for (SequenceOutputStream out:networkOS)
				if (format.equals("fasta"))
					fq.writeFasta(out);
				else
					fq.print(out);			//fq.print(out);
		}
	}

	void flush() throws IOException{
		sos.flush();
		if (networkOS != null){
			for (SequenceOutputStream out:networkOS)
				out.flush();
		}
	}


	public boolean readFastq2(String fileName) throws JapsaException, IOException{
		//LOG.info("Open " + fileName);
		try{					
			Fast5NPReader npReader = new Fast5NPReader(fileName);
			npReader.readFastq();
			if (getTimeStamp)
			    npReader.readTime();
			npReader.close();

			ArrayList<BaseCalledFastq> seqList = npReader.getFastqList();
			if(seqList == null || seqList.isEmpty())
				return false;
			else{
				for (BaseCalledFastq fq:seqList){
					if (fq.length() >= minLength){
						fq.setName((number?(getTotalFilesNumber() *3 + fq.type()) + "_":"") + fq.getName() +
                                (getTimeStamp? (" timestamp=" +(npReader.timeStamp/rps)):"")
                        );
						//do multiplexing here
						if(dmplx!=null)
							dmplx.clustering(fq);
						
						print(fq);						
						if (stats){						
							lengths.add(fq.length());
							double sumQual  = 0;
							for (int p = 0; p < fq.length(); p++){
								sumQual += (fq.getQualByte(p) - MIN_QUAL);
							}
							if (fq.isTwoDim()){
								lengths2D.add(fq.length());
								twoDCount ++;
								qual2D.add(sumQual/fq.length());
							}else if (fq.isComplement()){
								lengthsComp.add(fq.length());
								compCount ++;
								qualComp.add(sumQual/fq.length());
							}else if (fq.isTemplate()){
								lengthsTemp.add(fq.length());
								tempCount ++;
								qualTemp.add(sumQual/fq.length());								
							}
						}
					}
				}
			}

			//fileNumber ++;			
		}catch (JapsaException e){
			throw e;
		}catch (Exception e){
            LOG.error("Problem with reading " + fileName + ":" + e.getMessage());
			e.printStackTrace();
			return false;
		}		
		return true;
	}
	/*****************************************************************************/

	/**
	 * Read read sequence from a list of fast5 files.
	 * @param fileList
	 * @param gSos : output stream
	 * @param stats: print out statistics
	 * @throws IOException 
	 */
	HashSet<String> filesOK = new HashSet<String>(),
					filesSkipped = new HashSet<String>();
	public void readFast5() throws JapsaException, IOException{
		if (minLength < 1)
			minLength = 1;

        LOG.info("Start reading " + folder);


		//HashSet<String> filesDone = new HashSet<String>();

		while (wait){
			//Do main
			final long now = System.currentTimeMillis();

            LOG.info("Start reading  " + now );
			try{
				Files.walk(Paths.get(folder))
				//is a file
				.filter(Files::isRegularFile)
				//fail folder
				.filter(p -> {
					try{
						Path failFolderPath= Paths.get(folder+ File.separator + "fail");
						if(failFolderPath.toFile().isDirectory())
							return doFail || !Files.isSameFile(p.getParent(), failFolderPath);
						else 
							return true;
					}catch(IOException e){
						e.printStackTrace();
						return false;
					}
				})
				//fast5 file
				.filter(p -> p.toString().endsWith("fast5"))
				//age is old enough
				.filter(p -> {		        	
					try{					
						return now - Files.getLastModifiedTime(p).toMillis() > age;		        		
					}catch (IOException e1) {
						e1.printStackTrace();
						return false;
					}
				})
				//not read before
				.filter(p -> !filesOK.contains(p.toString()) && (exhaustive || !filesSkipped.contains(p.toString())) )
				//read
				.forEach(p -> {
				//	System.out.println(p);
					try {
						if (readFastq2(p.toString())){
							filesOK.add(p.toString());
							filesSkipped.remove(p.toString());
						}
						else
							filesSkipped.add(p.toString());
						
					} catch (JapsaException | IOException e) {
						e.printStackTrace();
					}		
	
					if(!wait)
						throw new BreakException("Stopping");
				});	

			}catch(BreakException e){
                LOG.info("Stop to read on directory " + folder);
			}
		/*******************************************************/
			if (!realtime)
				break;

			for (int x = 0; x < interval && wait; x++){
				try {
					Thread.sleep(1000);
				} catch (InterruptedException e) {					
					e.printStackTrace();
				}
			}

		}//while			
        LOG.info("EXITING");


		if (stats){
            LOG.info("Getting stats ... ");
			int [] ls = lengths.toArray();
			if (ls.length ==0){
                LOG.info("Open " + getTotalFilesNumber() + " files");
                LOG.info("Found 0 reads");
			}else{
				Arrays.sort(ls);

				long baseCount = 0;						
				for (int i = 0; i < ls.length; i++)
					baseCount += ls[i];

				double mean = baseCount / ls.length;
				double median = ls[ls.length/2];
				long sum = 0;
				int quantile1st = 0, quantile2nd = 0, quantile3rd = 0;
				for (int i = 0; i < ls.length; i++){
					sum += ls[i];
					if (quantile1st == 0 && sum >= baseCount / 4)
						quantile1st = i;

					if (quantile2nd == 0 && sum >= baseCount / 2)
						quantile2nd = i;

					if (quantile3rd == 0 && sum >= baseCount * 3/ 4)
						quantile3rd = i;
				}

                LOG.info("Open " + getTotalFilesNumber() + " files");
                LOG.info("Read count = " + ls.length + "(" + tempCount + " templates, " + compCount + " complements and " + twoDCount +"  2D)");
                LOG.info("Base count = " + baseCount);
                LOG.info("Longest read = " + ls[ls.length - 1] + ", shortest read = " + ls[0]);
                LOG.info("Average read length = " + mean);
                LOG.info("Median read length = " + median);
                LOG.info("Quantile first = " + ls[quantile1st] + " second = " + ls[quantile2nd] + " third = " + ls[quantile3rd]);

				if (qual2D.size() > 0){
					double sumQual = 0;
					double sumQualSq = 0;

					for (int i = 0; i < qual2D.size();i++){
						sumQual += qual2D.get(i);
						sumQualSq += qual2D.get(i) * qual2D.get(i);
					}

					double meanQual = sumQual / qual2D.size();
					double stdQual = Math.sqrt(sumQualSq / qual2D.size() - meanQual * meanQual);

                    LOG.info("Ave 2D qual " +meanQual + " " + qual2D.size() + " std = " + stdQual);
				}

				if (qualTemp.size() > 0){
					double sumQual = 0;
					double sumQualSq = 0;

					for (int i = 0; i < qualTemp.size();i++){
						sumQual += qualTemp.get(i);
						sumQualSq += qualTemp.get(i) * qualTemp.get(i);
					}

					double meanQual = sumQual / qualTemp.size();
					double stdQual = Math.sqrt(sumQualSq / qualTemp.size() - meanQual * meanQual);

                    LOG.info("Ave Temp qual " +meanQual + " " + qualTemp.size() + " std = " + stdQual);
				}	

				if (qualComp.size() > 0){
					double sumQual = 0;
					double sumQualSq = 0;


					for (int i = 0; i < qualComp.size();i++){
						sumQual += qualComp.get(i);
						sumQualSq += qualComp.get(i) * qualComp.get(i);
					}

					double meanQual = sumQual / qualComp.size();
					double stdQual = Math.sqrt(sumQualSq / qualComp.size() - meanQual * meanQual);

                    LOG.info("Ave Comp qual " + meanQual + " " + qualComp.size() + " std = " + stdQual);
				}	
			}
			printToFile("stats");
		}
	}

	public void printToFile(String prefix) throws IOException{
		if(prefix.length() < 1)
			prefix = "out";
		BufferedWriter lenFile = new BufferedWriter(new PrintWriter(prefix + ".len")),
				qualTempFile = new BufferedWriter(new PrintWriter(prefix + ".temp.qual")),
				qualCompFile = new BufferedWriter(new PrintWriter(prefix + ".comp.qual")),
				qual2DFile = new BufferedWriter(new PrintWriter(prefix + ".2d.qual"));
		for(int i=0; i < lengths.size(); i++){
			lenFile.write(lengths.get(i) + "\n");
		}
		for(int i=0; i < qualTemp.size(); i++){
			qualTempFile.write(new DecimalFormat("#0.000").format(qualTemp.get(i)) + "\n");
		}
		for(int i=0; i < qualComp.size(); i++){
			qualCompFile.write(new DecimalFormat("#0.000").format(qualComp.get(i)) + "\n");
		}
		for(int i=0; i < qual2D.size(); i++){
			qual2DFile.write(new DecimalFormat("#0.000").format(qual2D.get(i)) + "\n");
		}
		lenFile.close();
		qualTempFile.close();
		qualCompFile.close();
		qual2DFile.close();
	}
	
	public synchronized int getTotalFilesNumber(){
		return filesOK.size() + filesSkipped.size();
	}
	public synchronized int getOKFilesNumber(){
		return filesOK.size();
	}
	public synchronized int getSkippedFilesNumber(){
		return filesSkipped.size();
	}
	
}
