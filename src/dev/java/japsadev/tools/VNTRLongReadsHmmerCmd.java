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

/****************************************************************************
 *                           Revision History                                
 * 15/05/2014 - Minh Duc Cao: Started
 *  
 ****************************************************************************/
package japsadev.tools;

import japsa.bio.alignment.MultipleAlignment;
import japsa.bio.tr.TandemRepeat;
import japsa.bio.tr.TandemRepeatVariant;
import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.SequenceOutputStream;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.seq.XAFReader;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

/**
 * VNTR typing using long reads
 * 
 */

@Deployable(scriptName = "jsa.dev.longreadsHmmer", scriptDesc = "VNTR typing using long reads using hmmer")
public class VNTRLongReadsHmmerCmd extends CommandLine{
	static Alphabet alphabet = Alphabet.DNA16();

	public VNTRLongReadsHmmerCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());		

		addStdInputFile();
		addString("bamFile", null, "Name of the bam file", true);
		addString("output", "-",
				"Name of the output file, -  for stdout");
		addString("xafFile", null, "Name of the regions file in xaf",
				true);
		addString("prefix", "",
				"Prefix of temporary files, if not specified, will be automatically generated");
		addInt("flanking", 40, "Size of the flanking regions");
		addInt("qual", 0, "Minimum quality");
		addInt(
				"nploidy",
				2,
				"The ploidy of the genome 1 =  happloid, 2 = diploid. Currenly only support up to 2-ploidy");

		addString("msa", "kalign",
				"Name of the msa method, support kalign and clustalo");
		// cmdLine.addBoolean("contained", false,
		// "true: Reads contained in the region; false: reads overlap with the region");
		addStdHelp();
	}

	public static void main(String[] args) throws Exception,
	InterruptedException {
		/*********************** Setting up script ****************************/
		CommandLine cmdLine = new VNTRLongReadsHmmerCmd();		
		args = cmdLine.stdParseLine(args);	

		// Get options

		String prefix = cmdLine.getStringVal("prefix");

		if (prefix == null || prefix.length() == 0) {
			prefix = "p" + System.currentTimeMillis();
		}

		int flanking = cmdLine.getIntVal("flanking");
		if (flanking < 0)
			flanking = 0;

		int qual = cmdLine.getIntVal("qual");

		int np = cmdLine.getIntVal("nploidy");
		if (np > 2) {
			System.err
			.println("The program currenly only support haploid and diployd. Enter nploidy of 1 or 2");
			System.exit(1);
		}
		/**********************************************************************/

		String cmd;
		if ("kalign".equals(cmdLine.getStringVal("msa"))) {
			cmd = "kalign -gpo 60 -gpe 10 -tgpe 0 -bonus 0 -q -i " + prefix
					+ "i.fasta -o " + prefix + "o.fasta";
		} else {
			cmd = "clustalo --force -i " + prefix + "i.fasta -o " + prefix
					+ "o.fasta";
		}

		SequenceOutputStream outOS = SequenceOutputStream
				.makeOutputStream(cmdLine.getStringVal("output"));

		String[] headers = TandemRepeatVariant.SIMPLE_HEADERS;
		if (np > 1) {
			headers = TandemRepeatVariant.SIMPLE_HEADERS2;
		}

		TandemRepeatVariant.printHeader(outOS, headers);

		String strFile = cmdLine.getStringVal("xafFile");

		// TODO: make it multiple sequence
		Sequence seq = SequenceReader.getReader(cmdLine.getStringVal("input"))
				.nextSequence(alphabet);

		/**********************************************************************/
		//ArrayList<TandemRepeat> myList = TandemRepeat.readFromFile(
		//		SequenceReader.openFile(strFile), new ArrayList<String>());


		XAFReader xafReader = new XAFReader(strFile);
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader = SamReaderFactory.makeDefault().open(new File(cmdLine.getStringVal("bamFile")));	


		int _tIndex = 0;
		while (xafReader.next() != null){			
			_tIndex ++;
			//for (TandemRepeat str : myList) {
			TandemRepeat str = TandemRepeat.read(xafReader);
			System.out.println(_tIndex + "  " + str.getChr() + " " + str.getStart() + " " + str.getEnd());
			if (str.getPeriod() <= 4)
				continue;

			int start = str.getStart() - flanking;
			int end = str.getEnd() + flanking;

			if (end > seq.length())
				end = seq.length();

			if (start < 1)
				start = 1;

			SAMRecordIterator iter = reader.query(str.getParent(), start, end,
					false);

			int maxAlign = 300;
			MultipleAlignment ma = new MultipleAlignment(maxAlign, seq);

			while (iter.hasNext()) {
				SAMRecord rec = iter.next();
				//FIXME: this should be handdled by making sure sequence reads present
				if (rec.getReadLength() < 10)
					continue;

				// Check qualilty
				if (rec.getMappingQuality() < qual) {
					continue;
				}

				// Only reads that fully span the repeat and flankings
				if (rec.getAlignmentStart() > start)
					continue;
				if (rec.getAlignmentEnd() < end)
					continue;


				ma.addRead(rec);				
			}// while
			iter.close();
			// os.close();

			// seq.subSequence(start, end).writeFasta(prefix+"r.fasta");
			double var = 0;
			// double evidence = 0;

			TandemRepeatVariant trVar = new TandemRepeatVariant();
			trVar.setTandemRepeat(str);

			if (ma.printFasta(start, end, prefix + "i.fasta") >= 4) {
				Logging.info("Running " + cmd);
				Process process = Runtime.getRuntime().exec(cmd);
				process.waitFor();
				Logging.info("Done " + cmd);

				SequenceReader msaReader 
				= FastaReader.getReader(prefix	+ "o.fasta");
				ArrayList<Sequence> seqList = new ArrayList<Sequence>();
				Sequence nSeq = null;
				while ((nSeq = msaReader.nextSequence(Alphabet.DNA16())) != null) {
					seqList.add(nSeq);
				}

				String target = "pacbio_" + str.getID(); 
				//str.getChr()+"_"+str.getStart()+"_"+str.getEnd();

				if (np >= 2) {
					/*****************************************************************
					//Logging.info("Run hmmsearcg " + "hmmsearch -o "+target+".hmmsearch " + target + ".hmm " + prefix + "i.fasta");
					process = 
							Runtime.getRuntime().exec("hmmsearch -o "+target+".hmmsearch " + target + ".hmm " + prefix + "i.fasta");
					process.waitFor();
					//Read hmm results

					BufferedReader hmmReader = SequenceReader.openFile(target+".hmmsearch");
					HashMap<String, ReadInstance> instances = new HashMap<String, ReadInstance>();

					String hmmLine = "";					
					while ((hmmLine = hmmReader.readLine()) != null){
						String[] hmmToks = hmmLine.trim().split("  *");						
						//Logging.info(hmmLine + " " +hmmToks.length);
						if(hmmToks.length != 4)
							continue;


						if(!hmmToks[0].startsWith("SSS")){
							continue;
						}
						try{
							int hmmMatchStart = Integer.parseInt(hmmToks[3]),
									hmmMatchEnd   = Integer.parseInt(hmmToks[1]);

							ReadInstance inst = instances.get(hmmToks[0]);

							if (inst == null){
								inst = new ReadInstance(2, hmmToks[0]);
								inst.increase(0, 1) ;
								inst.increase(1,hmmMatchEnd - hmmMatchStart) ;							

								instances.put(hmmToks[0], inst);
							}else{
								inst.increase(0, 1) ;
								inst.increase(1,hmmMatchEnd - hmmMatchStart) ;	
							}
						}catch (Exception e){
							Logging.warn("Problems: " + e.getMessage() +" " +hmmLine);
						}
					}
					hmmReader.close();


					HashSet<String> set1  = new HashSet<String> ();



					Dataset data = new DefaultDataset();
					for (ReadInstance inst:instances.values()){					    
						data.add(inst);
					}
					if (data.size() > 1){
						Clusterer km = new KMeans(2);
						Dataset[] clusters = km.cluster(data);					

						for (int i = 0; i < clusters[0].size();i++){
							set1.add(clusters[0].get(i).toString());
						}
					}


					ArrayList<Sequence> seqList1 = new ArrayList<Sequence>(),
							seqList2 = new ArrayList<Sequence>();

					for (Sequence s:seqList){
						if (set1.contains(s.getName())){
							seqList1.add(s);
						}else{
							seqList2.add(s);
						}
					}

					int llength = seqList.get(0).length();


					System.err.print("#Call :" + _tIndex + " 1:");
					for (int s = 0; s < seqList1.size(); s++) {
						System.err.print(seqList1.get(s).getName() + ",");
					}
					System.err.println();


					System.err.print("#Call :" + _tIndex + " 2:");
					for (int s = 0; s < seqList2.size(); s++) {
						System.err.print(seqList2.get(s).getName() + ",");
					}
					System.err.println();


					int gaps = call(seqList1);					

					var = (llength - gaps - end + start) * 1.0
							/ str.getPeriod();

					trVar.setVar(var);

					gaps = call(seqList2);					
					var = (llength - gaps - end + start) * 1.0
							/ str.getPeriod();

					trVar.setVar2(var);

					trVar.addEvidence(seqList1.size());
					trVar.addEvidence2(seqList2.size());

					// trVar.setVar2(var);
					// trVar.addEvidence(seqList.size());

					 /*****************************************************************/
				} else {// nploidy ==1
					int llength = seqList.get(0).length();

					int gaps = call(seqList);

					var = (llength - gaps - end + start) * 1.0
							/ str.getPeriod();

					trVar.setVar(var);
					trVar.addEvidence(seqList.size());
				}

			}// if

			Process process = 
					Runtime.getRuntime().exec("rm -f " + prefix + _tIndex + "i.fasta");
			process.waitFor();

			process = 
					Runtime.getRuntime().exec("cp " + prefix + "i.fasta " + prefix + _tIndex + "i.fasta");
			process.waitFor();			

			outOS.print(trVar.toString(headers));
			outOS.print('\n');
			ma.printAlignment(start, end);
			ma.reduceAlignment(start, end). printAlignment(start, end);
		}// for

		reader.close();
		outOS.close();

	}
	/*****************************************************************
	 * Temporary commented out for independence of javaml
	 * 
	static class ReadInstance extends DenseInstance{

		private static final long serialVersionUID = 1L;
		String readName;

		public ReadInstance(int nAtt, String name) {
			super(nAtt);
			readName = name;
		}		

		public void increase(int attNo, double added){
			put(attNo,value(attNo) + added);
		}

		public String toString(){
			return readName;
		}
	}
/*****************************************************************/

	/**
	 * 
	 * @param seqList
	 * @param start
	 *            : the start index of the list (inclusive)
	 * @param end
	 *            : the end index of the list (exclusive)
	 */
	static int call(ArrayList<Sequence> seqList, int indexStart, int indexEnd) {
		if (indexEnd <= indexStart)
			return 0;	

		// Get consensus
		int gaps = 0;
		Sequence nSeq = new Sequence(Alphabet.DNA6(), seqList.get(0).length(),
				"consensus");
		int[] votes = new int[6];
		for (int i = 0; i < nSeq.length(); i++) {			
			Arrays.fill(votes, 0);
			for (int s = indexStart; s < indexEnd; s++) {
				votes[seqList.get(s).symbolAt(i)]++;
			}
			byte best = 0;
			for (byte b = 1; b < 6; b++)
				if (votes[b] > votes[best])
					best = b;

			nSeq.setBase(i, best);
			if (best == 5)
				gaps++;
		}// for
		return gaps;
	}

	static int call(ArrayList<Sequence> seqList) {
		return call(seqList,0,seqList.size());

	}
	/************************************************
	static int[][] distance(ArrayList<Sequence> seqList) {
		int[][] dis = new int[seqList.size()][seqList.size()];

		for (int i = seqList.size() - 1; i > 1; i--) {
			for (int j = i - 1; j >= 0; j--) {
				dis[i][j] = dis[j][i] = distance(seqList.get(i), seqList.get(j));
			}
		}

		return dis;
	}

	static int distance(Sequence s1, Sequence s2) {
		int dis = 0;
		for (int i = 0; i < s1.length(); i++) {
			if (s1.getBase(i) != s2.getBase(i)) {
				if (s1.getBase(i) == DNA.GAP || s2.getBase(i) == DNA.GAP)
					dis += 1;
				else
					dis += 0;
			}
		}
		return dis;
	}
/************************************************/
}
