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
 * 18/10/2013 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsadev.tools;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import japsa.bio.hts.scaffold.AlignmentRecord;
import japsa.bio.hts.scaffold.Contig;



@Deployable(
		scriptName = "jsa.dev.flankDetect",
		scriptDesc = "Detect flanking sequences from both ends of nanopore reads"
		)
public class FlankSeqsDetectorCmd extends CommandLine{
	public FlankSeqsDetectorCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("flankFile",null,"Flank sequences file, maximum 2 sequences",true);
		addString("bamFile",null,"Bam file",true);
		addDouble("qual", 1, "Mininum quality");
		addStdHelp();	
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{
		/*********************** Setting up script ****************************/		 
		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new FlankSeqsDetectorCmd();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/
		String flankSeqsFile= cmdLine.getStringVal("flankFile");
		String bamFile = cmdLine.getStringVal("bamFile");
		double qual = cmdLine.getDoubleVal("qual");
		
		SequenceReader seqReader = SequenceReader.getReader(flankSeqsFile);
		Sequence seq;
		ArrayList<Contig> flankSeqs = new ArrayList<>(); 
		int index=0;
		while ((seq = seqReader.nextSequence(Alphabet.DNA())) != null)
			flankSeqs.add(new Contig(index++,seq));
		
		seqReader.close();
		
		if(flankSeqs.size() > 2){
			System.err.println("More than 2 sequences!");
			System.exit(1);
		}
		
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader = SamReaderFactory.makeDefault().open(new File(bamFile));
		SAMRecordIterator iter = reader.iterator();

		SAMRecord rec;
		AlignmentRecord curAlnRec;
		String curReadName = "";
		HashMap<String, ArrayList<AlignmentRecord>> map = new HashMap<>();
		
		while (iter.hasNext()) {
			rec = iter.next();			
			curReadName=rec.getReadName();

			if (rec.getReadUnmappedFlag() || rec.getMappingQuality() < qual){	
				if(!map.containsKey(curReadName))
					map.put(curReadName, new ArrayList<>());
				continue;	
			}
			Contig flk = flankSeqs.get(rec.getReferenceIndex());
			curAlnRec=new AlignmentRecord(rec, flk);
			if(curAlnRec.refEnd-curAlnRec.refStart < .5*flk.length())
				continue;
			//not too far from the tip of read
			else if(Math.min(-curAlnRec.readAlignmentEnd()+curAlnRec.readLength, curAlnRec.readAlignmentStart()) > 50){
				continue;
			}
			
			
			ArrayList<AlignmentRecord> recs = new ArrayList<>();
			recs.add(curAlnRec);
			if(!map.containsKey(curReadName) || map.get(curReadName).isEmpty()){
				map.put(curReadName, recs);
			}else if(map.get(curReadName).size()==1 ){
				AlignmentRecord prev = map.get(curReadName).get(0);
				//not overlap each other
				if((prev.readStart-curAlnRec.readStart)*(prev.readEnd-curAlnRec.readStart) <= 0 
					|| (prev.readStart-curAlnRec.readEnd)*(prev.readEnd-curAlnRec.readEnd) <= 0 )
					continue;
					//the in-between must longer than 100bp
				if(prev.readAlignmentEnd()-prev.readAlignmentStart() + curAlnRec.readAlignmentEnd()-curAlnRec.readAlignmentStart() > curAlnRec.readLength -100)
					continue;
				
				map.get(curReadName).add(curAlnRec);
				
			}else{
				System.err.println("Something wrong!");
				System.exit(1);
			}				
		
			

		}// while
		iter.close();
		/**********************************************************************/
		int totReadNum = map.keySet().size();
		ArrayList<String>  	flank0 = new ArrayList<String>(),
							flank1 = new ArrayList<String>(),
							flank2 = new ArrayList<String>();
		map.keySet().stream().forEach(r->{
			if(map.get(r).isEmpty())
				flank0.add(r);
			else if(map.get(r).size()==1)
				flank1.add(r);
			else if(map.get(r).size()==2)
				flank2.add(r);
		});
		System.out.println( "Total number of reads: " + totReadNum);
		System.out.println("Number of reads with 0 flank sequences: " + flank0.size());
		flank0.stream().forEach(r->System.out.println(r));
		System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
		
		System.out.println("Number of reads with 1 flank sequences: " + flank1.size());
		flank1.stream().forEach(r->System.out.println(r));
		System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
		
		System.out.println("Number of reads with 2 flank sequences: " + flank2.size());
		flank2.stream().forEach(r->System.out.println(r));
		System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
	}

}
