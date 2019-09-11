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

package japsa.tools.bio.np;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

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
		scriptName = "jsa.np.flankDetect",
		scriptDesc = "Detect flanking sequences from both ends of nanopore reads"
		)
public class FlankSeqsDetectorCmd extends CommandLine{
    private static final Logger LOG = LoggerFactory.getLogger(FlankSeqsDetectorCmd.class);

	public FlankSeqsDetectorCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("flankFile",null,"Flank sequences file, maximum 2 sequences",true);
		addString("refFile","","Reference sequences");
		addString("bamFile",null,"Bam file",true);
		addDouble("qual", 1, "Mininum quality");
		addInt("insert", 10, "Minimum length of insert sequence in-between 2 flanking sequences");
		addInt("tips", 20, "Maximum percentage of the overhangs compared to the corresponding flanking sequence");
		addDouble("cover", 80, "Mininum percentage of flank sequence coverage for a valid alignment");

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
		String 	flankSeqsFile= cmdLine.getStringVal("flankFile"),
				refSeqsFile= cmdLine.getStringVal("refFile");
		String bamFile = cmdLine.getStringVal("bamFile");
		double 	qual = cmdLine.getDoubleVal("qual"),
				flkCov = cmdLine.getDoubleVal("cover");
		int insertLength = cmdLine.getIntVal("insert"),
			tipsPercentage = cmdLine.getIntVal("tips");
		
		SequenceReader seqReader = SequenceReader.getReader(flankSeqsFile);
		Sequence seq;
		HashMap<String, Contig> refSeqs = new HashMap<>();
		ArrayList<Contig> flankSeqs=new ArrayList<>();
							
		int index=0;
		while ((seq = seqReader.nextSequence(Alphabet.DNA())) != null)
			flankSeqs.add(new Contig(index++,seq));
		seqReader.close();
		
		if(!refSeqsFile.isEmpty()){
			seqReader = SequenceReader.getReader(refSeqsFile);
			index=0;
			while ((seq = seqReader.nextSequence(Alphabet.DNA())) != null)
				refSeqs.put(seq.getName(), new Contig(index++,seq));
			seqReader.close();
		}
		
		if(flankSeqs.size() > 2){
			System.err.println("More than 2 sequences!");
			System.exit(1);
		}
		
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader = SamReaderFactory.makeDefault().open(new File(bamFile));
		SAMRecordIterator iter = reader.iterator();

		SAMRecord curSAMRecord=null;
		AlignmentRecord curAlnRecord=null;
		FlankRecord fr=null;
		String readID = "";

		while (iter.hasNext()) {
	
			try {
				curSAMRecord = iter.next();
			}catch(Exception e) {
				LOG.warn("Ignore one faulty SAM record: \n {}", e.getMessage());
				continue;
			}
			
			if (curSAMRecord.getReadUnmappedFlag() || curSAMRecord.getMappingQuality() < qual || curSAMRecord.isSecondaryAlignment()){		
				LOG.info("Ignore unmapped, low-quality or not primarily map record from {}...", curSAMRecord.getReadName());
				continue;		
			}
								
			if (!readID.equals(curSAMRecord.getReadName())){
				//output prev
				if(fr!=null)
					System.out.println(fr.printJunctionOnRef());
					
				//update for next
				readID = curSAMRecord.getReadName();
				fr=new FlankRecord(readID);
			}
			//adding info to current FlankRecord
			Contig ctg=null;
			String refName = curSAMRecord.getReferenceName();
			if(refSeqs.containsKey(refName)){
				ctg=refSeqs.get(refName);
				curAlnRecord=new AlignmentRecord(curSAMRecord, ctg);
				//increase score if primary
				if(!curSAMRecord.getSupplementaryAlignmentFlag())
					curAlnRecord.score+=100;

				if(curAlnRecord.readAlignmentEnd()-curAlnRecord.readAlignmentStart() < insertLength){
					LOG.info("Ignore integration size too short: {}", curAlnRecord.toString());
					continue;
				}
				if(fr.refRec==null||fr.refRec.score < curAlnRecord.score)
					fr.refRec=curAlnRecord;
				
			}else{
				
				if(flankSeqs.get(0).getName().equals(refName)){
					ctg=flankSeqs.get(0);
					curAlnRecord=new AlignmentRecord(curSAMRecord, ctg);
					if(curAlnRecord.refEnd-curAlnRecord.refStart < (double)flkCov*ctg.length()/100.0)
						continue;
					//not too far from the tip of read
					else if(Math.min(-curAlnRecord.readAlignmentEnd()+curAlnRecord.readLength, curAlnRecord.readAlignmentStart()) > (double)ctg.length()*tipsPercentage/100.0){
						continue;
					}
					if(fr.f0Rec==null||fr.f0Rec.score < curAlnRecord.score)
						fr.f0Rec=curAlnRecord;
					
				}else if(flankSeqs.size()>1 && flankSeqs.get(1).getName().equals(refName)){
					ctg=flankSeqs.get(1);
					curAlnRecord=new AlignmentRecord(curSAMRecord, ctg);
					if(curAlnRecord.refEnd-curAlnRecord.refStart < (double)flkCov*ctg.length()/100.0)
						continue;
					//not too far from the tip of read
					else if(Math.min(-curAlnRecord.readAlignmentEnd()+curAlnRecord.readLength, curAlnRecord.readAlignmentStart()) > (double)ctg.length()*tipsPercentage/100.0){
						continue;
					}
					if(fr.f1Rec==null||fr.f1Rec.score < curAlnRecord.score)
						fr.f1Rec=curAlnRecord;
				}else{
					LOG.info("Flank not found: {} != {}!", refName, flankSeqs.get(0).getName());
					System.exit(1);
				}
					
			}

			
		}// while
		iter.close();
		/**********************************************************************/
	}
	

//	public static void main(String[] args) throws IOException{
//		/*********************** Setting up script ****************************/		 
//		/*********************** Setting up script ****************************/		
//		CommandLine cmdLine = new FlankSeqsDetectorCmd();		
//		args = cmdLine.stdParseLine(args);
//		/**********************************************************************/
//		String flankSeqsFile= cmdLine.getStringVal("flankFile");
//		String bamFile = cmdLine.getStringVal("bamFile");
//		double 	qual = cmdLine.getDoubleVal("qual"),
//				flkCov = cmdLine.getDoubleVal("cover");
//		int insertLength = cmdLine.getIntVal("insert"),
//			tipsPercentage = cmdLine.getIntVal("tips");
//		
//		SequenceReader seqReader = SequenceReader.getReader(flankSeqsFile);
//		Sequence seq;
//		ArrayList<Contig> flankSeqs = new ArrayList<>(); 
//		int index=0;
//		while ((seq = seqReader.nextSequence(Alphabet.DNA())) != null)
//			flankSeqs.add(new Contig(index++,seq));
//		
//		seqReader.close();
//		
//		if(flankSeqs.size() > 2){
//			System.err.println("More than 2 sequences!");
//			System.exit(1);
//		}
//		
//		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
//		SamReader reader = SamReaderFactory.makeDefault().open(new File(bamFile));
//		SAMRecordIterator iter = reader.iterator();
//
//		SAMRecord rec;
//		AlignmentRecord curAlnRec;
//		String curReadName = "";
//		HashMap<String, HashMap<Contig, AlignmentRecord>> map = new HashMap<>();
////		HashMap<String, Sequence> readsMap = new HashMap<>(); 
//		while (iter.hasNext()) {
//			rec = iter.next();			
//			curReadName=rec.getReadName();
////			if(!readsMap.containsKey(curReadName)){
////				
////			}
//			
//			if (rec.getReadUnmappedFlag() || rec.getMappingQuality() < qual || rec.isSecondaryOrSupplementary()){	
//				if(!map.containsKey(curReadName))
//					map.put(curReadName, new HashMap<>());
//				continue;	
//			}
//			Contig flk = flankSeqs.get(rec.getReferenceIndex());
//			curAlnRec=new AlignmentRecord(rec, flk);
//			if(curAlnRec.refEnd-curAlnRec.refStart < (double)flkCov*flk.length()/100.0)
//				continue;
//			//not too far from the tip of read
//			else if(Math.min(-curAlnRec.readAlignmentEnd()+curAlnRec.readLength, curAlnRec.readAlignmentStart()) > (double)flk.length()*tipsPercentage/100.0){
//				continue;
//			}
//					
//			HashMap<Contig, AlignmentRecord> data;
//			if(!map.containsKey(curReadName) || map.get(curReadName).isEmpty()){
//				data=new HashMap<>();
//				data.put(flk, curAlnRec);
//				map.put(curReadName,  data);
//				
//			}else{
//				data=map.get(curReadName);
//				if(!data.containsKey(flk)){
//					data.put(flk, curAlnRec);
//				}else if(data.get(flk).score < curAlnRec.score){
//					data.replace(flk, curAlnRec);
//				}
//			}				
//		
//			
//
//		}// while
//		iter.close();
//		/**********************************************************************/
//		int totReadNum = map.keySet().size();
//		ArrayList<String>  	flank0 = new ArrayList<String>(),
//							flank1_0 = new ArrayList<String>(),
//							flank1_1 = new ArrayList<String>(),
//							flank2 = new ArrayList<String>();
//		map.keySet().stream().forEach(r->{
//			if(map.get(r).isEmpty())
//				flank0.add(r);
//			else if(map.get(r).keySet().size()==1){
//				if(map.get(r).containsKey(flankSeqs.get(0))){
//					AlignmentRecord alg=map.get(r).get(flankSeqs.get(0));
//					int totAlgFlank=alg.readAlignmentEnd()-alg.readAlignmentStart();
////					if( totAlgFlank > alg.readLength -insertLength)
////						return;
//					System.out.printf("Found read %s with only 1 flank sequence %s, insert length=%d\n", r, flankSeqs.get(0).getName(), (alg.readLength-totAlgFlank));
//					flank1_0.add(r);
//					
//				}
//				else if(map.get(r).containsKey(flankSeqs.get(1))){
//					AlignmentRecord alg=map.get(r).get(flankSeqs.get(1));
//					int totAlgFlank=alg.readAlignmentEnd()-alg.readAlignmentStart();
////					if( totAlgFlank > alg.readLength -insertLength)
////						return;
//					System.out.printf("Found read %s with only 1 flank sequence %s, insert length=%d\n", r, flankSeqs.get(1).getName(), (alg.readLength-totAlgFlank));
//					flank1_1.add(r);
//					
//				}
//				
//			}
//			else if(map.get(r).keySet().size()==2){
//				AlignmentRecord aln0 = map.get(r).get(flankSeqs.get(0)),
//								aln1 = map.get(r).get(flankSeqs.get(1));
//				if((aln0.readStart-aln1.readStart)*(aln0.readEnd-aln1.readStart) <= 0 
//				|| (aln0.readStart-aln1.readEnd)*(aln0.readEnd-aln1.readEnd) <= 0 )
//					return;
//				//the in-between must longer than insertLength
//				int totAlgFlank=aln0.readAlignmentEnd()-aln0.readAlignmentStart() + aln1.readAlignmentEnd()-aln1.readAlignmentStart();
//				if( totAlgFlank > aln1.readLength -insertLength)
//					return;
//				
//				System.out.printf("Found read %s with both flank sequences, insert length=%d\n", r, (aln1.readLength-totAlgFlank));
//				flank2.add(r);
//			}
//		});
//		System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
//		System.out.println( "Total number of reads: " + totReadNum);
//		System.out.println("Number of reads with 0 flank sequences: " + flank0.size());
//		flank0.stream().forEach(r->System.out.println(r));
//		System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
//		
//		System.out.printf("Number of reads with only 1 flank sequence %s: %d \n" , flankSeqs.get(0).getName(), flank1_0.size());
//		flank1_0.stream().forEach(r->System.out.println(r));
//		System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
//		
//		System.out.printf("Number of reads with only 1 flank sequence %s: %d \n" , flankSeqs.get(1).getName(), flank1_1.size());
//		flank1_1.stream().forEach(r->System.out.println(r));
//		System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
//		
//		System.out.println("Number of reads with 2 flank sequences: " + flank2.size());
//		flank2.stream().forEach(r->System.out.println(r));
//		System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
//	}


}

class FlankRecord{
	String readID;
	AlignmentRecord f0Rec, f1Rec, refRec;
	
	FlankRecord(String readID){
		this.readID=readID;
		f0Rec=f1Rec=refRec=null;
	}
	
	public String toString(){
		String retval = readID+"\t";
		if(f0Rec!=null)
			retval+=f0Rec.readStart+"\t"+f0Rec.readEnd+"\t";
		else retval+="-1\t-1\t";
		
		if(refRec!=null)
			retval+=refRec.readStart+"\t"+refRec.readEnd+"\t";
		else retval+="-1\t-1\t";
		
		if(f1Rec!=null)
			retval+=f1Rec.readStart+"\t"+f1Rec.readEnd+"\t";
		else retval+="-1\t-1\t";
		
		int trueFlank=-1; 
		if(f0Rec!=null && refRec!=null){
			int t0=(refRec.readStart-f0Rec.readStart)*(refRec.readStart-f0Rec.readEnd),
				t1=(refRec.readEnd-f0Rec.readStart)*(refRec.readEnd-f0Rec.readEnd);
			
			trueFlank=(t0<t1?refRec.readStart:refRec.readEnd);
		}
		return retval+trueFlank;
	} 
	
	public String printJunctionOnRef(){
		String retval = readID+"\t";
		
		int trueFlank=-1;
		if(refRec!=null){			
			if(f0Rec!=null){
				int t0=(refRec.readStart-f0Rec.readStart)*(refRec.readStart-f0Rec.readEnd),
					t1=(refRec.readEnd-f0Rec.readStart)*(refRec.readEnd-f0Rec.readEnd);
				
				trueFlank=(t0<t1?refRec.refStart:refRec.refEnd);
			}
			retval+=refRec.contig.getName()+"\t"+refRec.refStart+"\t"+refRec.refEnd+"\t"+refRec.strand+"\t";
			
		}else{
			retval+="NA\tNA\tNA\tNA\t";

		}
		

		return retval+(trueFlank==-1?"NA":trueFlank);
	}
}
