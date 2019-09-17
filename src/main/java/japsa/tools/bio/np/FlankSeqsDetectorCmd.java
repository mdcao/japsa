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
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.math3.ml.clustering.Cluster;
import org.apache.commons.math3.ml.clustering.DBSCANClusterer;
import org.apache.commons.math3.ml.clustering.DoublePoint;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

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
		addInt("distance", 3, "Distance for DBSCAN clustering algorithm.");
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
			distance = cmdLine.getIntVal("distance"),
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
		List<FlankRecord> records=new ArrayList<>();
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
				if(fr!=null){
//					System.out.println(fr.printJunctionOnRef());
					records.add(new FlankRecord(fr));
				}
					
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

				if(curAlnRecord.readAlignmentEnd()-curAlnRecord.readAlignmentStart() < insertLength){
					LOG.info("Ignore integration size too short: {}", curAlnRecord.toString());
					continue;
				}
				if(fr.refRec==null||fr.refRec.qual < curAlnRecord.qual)
					fr.refRec=curAlnRecord;
				else if(fr.refRec.qual == curAlnRecord.qual){
					LOG.info("Ignore read with confusing alignment to {}:\n {} \nvs\n {}", refName, curAlnRecord.toString(), fr.refRec);
					fr.refRec=null;
					continue;
				}
				
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
					if(fr.f0Rec==null||fr.f0Rec.qual < curAlnRecord.qual)
						fr.f0Rec=curAlnRecord;
					else if(fr.f0Rec.qual == curAlnRecord.qual){
						LOG.info("Ignore read with confusing alignment to {}:\n {} \nvs\n {}", refName, curAlnRecord.toString(), fr.refRec);
						fr.f0Rec=null;
						continue;
					}
					
				}else if(flankSeqs.size()>1 && flankSeqs.get(1).getName().equals(refName)){
					ctg=flankSeqs.get(1);
					curAlnRecord=new AlignmentRecord(curSAMRecord, ctg);
					if(curAlnRecord.refEnd-curAlnRecord.refStart < (double)flkCov*ctg.length()/100.0)
						continue;
					//not too far from the tip of read
					else if(Math.min(-curAlnRecord.readAlignmentEnd()+curAlnRecord.readLength, curAlnRecord.readAlignmentStart()) > (double)ctg.length()*tipsPercentage/100.0){
						continue;
					}
					if(fr.f1Rec==null||fr.f1Rec.qual < curAlnRecord.qual)
						fr.f1Rec=curAlnRecord;
					else if(fr.f1Rec.qual == curAlnRecord.qual){
						LOG.info("Ignore read with confusing alignment to {}:\n {} \nvs\n {}", refName, curAlnRecord.toString(), fr.refRec);
						fr.f1Rec=null;
						continue;
					}
				}else{
					LOG.info("Flank not found: {} != {}!", refName, flankSeqs.get(0).getName());
					System.exit(1);
				}
					
			}

			
		}// while
		iter.close();
		/**********************************************************************
		 * DBSCAN clustering
		 */
		
		List<DoublePoint> points = new ArrayList<DoublePoint>();
		for(int i=0;i<records.size();i++){
			FlankRecord frec = records.get(i);
			frec.calculateTF();
			if(frec.trueFlank>0)
				points.add(new DoublePoint(new double[]{frec.trueFlank, new Double(i)}));
		}
		
		DBSCANClusterer dbscan = new DBSCANClusterer(distance, 0, (a,b)->Math.abs(a[0]- b[0]));
		List<Cluster<DoublePoint>> cluster = dbscan.cluster(points);
		for(int i=0;i<cluster.size();i++) {
			Cluster<DoublePoint> c=cluster.get(i);
			TFCluster tfCluster=new TFCluster();
			for(DoublePoint p:c.getPoints()){ 
				records.get((int)p.getPoint()[1]).setTFCluster(tfCluster);
				tfCluster.add((int)p.getPoint()[0]);
			}

		}
		
		for(FlankRecord rec:records)
			System.out.println(rec.print());
	}



}
class TFCluster{
	private static final AtomicInteger count = new AtomicInteger(0); 
	private final int ID;

	List<Integer> values;
	
	TFCluster(){
		ID=count.incrementAndGet();
		values=new ArrayList<Integer>();
	}
	public int getId(){
		return ID;
	}
	public void add(int value){
		values.add(value);
	}
	public String toString(){
		//TODO: return stats as well
		DescriptiveStatistics stats = new DescriptiveStatistics(values.stream().mapToDouble(i->i).toArray());

		return getId()+"\t"+values.size()+"\t"+(int)stats.getMin()+"\t"+(int)stats.getMax()+"\t"+Math.round(stats.getMean())+"\t"+stats.getKurtosis();
	}
}
class FlankRecord{
	String readID;
	AlignmentRecord f0Rec, f1Rec, refRec;
	int trueFlank;
	TFCluster cluster;
	FlankRecord(String readID){
		this.readID=readID;
		f0Rec=f1Rec=refRec=null;
		trueFlank=-1;
		cluster=null;
	}
	FlankRecord(FlankRecord rec){
		readID=rec.readID;
		f0Rec=rec.f0Rec;
		f1Rec=rec.f1Rec;
		refRec=rec.refRec;
		trueFlank=rec.trueFlank;
		cluster=rec.cluster;
	}
	public void setTFCluster(TFCluster cluster){
		this.cluster=cluster;
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
		
		return retval+trueFlank;
	} 
	//call to calculate true flank
	public void calculateTF(){
		if(refRec!=null && f0Rec!=null){
			int t0=(refRec.readStart-f0Rec.readStart)*(refRec.readStart-f0Rec.readEnd),
					t1=(refRec.readEnd-f0Rec.readStart)*(refRec.readEnd-f0Rec.readEnd);
				
				trueFlank=(t0<t1?refRec.refStart:refRec.refEnd);
		}
	}
	public String print(){
		String retval = readID+"\t";
		
		if(refRec!=null)			
			retval+=refRec.contig.getName()+"\t"+refRec.refStart+"\t"+refRec.refEnd+"\t"+refRec.strand+"\t";		
		else
			retval+="-1\t-1\t-1\t-1\t";		

		return retval+""+trueFlank+"\t"+(cluster==null?"-1\t-1\t-1\t-1\t-1\t-1":cluster.toString());
	}
}
