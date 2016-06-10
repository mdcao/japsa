package japsa.bio.hts.scaffold;
import htsjdk.samtools.SAMRecord;

import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;

import japsa.bio.np.RealtimeAnalysis;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.util.Logging;
//Simulate fastq realtime generator: jsa.np.timeEmulate -i <input> -output -
public class RealtimeScaffolding {
	RealtimeScaffolder scaffolder;
	public ScaffoldGraphDFS graph;
	int currentReadCount = 0;
	long currentBaseCount = 0;	
	
	public RealtimeScaffolding(String seqFile, String genesFile, String resistFile, String isFile, String oriFile, String output)throws IOException, InterruptedException{
		scaffolder = new RealtimeScaffolder(this, output);		
		graph = new ScaffoldGraphDFS(seqFile, genesFile, resistFile, isFile, oriFile);
	}
	
	public void scaffolding(String bamFile, int readNumber, int timeNumber, double minCov, int qual) 
							throws IOException, InterruptedException{
		scaffolder.setReadPeriod(readNumber);
		scaffolder.setTimePeriod(timeNumber * 1000);

		Logging.info("Scaffolding ready at " + new Date());

		//...
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		
		SamReader reader;
		if ("-".equals(bamFile))
			reader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
		else
			reader = SamReaderFactory.makeDefault().open(new File(bamFile));	

		SAMRecordIterator iter = reader.iterator();

		String readID = "";
		ReadFilling readFilling = null;
		ArrayList<AlignmentRecord> samList = null;// alignment record of the same read;		

		Thread thread = new Thread(scaffolder);
		thread.start();	
		while (iter.hasNext()) {
			SAMRecord rec = iter.next();

			if (rec.getReadUnmappedFlag() || rec.getMappingQuality() < qual){		
				if (!readID.equals(rec.getReadName())){
					readID = rec.getReadName();
					synchronized(this){
						currentReadCount ++;
						currentBaseCount += rec.getReadLength();
					}
				}
				continue;		
			}
			AlignmentRecord myRec = new AlignmentRecord(rec, graph.contigs.get(rec.getReferenceIndex()));

			if (readID.equals(myRec.readID)) {				
	
				if (myRec.useful){				
					for (AlignmentRecord s : samList) {
						if (s.useful){				
							//...update with synchronized
							synchronized(this.graph){
								graph.addBridge(readFilling, s, myRec, minCov);
								//Collections.sort(graph.bridgeList);
							}
						}
					}
				}
			} else {
				samList = new ArrayList<AlignmentRecord>();
				readID = myRec.readID;	
				readFilling = new ReadFilling(new Sequence(Alphabet.DNA5(), rec.getReadString(), "R" + readID), samList);	
				synchronized(this){
					currentReadCount ++;
					currentBaseCount += rec.getReadLength();
				}
			}
				
			samList.add(myRec);

		}// while
		scaffolder.stopWaiting();
		thread.join();
		iter.close();
		reader.close();		

	}
	public static class RealtimeScaffolder extends RealtimeAnalysis{
		RealtimeScaffolding scaffolding;
		public SequenceOutputStream outOS;
		RealtimeScaffolder(RealtimeScaffolding scf, String output)  throws IOException{
			scaffolding = scf;
			outOS = SequenceOutputStream.makeOutputStream(output);
			outOS.print("time\tstep\treads\tbases\tscaffolds\n");
		}

		@Override
		protected void close() {
			// TODO Auto-generated method stub
			try{
				outOS.close();
			}catch (Exception e){
				e.printStackTrace();
			}
		}

		@Override
		protected void analysis() {
			long step = (lastTime - startTime)/1000;//convert to second	
			scaffolding.graph.connectBridges();
			int scfCount = 0,
				cirCount = 0;
			for (int i = 0; i < scaffolding.graph.scaffolds.length;i++){
				if (scaffolding.graph.scaffolds[i].size() > 0){
					int len = scaffolding.graph.scaffolds[i].getLast().rightMost() - scaffolding.graph.scaffolds[i].getFirst().leftMost();
					if(scaffolding.graph.scaffolds[i].closeBridge != null){
						cirCount++;
						scfCount++;
						continue;
					}
					if (scaffolding.graph.contigs.get(i).head == i 
						&& !ScaffoldGraph.isRepeat(scaffolding.graph.contigs.get(i))
						&& len > ScaffoldGraph.maxRepeatLength)				
						scfCount++;
				}
			}
			try {
				//scaffolding.graph.printSequences();
				scaffolding.graph.printRT(scaffolding.currentBaseCount);
				outOS.print(timeNow + "\t" + step + "\t" + lastReadNumber + "\t" + scaffolding.currentBaseCount + "\t" + scfCount + "\t" + cirCount + "\t" + scaffolding.graph.getN50());
				outOS.println();
				outOS.flush();
			} catch (IOException e) {
				e.printStackTrace();
			}			
		}

		@Override
		protected int getCurrentRead() {
			// TODO Auto-generated method stub
			return scaffolding.currentReadCount;
		}
		
	}
}
