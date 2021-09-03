package japsa.bio.hts.scaffold;


import java.io.File;
import java.io.IOException;
import java.lang.ProcessBuilder.Redirect;
import java.util.ArrayList;
import java.util.Date;

import japsa.bio.np.RealtimeAnalysis;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;


//Simulate fastq realtime generator: jsa.np.timeEmulate -i <input> -output -
public class RealtimeScaffolding {
	private static final Logger LOG = LoggerFactory.getLogger(RealtimeScaffolding.class);

	RealtimeScaffolder scaffolder;
	public ScaffoldGraphDFS graph;
	int currentReadCount = 0;
	long currentBaseCount = 0;	

	public RealtimeScaffolding(String seqFile, String genesFile, String resistFile, String isFile, String oriFile, String output)throws IOException, InterruptedException{
		scaffolder = new RealtimeScaffolder(this, output);		
		graph = new ScaffoldGraphDFS(seqFile, genesFile, resistFile, isFile, oriFile);
	}


	/**
	 * MDC tried to include BWA as part
	 * @param inFile
	 * @param readNumber
	 * @param timeNumber
	 * @param minCov
	 * @param qual
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void scaffolding2(String inFile, int readNumber, int timeNumber, double minCov, int qual, String format, String bwaExe, int bwaThread, String bwaIndex) 
			throws IOException, InterruptedException{
		scaffolder.setReadPeriod(readNumber);
		scaffolder.setTimePeriod(timeNumber * 1000);

		LOG.info("Scaffolding ready at " + new Date());

		//...
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader = null;

		Process bwaProcess = null;

		if (format.endsWith("am")){//bam or sam
			if ("-".equals(inFile))
				reader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
			else
				reader = SamReaderFactory.makeDefault().open(new File(inFile));	
		}else{
			LOG.info("Starting bwa  at " + new Date());
			ProcessBuilder pb = null;
			if ("-".equals(inFile)){
				pb = new ProcessBuilder(bwaExe, 
						"mem",
						"-t",
						"" + bwaThread,
						"-k11",
						"-W20",
						"-r10",
						"-A1",
						"-B1",
						"-O1",
						"-E1",
						"-L0",
						"-a",
						"-Y",
						"-K",
						"20000",
						bwaIndex,
						"-"
						).
						redirectInput(Redirect.INHERIT);
			}else{
				pb = new ProcessBuilder(bwaExe, 
						"mem",
						"-t",
						"" + bwaThread,
						"-k11",
						"-W20",
						"-r10",
						"-A1",
						"-B1",
						"-O1",
						"-E1",
						"-L0",
						"-a",
						"-Y",
						"-K",
						"20000",
						bwaIndex,
						inFile
						);
			}

			bwaProcess  = pb.redirectError(ProcessBuilder.Redirect.to(new File("/dev/null"))).start();

			LOG.info("bwa started x");

			//SequenceReader seqReader = SequenceReader.getReader(inFile);

			//SequenceOutputStream 
			//outStrs = new SequenceOutputStream(bwaProcess.getOutputStream());
			//LOG.info("set up output from bwa");

			//Start a new thread to feed the inFile into bwa input			
			//Thread thread = new Thread(){
			//	public void run(){
			//		Sequence seq;
			//		Alphabet dna = Alphabet.DNA16();
			//		try {
			//			LOG.info("Thread to feed bwa started");
			//			while ( (seq = seqReader.nextSequence(dna)) !=null){
			//				seq.writeFasta(outStrs);
			//			}
			//			outStrs.close();//as well as signaling
			//			seqReader.close();
			//		} catch (IOException e) {						//

			//		}finally{

			//		}
			//	}
			//};

			//thread.start();
			reader = SamReaderFactory.makeDefault().open(SamInputResource.of(bwaProcess.getInputStream()));

		}
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

		if (bwaProcess != null){
			bwaProcess.waitFor();
		}

	}

	public void scaffolding(String bamFile, int readNumber, int timeNumber, double minCov, int qual) 
			throws IOException, InterruptedException{
		scaffolder.setReadPeriod(readNumber);
		scaffolder.setTimePeriod(timeNumber * 1000);

		LOG.info("Scaffolding ready at " + new Date());

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
		}

		@Override
		protected void close() {
			//if SPAdes assembly graph is involved
			if(Contig.hasGraph()){
				ContigBridge.forceFilling();
				analysis(0);
			}

			try{
				outOS.close();
			}catch (Exception e){
				e.printStackTrace();
			}
		}

		@Override
		protected void analysis(int j) {
			if(j>0) throw new RuntimeException("!!");
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
				// This function is for the sake of real-time annotation experiments being more readable
				//scaffolding.graph.printRT(scaffolding.currentBaseCount);
				scaffolding.graph.printSequences();
				outOS.print("Time |\tStep |\tRead count |\tBase count|\tNumber of scaffolds|\tCircular scaffolds |\tN50 | \tBreaks (maxlen)\n");
				outOS.print(lastTime + " |\t" + step + " |\t" + lastReadNumber + " |\t" + scaffolding.currentBaseCount + " |\t" + scfCount 
						+ " |\t" + cirCount + " |\t" + scaffolding.graph.getN50() + " |\t" + scaffolding.graph.getGapsInfo());

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

		@Override
		protected int numSources() {
		return 1;
		}

	}
}
