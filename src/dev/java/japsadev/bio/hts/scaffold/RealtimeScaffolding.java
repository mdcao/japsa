package japsadev.bio.hts.scaffold;
import htsjdk.samtools.SAMRecord;

import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.io.IOException;
import java.lang.ProcessBuilder.Redirect;
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


	/**
	 * MDC tried to include BWA as part
	 * @param bamFile
	 * @param readNumber
	 * @param timeNumber
	 * @param minCov
	 * @param qual
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void scaffolding(String inFile, int readNumber, int timeNumber, double minCov, int qual, String format, String bwaExe, int bwaThread, String bwaIndex) 
			throws IOException, InterruptedException{
		scaffolder.setReadPeriod(readNumber);
		scaffolder.setTimePeriod(timeNumber * 1000);

		Logging.info("Scaffolding ready at " + new Date());

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
			Logging.info("Starting BWA at " + new Date());
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

			Logging.info("BWA started!");			

			//SequenceReader seqReader = SequenceReader.getReader(inFile);

			//SequenceOutputStream 
			//outStrs = new SequenceOutputStream(bwaProcess.getOutputStream());
			//Logging.info("set up output from bwa");

			//Start a new thread to feed the inFile into bwa input			
			//Thread thread = new Thread(){
			//	public void run(){
			//		Sequence seq;
			//		Alphabet dna = Alphabet.DNA16();
			//		try {
			//			Logging.info("Thread to feed bwa started");
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
		AlignmentRecord myRec = null;
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
			myRec = new AlignmentRecord(rec, graph.contigs.get(rec.getReferenceIndex()));

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
		AlignmentRecord myRec = null;
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
			myRec = new AlignmentRecord(rec, graph.contigs.get(rec.getReferenceIndex()));

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
				analysis();
			}

			try{
				//print for the last time if needed
				if(!ScaffoldGraph.updateGenome)
					scaffolding.graph.printSequences(true);
				
				outOS.close();
			}catch (Exception e){
				e.printStackTrace();
			}
		}

		@Override
		protected void analysis() {
			long step = (lastTime - startTime)/1000;//convert to second	
			ScaffoldGraph sg = scaffolding.graph;
			synchronized(sg){
				sg.connectBridges();
		
				try {
					// This function is for the sake of real-time annotation experiments being more readable
					//scaffolding.graph.printRT(scaffolding.currentBaseCount);
					sg.printSequences(ScaffoldGraph.updateGenome);
					outOS.print("Time |\tStep |\tRead count |\tBase count|\tNumber of scaffolds|\tCircular scaffolds |\tN50 | \tBreaks (maxlen)\n");
					outOS.print(timeNow + " |\t" + step + " |\t" + lastReadNumber + " |\t" + scaffolding.currentBaseCount + " |\t" + sg.getNumberOfContigs() 
							+ " |\t" + sg.getNumberOfCirculars() + " |\t" + sg.getN50() + " |\t" + sg.getGapsInfo());
	
					outOS.println();
					outOS.flush();
				} catch (IOException e) {
					e.printStackTrace();
				}			
			}
		}

		@Override
		protected int getCurrentRead() {
			// TODO Auto-generated method stub
			return scaffolding.currentReadCount;
		}

	}
}
