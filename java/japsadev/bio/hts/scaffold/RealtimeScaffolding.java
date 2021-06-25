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
	public ScaffoldGraph graph;
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
	public void scaffoldingWithBWA(String inFile, int readNumber, int timeNumber, double minCov, int qual, String format, String bwaExe, int bwaThread, String bwaIndex) 
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

		String curReadID = "";
		Sequence fullSeq=null;
		ReadFilling readFilling = null;
		AlignmentRecord curAlnRec = null;
		ArrayList<AlignmentRecord> samList = null;// alignment record of the same read;		

		Thread thread = new Thread(scaffolder);
		thread.start();	
		
		SAMRecord rec;
		while (iter.hasNext()) {
			rec = iter.next();
			
			if (rec.getReadUnmappedFlag() || rec.getMappingQuality() < qual){		
				if (!curReadID.equals(rec.getReadName())){
					curReadID=rec.getReadName();
					fullSeq = new Sequence(Alphabet.DNA5(), rec.getReadString(), "R" + curReadID);
					synchronized(this){
						currentReadCount ++;
						currentBaseCount += rec.getReadLength();
					}
				}
				continue;		
			}
			
			Contig tmp = graph.contigs.get(rec.getReferenceIndex());
			if(tmp==null){
				Logging.error("Contig " + rec.getReferenceIndex() + " doesn't exist!");
				System.exit(1);
			}
			curAlnRec = new AlignmentRecord(rec, tmp);
//			System.out.println("Processing record of read " + rec.getReadName() + " and ref " + rec.getReferenceName() + (myRec.useful?": useful ":": useless ") + myRec);

			if (curReadID.equals(curAlnRec.readID) && fullSeq==null) { //meaning the readFilling is already set				

				if (curAlnRec.useful){				
					for (AlignmentRecord alnRec : samList) {
						if (alnRec.useful){				
							//...update with synchronized
							synchronized(this.graph){
								graph.addBridge(readFilling, alnRec, curAlnRec, minCov);
								//Collections.sort(graph.bridgeList);
							}
						}
					}
				}
			} else {
				samList = new ArrayList<AlignmentRecord>();
				curReadID = curAlnRec.readID;	
				if(fullSeq==null || !fullSeq.getName().equals("R"+curReadID))
//					fullSeq=rec.getReadBases();
					fullSeq = new Sequence(Alphabet.DNA5(), rec.getReadString(), "R" + curReadID);
				
				readFilling = new ReadFilling(fullSeq, samList);
				fullSeq=null;

				synchronized(this){
					currentReadCount ++;
					currentBaseCount += rec.getReadLength();
				}
			}

			samList.add(curAlnRec);
//			readFilling.print();

		}// while
		scaffolder.stopWaiting();
		thread.join();
		iter.close();
		reader.close();

		if (bwaProcess != null){
			bwaProcess.waitFor();
		}

	}

	/**
	 * SHN modified the default aligner to minimap2
	 * @param bamFile
	 * @param readNumber
	 * @param timeNumber
	 * @param minCov
	 * @param qual
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void scaffoldingWithMinimap2(String inFile, int readNumber, int timeNumber, double minCov, int qual, String format, String mm2Preset, int mm2Threads, String mm2Index) 
			throws IOException, InterruptedException{
		scaffolder.setReadPeriod(readNumber);
		scaffolder.setTimePeriod(timeNumber * 1000);

		Logging.info("Scaffolding ready at " + new Date());

		//...
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader = null;

		Process mm2Process = null;

		if (format.endsWith("am")){//bam or sam
			if ("-".equals(inFile))
				reader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
			else
				reader = SamReaderFactory.makeDefault().open(new File(inFile));	
		}else{
			Logging.info("Starting alignment by minimap2 at " + new Date());
			ProcessBuilder pb = null;
			if ("-".equals(inFile)){
				pb = new ProcessBuilder("minimap2", 
						"-t",
						"" + mm2Threads,
						"-ax",
						mm2Preset,
						"-I",
						"40g",
//						"-K",
//						"200M",
						mm2Index,
						"-"
						).
						redirectInput(Redirect.INHERIT);
			}else{
				pb = new ProcessBuilder("minimap2", 
						"-t",
						"" + mm2Threads,
						"-ax",
						mm2Preset,
						"-I",
						"40g",
//						"-K",
//						"200M",
						mm2Index,
						inFile
						);
			}

			mm2Process  = pb.redirectError(ProcessBuilder.Redirect.to(new File("/dev/null"))).start();

			Logging.info("minimap2 started!");			

			reader = SamReaderFactory.makeDefault().open(SamInputResource.of(mm2Process.getInputStream()));

		}
		SAMRecordIterator iter = reader.iterator();

		String curReadID = "";
		Sequence fullSeq=null;
		ReadFilling readFilling = null;
		AlignmentRecord curAlnRec = null;
		ArrayList<AlignmentRecord> samList = null;// alignment record of the same read;		

		Thread thread = new Thread(scaffolder);
		thread.start();	
		
		SAMRecord rec;
		while (iter.hasNext()) {
			rec = iter.next();
			
			if (rec.getReadUnmappedFlag() || rec.getMappingQuality() < qual){		
				if (!curReadID.equals(rec.getReadName())){
					curReadID=rec.getReadName();
					fullSeq = new Sequence(Alphabet.DNA5(), rec.getReadString(), "R" + curReadID);
					synchronized(this){
						currentReadCount ++;
						currentBaseCount += rec.getReadLength();
					}
				}
				continue;		
			}
			
			Contig tmp = graph.contigs.get(rec.getReferenceIndex());
			if(tmp==null){
				Logging.error("Contig " + rec.getReferenceIndex() + " doesn't exist!");
				System.exit(1);
			}
			curAlnRec = new AlignmentRecord(rec, tmp);
//			System.out.println("Processing record of read " + rec.getReadName() + " and ref " + rec.getReferenceName() + (myRec.useful?": useful ":": useless ") + myRec);

			if (curReadID.equals(curAlnRec.readID) && fullSeq==null) { //meaning the readFilling is already set				

				if (curAlnRec.useful){				
					for (AlignmentRecord alnRec : samList) {
						if (alnRec.useful){				
							//...update with synchronized
							synchronized(this.graph){
								graph.addBridge(readFilling, alnRec, curAlnRec, minCov);
								//Collections.sort(graph.bridgeList);
							}
						}
					}
				}
			} else {
				samList = new ArrayList<AlignmentRecord>();
				curReadID = curAlnRec.readID;	
				if(fullSeq==null || !fullSeq.getName().equals("R"+curReadID))
//					fullSeq=rec.getReadBases();
					fullSeq = new Sequence(Alphabet.DNA5(), rec.getReadString(), "R" + curReadID);
				
				readFilling = new ReadFilling(fullSeq, samList);
				fullSeq=null;

				synchronized(this){
					currentReadCount ++;
					currentBaseCount += rec.getReadLength();
				}
			}

			samList.add(curAlnRec);
//			readFilling.print();

		}// while
		scaffolder.stopWaiting();
		thread.join();
		iter.close();
		reader.close();

		if (mm2Process != null){
			mm2Process.waitFor();
		}

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
					scaffolding.graph.printSequences(true,false);
				
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
					scaffolding.graph.printRT(scaffolding.currentBaseCount);
					sg.printSequences(ScaffoldGraph.updateGenome,false);
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
