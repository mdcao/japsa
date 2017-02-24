package japsadev.bio.np.phage;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.lang.ProcessBuilder.Redirect;
import java.util.ArrayList;
import java.util.Date;
import java.util.stream.Stream;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.HTSUtilities;
import japsa.util.Logging;

public class SequenceExtractor {
	static final int FLANKING=100;
	//Sequence plasmid;
	String plasmidFile; //plasmid fasta/q file that are already indexed by bwa
	int s5, e5,
		s3, e3;
	String currentName = "";
	String bwaExe = "bwa";
	
	public SequenceExtractor(String seqFile, String bwaExe, int e5, int s3) throws IOException{
		plasmidFile=seqFile;
		this.e5=e5;
		this.s5=this.e5-FLANKING;
		this.s3=s3;
		this.e3=this.s3+FLANKING;

	}
	
	public Sequence extractInsertSequence(Sequence read){
		Sequence insert = null;
		
		Thread bwa = new Thread(){
			public void run(){
		        Process bwaProcess = null;
		        try {
		    		ProcessBuilder pb  = new ProcessBuilder(bwaExe, 
		    				"mem",
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
		    				plasmidFile,
		    				"-"
		    				);
		    		bwaProcess = pb.redirectError(Redirect.to(new File("/dev/null"))).start();
		    		SequenceOutputStream toBWA = new SequenceOutputStream(bwaProcess.getOutputStream());
		    		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);

		    		read.writeFasta(toBWA);

		    		SamReader reader = SamReaderFactory.makeDefault().open(SamInputResource.of(bwaProcess.getInputStream())); // take the output from BWA
		    		SAMRecordIterator iter = reader.iterator();
		    		
		    		boolean direction = true;
		    		int readFront = 0, readBack = 0;
					bwaProcess.waitFor();

					while (iter.hasNext()){
						SAMRecord record = iter.next();
						if (record.getReadString().length() < 10){
							System.out.println("== " + record.getReadName());
							continue;//while
						}
						
						if (!record.getReadName().equals(currentName)){
							currentName = record.getReadName();
							direction = true;
							readFront = 0;
							readBack = 0;
						}
						System.out.println("Read " + currentName + " length " + record.getReadLength() + ": (" + record.getAlignmentStart() + ", " + record.getAlignmentEnd() + ")");
						
					}				
					toBWA.close();
					reader.close();
		        }catch (IOException e) {
		            e.printStackTrace();
		        }catch (InterruptedException e) {
		            return;
		        }
				
			}
		};
		bwa.start();

		
		return insert;
	}
	
	public static void main(String[] args){
		try {
			SequenceExtractor vector = new SequenceExtractor("/home/s.hoangnguyen/Projects/Phage/plasmid.fasta","bwa",1658,2735);
//			SequenceReader.readAll("/home/s.hoangnguyen/Projects/Phage/2d_1.fasta",Alphabet.DNA())
//			.stream()
//			.map(e -> vector.extractInsertSequence(e))
//			.forEach(e -> System.out.println(e==null));
			ArrayList<Sequence> list = SequenceReader.readAll("/home/s.hoangnguyen/Projects/Phage/test.fasta",Alphabet.DNA());
			for(Sequence e:list)
				vector.extractInsertSequence(e);
				
				
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
	}
	
 }
