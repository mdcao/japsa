package japsadev.bio.np.phage;

import java.io.File;
import java.io.IOException;
import java.lang.ProcessBuilder.Redirect;
import java.util.Date;

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

public class VectorSequence {
	static final int FLANKING=100;
	//Sequence plasmid;
	String plasmidFile; //plasmid fasta/q file that are already indexed by bwa
	int s5, e5,
		s3, e3;
	
	public VectorSequence(String seqFile, int e5, int s3) throws IOException{
//		SequenceReader reader = SequenceReader.getReader(seqFile);
//		plasmid = reader.nextSequence(Alphabet.DNA());
		plasmidFile=seqFile;
		this.e5=e5;
		this.s5=this.e5-FLANKING;
		this.s3=s3;
		this.e3=this.s3+FLANKING;
//		reader.close();
	}
	
	public Sequence extractInsertSequence(Sequence read) throws IOException{
		Sequence insert = null;

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		Process bwaProcess = null;
		ProcessBuilder pb = null;
		String bwaExe="bwa";
		
		read.writeFasta("-");
		
		pb = new ProcessBuilder(bwaExe, 
				"mem",
				"-t4",
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
				)
				.redirectInput(Redirect.INHERIT);
	
		bwaProcess  = pb.redirectError(ProcessBuilder.Redirect.to(new File("/dev/null"))).start();

		//Logging.info("bwa started x");	
		SamReader reader = SamReaderFactory.makeDefault().open(SamInputResource.of(bwaProcess.getInputStream()));
		String currentName = "";
		boolean direction = true;
		int readFront = 0, readBack = 0;
		SequenceOutputStream outFile = SequenceOutputStream.makeOutputStream("out.fasta");
		SAMRecordIterator iter = reader.iterator();
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
			if (record.getAlignmentStart() <= s5 && record.getAlignmentEnd() >= e5){
				int  [] refPositions = {s5, e5}; 
				int [] pos = HTSUtilities.positionsInRead(record, refPositions);
				if (pos[0] > 0 || pos[1] > 0){
					//System.out.printf("%5d %s %s %d %d %b\n", count, record.getReadName(), "FRONT", pos[0], pos[1], record.getReadNegativeStrandFlag());
					readFront = pos[0];
					if (readBack > 0 && direction == record.getReadNegativeStrandFlag()){
						if (readBack < readFront){
							System.err.printf("Bugger 1\n");							
						}else{
							String readSub = record.getReadString().substring(readFront,readBack);
							Sequence rs = new Sequence(Alphabet.DNA16(), readSub, record.getReadName());
							rs.writeFasta(outFile);													
						}
						direction = record.getReadNegativeStrandFlag();
					}
				}
					
			}
			
			if (record.getAlignmentStart() <= s3 && record.getAlignmentEnd() >= e3){
				int  [] refPositions = {s3, e3}; 
				int [] pos = HTSUtilities.positionsInRead(record, refPositions);
				if (pos[0] > 0 || pos[1] > 0){
					//System.out.printf("%5d %s %s %d %d %b\n", count, record.getReadName(), "FRONT", pos[0], pos[1], record.getReadNegativeStrandFlag());
					readBack = pos[0];
					if (readFront > 0 && direction == record.getReadNegativeStrandFlag()){
						if (readBack < readFront){
							System.err.printf("Bugger 2\n");							
						}else{
							String readSub = record.getReadString().substring(readFront,readBack);
							Sequence rs = new Sequence(Alphabet.DNA16(), readSub, record.getReadName());
							rs.writeFasta(outFile);													
						}
						direction = record.getReadNegativeStrandFlag();
					}
				}				
			}		
		}
		reader.close();
		return insert;
	}
}
