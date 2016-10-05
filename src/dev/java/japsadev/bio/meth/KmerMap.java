package japsadev.bio.meth;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.FastqSequence;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.seq.nanopore.Fast5DetailReader;
import japsa.seq.nanopore.Fast5DetailReader.BaseCallAlignment2D;
import japsa.seq.nanopore.Fast5DetailReader.BaseCallEvents;
import japsa.seq.nanopore.Fast5NPReader.BaseCalledFastq;
import japsa.util.JapsaException;
import japsa.util.Logging;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class KmerMap {
	HashMap<String,ArrayList<Kmer>> methKmer = new HashMap<String,ArrayList<Kmer>>(),
			unmethKmer = new HashMap<String,ArrayList<Kmer>>(),
			finalList = new HashMap<String, ArrayList<Kmer>>();

	HashMap<String, Kmer> 	kmerList = new HashMap<String, Kmer>(); // kmer string -> Kmer object

	HashMap<String, SAMRecord> alignment = new HashMap<String, SAMRecord>(); // read name -> SAM record
	public KmerMap(String gffFile, String methID) throws IOException{
		BufferedReader gffReader = new BufferedReader(new FileReader(gffFile));
		String s;	
		Pattern tab = Pattern.compile("\t"),
				acgt = Pattern.compile("context=([ACGT]+);");
		for (s = gffReader.readLine(); null != s; s = gffReader.readLine()) {
			s = s.trim();
			if (s.length() > 0) {
				if (s.charAt(0) != '#'){
					String[] line = tab.split(s);
					String type=line[2].trim();
					if(type.equalsIgnoreCase(methID)){
						//String	name=line[0].split("\\|")[0];
						String	name=line[0];
						int	basePos = Integer.parseInt(line[3])-1;
						char strand = line[6].charAt(0);
						Matcher	des = acgt.matcher(line[8]);
						String context="";
						if(des.find())
							context=des.group(1);
						for(int i=16;i<21;i++){
							String kmer = context.substring(i, i+5);

							if(strand == '+'){
								Kmer aKmer = new Kmer(true, true, name, basePos+i-20);
								aKmer.setSeq(kmer);
								aKmer.setDesc(methID + ":" + (21-i) + "A");
								if(methKmer.containsKey(kmer))
									methKmer.get(kmer).add(aKmer);
								else{
									ArrayList<Kmer> list = new ArrayList<Kmer>();
									list.add(aKmer);
									methKmer.put(kmer, list);
								}
								kmerList.put(aKmer.toString(), aKmer);
							}
							else{
								Kmer aKmer = new Kmer(true, false, name, basePos+16-i);
								aKmer.setSeq(kmer);
								aKmer.setDesc(methID + ":" + (21-i) + "A");
								if(methKmer.containsKey(kmer))
									methKmer.get(kmer).add(aKmer);
								else{
									ArrayList<Kmer> list = new ArrayList<Kmer>();
									list.add(aKmer);
									methKmer.put(kmer, list);
								}
								kmerList.put(aKmer.toString(), aKmer);
							}
						}
					}
				}
			}
		}

		gffReader.close();

	}
	public void scanUnmethylate(String refFile) throws IOException{
		SequenceReader reader = SequenceReader.getReader(refFile);
		Sequence seq;
		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
			String name = seq.getName();
			for(int i=0; i< seq.length()-5; i++){
				Sequence kmer = seq.subSequence(i, i+5);
				if(methKmer.containsKey(kmer.toString())){			
					Kmer tmp = new Kmer(false, true, name, i);
					tmp.setSeq(kmer.toString());
					if(!methKmer.get(kmer.toString()).contains(tmp))
						if(unmethKmer.containsKey(kmer.toString()))
							unmethKmer.get(kmer.toString()).add(tmp);
						else{
							ArrayList<Kmer> list = new ArrayList<Kmer>();
							list.add(tmp);
							unmethKmer.put(kmer.toString(), list);
						}
					if(!kmerList.containsKey(tmp.toString()))
						kmerList.put(tmp.toString(), tmp);
				}

				kmer = Alphabet.DNA.complement(kmer);
				if(methKmer.containsKey(kmer.toString())){			
					Kmer tmp = new Kmer(false, false, name, i);
					tmp.setSeq(kmer.toString());
					if(!methKmer.get(kmer.toString()).contains(tmp))
						if(unmethKmer.containsKey(kmer.toString()))
							unmethKmer.get(kmer.toString()).add(tmp);
						else{
							ArrayList<Kmer> list = new ArrayList<Kmer>();
							list.add(tmp);
							unmethKmer.put(kmer.toString(), list);
						}
					if(!kmerList.containsKey(tmp.toString()))
						kmerList.put(tmp.toString(), tmp);
				}	
			}
		}
		reader.close();
	}
	public void scanAlignment(String samFile) throws IOException{
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);

		SamReader reader;
		if ("-".equals(samFile))
			reader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
		else
			reader = SamReaderFactory.makeDefault().open(new File(samFile));	

		SAMRecordIterator iter = reader.iterator();

		String readID = "";
		while (iter.hasNext()) {
			SAMRecord record = iter.next();
			if (record.getReadUnmappedFlag())
				continue;
			if (record.getMappingQuality() < 1) // or else?
				continue;
			readID = record.getReadName();
			if(!readID.contains("twodimentional"))
				continue;

			alignment.put(readID.substring(readID.indexOf("_")+1), record);

		}// while
		iter.close();
		reader.close();			
	}

	public static int move(String[] kmer, int begin, int step){
		int index = begin;
		while(step-- > 0){
			if(index >= kmer.length-1)
				return kmer.length-1;
			while(kmer[index].equals(kmer[++index]));
		}
		return index;
	}
	public void scanHDF5(String folderName) throws OutOfMemoryError{
		File mainFolder = new File(folderName);
		File [] fileList = mainFolder.listFiles();
		Logging.info("Reading in folder " + mainFolder.getAbsolutePath());
		if (fileList!=null){
			for (File f:fileList){
				//directory
				if (!f.isFile())
					continue;//for						

				if (!f.getName().endsWith("fast5"))
					continue;//for						
				String sPath = f.getAbsolutePath();

				try{
					Fast5DetailReader npReader = new Fast5DetailReader(sPath);
					npReader.readData();
					npReader.readFastq();
					ArrayList<BaseCalledFastq> seqList = npReader.getFastqList();
					FastqSequence fastq = null;
					for (BaseCalledFastq bfcq:seqList){
						if (bfcq.isTwoDim()){
							fastq = bfcq;
							break;
						}
					}					
					if(fastq==null){
						npReader.close();
						continue;
					}else{
						String readName = fastq.getName().split(" ")[0];
						SAMRecord record;
						if(alignment.containsKey(readName))
							record = alignment.get(readName);
						else{
							npReader.close();
							continue;
						}
						System.out.print("Reading "+ readName);
						// now read the whole HDF5 file
						npReader.readData();

						BaseCallAlignment2D align2d = npReader.getBcAlignment2D();

						BaseCallEvents temp = npReader.getBcTempEvents();
						if(align2d == null || temp == null){
							System.out.println("...ignored!");
							npReader.close();
							continue;
						}
						System.out.println();
						npReader.close();

						String [] kmer2d = align2d.getKmer();
						int curPos = 0;

						char strand = record.getReadNegativeStrandFlag()?'-':'+';
						String refID = record.getReferenceName();
						Cigar cigar = record.getCigar();

						int posOnRef = record.getAlignmentStart();
						for (final CigarElement e : cigar.getCigarElements()) {
							final int  length = e.getLength();
							switch (e.getOperator()) {
							case H :
							case S :					
							case P :
							case I :				
								curPos = move(kmer2d, curPos, length);
								break;
							case M ://match or mismatch				
							case EQ://match
							case X ://mismatch
								if(length < 5){
									posOnRef += length;
									curPos = move(kmer2d, curPos, length);
									break;
								}else{
									for(int i = 0; i < length-4; i++){
										posOnRef++; curPos = move(kmer2d, curPos, 1);
										String hash = new String(refID+":"+posOnRef+":"+strand);
										if(kmerList.containsKey(hash)){
											Kmer sigKmer = kmerList.get(hash);
											String seq = sigKmer.getSeq();
											if(kmer2d[curPos].equals(seq)){
												System.out.println("...hit " + hash + "\t" + kmer2d[curPos] + " vs. " + seq +  " " + sigKmer.getDesc());
												sigKmer.setSignal(temp.mean()[curPos], temp.stdv()[curPos]);

												if(finalList.containsKey(seq))
													finalList.get(seq).add(sigKmer);
												else{
													ArrayList<Kmer> list = new ArrayList<Kmer>();
													list.add(sigKmer);
													finalList.put(seq, list);
												}
											}else{
												System.out.println("...mismatch => ignored: " + hash);
											}
										}else{
											System.out.println("...no hit " + hash);
										}
									}
									break;
								}
							case D :
							case N :	
								posOnRef += length;
								break;								
							default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
							}//casse
						}//for	

					}
				}catch (JapsaException e){
					e.printStackTrace();
					continue;
				}catch (Exception e){
					Logging.error("Problem with reading " + sPath + ":" + e.getMessage());
					e.printStackTrace();
					continue;
				}	


			}//for			
		}//if
		else{
			Logging.info("Folder " + mainFolder.getAbsolutePath()  + " does not exist, are you sure this is the right folder?");					
		}

	}

	public static void statsHDF5(String folderName) throws OutOfMemoryError, IOException{
		BufferedWriter statFile = new BufferedWriter(new PrintWriter("kmerStats.out"));
		File mainFolder = new File(folderName);
		File [] fileList = mainFolder.listFiles();
		Logging.info("Reading in folder " + mainFolder.getAbsolutePath());
		if (fileList!=null){
			for (File f:fileList){
				//directory
				if (!f.isFile())
					continue;//for						

				if (!f.getName().endsWith("fast5"))
					continue;//for						
				String sPath = f.getAbsolutePath();

				try{
					Fast5DetailReader npReader = new Fast5DetailReader(sPath);
					npReader.readData();

					ArrayList<BaseCalledFastq> seqList = npReader.getFastqList();
					FastqSequence fastq = null;
					for (BaseCalledFastq bfcq:seqList){
						if (bfcq.isTwoDim()){
							fastq = bfcq;
							break;
						}
					}
					if(fastq==null){
						npReader.close();
						continue;
					}else{
						String readName = fastq.getName().split(" ")[0];
						System.out.print("Reading "+ readName);
						// now read the whole HDF5 file
						npReader.readData();

						BaseCallEvents 	temp = npReader.getBcTempEvents(),
								comp = npReader.getBcCompEvents();
						if(comp == null || temp == null){
							System.out.println("...ignored!");
							npReader.close();
							continue;
						}
						System.out.println();
						npReader.close();
						long [] move = temp.getMove();
						double[] 	mean = temp.mean(),
								stdv = temp.stdv();
						//long[]			modelLv = temp.modelLv();
						String[] modelState = temp.modelState();
						for(int i = 0; i < move.length; i++){
							if(move[i] == 1)
								statFile.write(modelState[i] + " " + "temp" + " " + mean[i] + " " + stdv[i] + "\n");
						}

						move = comp.getMove();
						mean = comp.mean();
						stdv = comp.stdv();
						//modelLv = comp.modelLv();
						modelState = comp.modelState();
						for(int i = 0; i < move.length; i++){
							if(move[i] == 1)
								statFile.write(modelState[i] + " " + "comp" + " " + mean[i] + " " + stdv[i]  + "\n");
						}
					}
				}catch (JapsaException e){
					e.printStackTrace();
					continue;
				}catch (Exception e){
					Logging.error("Problem with reading " + sPath + ":" + e.getMessage());
					e.printStackTrace();
					continue;
				}	


			}//for			
		}//if
		else{
			Logging.info("Folder " + mainFolder.getAbsolutePath()  + " does not exist, are you sure this is the right folder?");					
		}
		statFile.close();

	}

	public void print() throws IOException{
		BufferedWriter methFile = new BufferedWriter(new PrintWriter("methylated.out")),
				unmethFile = new BufferedWriter(new PrintWriter("unmethylated.out"));
		// Get a set of the entries
		Set<Entry<String,ArrayList<Kmer>>> set = finalList.entrySet();
		// Get an iterator
		Iterator<Entry<String, ArrayList<Kmer>>> ite = set.iterator();
		// Display elements
		while(ite.hasNext()) {
			Map.Entry<String,ArrayList<Kmer>> me = (Map.Entry<String,ArrayList<Kmer>>)ite.next();
			String kmer = me.getKey(),
					desc = "";
			ArrayList<Kmer> list = me.getValue();

			int index=0;
			while(index < list.size() && !list.get(index).isMeth) index++;
			desc = (index < list.size()) ? list.get(index).getDesc():""; 
			methFile.write(">"+ kmer + "\t" + desc + "\n");
			unmethFile.write(">"+ kmer + "\n");
			for(Kmer km:list)
				if(km.isMeth)
					methFile.write(km + "\t" + km.getMeanSignal() + "\t" + km.getStdvSignal() + "\n");
				else
					unmethFile.write(km + "\t" + km.getMeanSignal() + "\t" + km.getStdvSignal() + "\n");
		}
		methFile.close();
		unmethFile.close();

	}
}
