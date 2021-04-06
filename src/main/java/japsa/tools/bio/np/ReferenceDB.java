package japsa.tools.bio.np;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPOutputStream;

import japsa.bio.phylo.CSSProcessCommand;
import japsa.bio.phylo.GetTaxonID;
import japsa.bio.phylo.NCBITree;
import japsa.bio.phylo.Trie;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.tools.seq.SequenceUtils;
import pal.tree.Node;

public class ReferenceDB{
	// this class encapsulates everything required for a referenceDB
	public File refFile, modDB,speciesIndex;
	
	String dbs;
	public NCBITree tree;
	
	public HashMap<String, Integer> seq2Species = new HashMap<String, Integer>();
	public HashMap<String, Integer> species2Len = new HashMap<String, Integer>();
	public HashMap<String, Integer> species2Index = new HashMap<String, Integer>();
	//HashMap<String, SpeciesCount> species2Count = new HashMap<String, SpeciesCount>();
	public ArrayList<String> speciesList = new ArrayList<String>(); 
	//speciesFile should be local
	public void pretyping() throws IOException{
		Map<String, String> seq2Species1 = new HashMap<String, String>();
		readSpeciesIndex(speciesIndex.getAbsolutePath(), seq2Species1, species2Len, true);
		for(Iterator<String> it = seq2Species1.keySet().iterator(); it.hasNext();){
			String seq = it.next();
			String spec = seq2Species1.get(seq);
			Integer index = species2Index.get(spec);
			if(index==null){
				index = speciesList.size();
				speciesList.add(spec);
				species2Index.put(spec,index);
			}
			seq2Species.put(seq, index);
		}
	}
	
	public static void readSpeciesIndex(String indexFile, Map<String, String> seq2Species1, Map<String, Integer> seq2Len, boolean splitPlasmids)throws IOException{
		BufferedReader indexBufferedReader = GetTaxonID.getBR(new File(indexFile));
		String line = "";
		while ( (line = indexBufferedReader.readLine())!=null){
			if (line.startsWith("#"))
				continue;


			String sp=null,seq=null;
				
			String [] toks = line.split("\t");
			if(toks.length < 2){
				System.err.println("Illegal speciesIndex file!");
				System.exit(1);
			}
			//	boolean plasmid = splitPlasmids && toks[0].indexOf("plasmid")>=0;
			sp=toks[0].trim();
			//if(plasmid){
			//	sp = sp+".plasmid";
			//}
			seq=toks[1].split("\\s+")[0];
			seq2Len.put(seq,Integer.parseInt(toks[3]));
			//System.err.println("putting: "+sp);
			if (seq2Species1.put(seq, sp) != null)
				throw new RuntimeException("sequence " + seq +" presents multiple time");
//			else
//				LOG.info("==>adding " + seq + " to " + sp);
			
					
		}//while
	}
	
	String treef;
	
	public ReferenceDB(File dbdir)  throws IOException{
		this(dbdir,new File(dbdir.getParentFile()+"/taxdump"));
	}
	File taxaDir;
	public ReferenceDB(File dbdir, File taxaDir)  throws IOException{
		this.dbs  = dbdir.getName();
		this.taxaDir = taxaDir;
		//File dbPath = dbdir.getParentFile();
		//File dbdir = new File(dbPath+"/"+dbs);
		 refFile = new File(dbdir, "genomeDB.fna.gz");
		speciesIndex=Trie.getIndexFile(taxaDir, refFile);
		 treef = CSSProcessCommand.getTree(taxaDir,dbdir, speciesIndex, false).getAbsolutePath();
		boolean useTaxaAsSlug=false;
	//	String treef_mod = treef+".mod";
		tree = new NCBITree(new File(treef), useTaxaAsSlug);
		tree.addSpeciesIndex(speciesIndex);
		//tree.print(new File(treef_mod)); // this prints out to tree
//		if(true)System.exit(0);
		// treef = dbPath+"/"+dbs+"/"+ "commontree.txt.css.mod";
		 modDB = new File("./db");
		
		this.pretyping();
	}
	
	
	public ReferenceDB update(File speciesFile)  throws IOException{
		//long last_m = speciesFile.lastModified();
		Collection<String> targetSpecies = SequenceUtils.getReadList(speciesFile.getAbsolutePath(),false);
			 File newDB = new File(modDB, dbs+"_"+speciesFile.getName());//+"."+last_m);
			 newDB.mkdirs();
			  Alphabet alphabet = Alphabet.DNA();
				
				Set	<Node> targets = new HashSet<Node>();
				Iterator<String> it = targetSpecies.iterator();
				while(it.hasNext()){
					String nme = it.next();
					Node n = tree.getNode(nme);
					if(n!=null) targets.add(n);
				}
				SequenceReader reader = SequenceReader.getReader(refFile.getAbsolutePath());
				 
				File refFileOut = new File(newDB,"genomeDB.fna.gz");
				File indexFileOut = new File(newDB, "genomeDB.fna.gz.index.txt.gz" );
				SequenceOutputStream sos = new SequenceOutputStream(new GZIPOutputStream(new FileOutputStream(refFileOut)));
				int printed =0;
				 PrintWriter pw = new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(indexFileOut))));
				while (true){
					Sequence genome = reader.nextSequence(alphabet);
					if (genome == null)break;
					Integer ind = seq2Species.get(genome.getName());
					String nme = this.speciesList.get(ind);
				
					if(targetSpecies.contains(nme)){
						//seq2Species1.put(genome.getNa, value)
						genome.writeFasta(sos);
						Integer taxon =(Integer) this.tree.getNode(nme).getIdentifier().getAttribute("taxon");
						pw.println(nme+"\t"+genome.getName()+"\t"+taxon+"\t"+genome.length());
						printed++;
					}else{
						seq2Species.remove(genome.getName());
						if(false){
							Node node = tree.getNode(nme);
							Node parent = node;
							inner: while(parent!=null){
								if(targets.contains(node)){
									genome.writeFasta(sos);
									printed++;
									break inner;
								}
						}
						}
					}
				}
				if(printed==0) throw new RuntimeException("none extracted");

				sos.close();
				pw.close();
				return new ReferenceDB(newDB, this.taxaDir);
	}
	
}