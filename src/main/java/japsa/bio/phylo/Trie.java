package japsa.bio.phylo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import pal.misc.Identifier;
import pal.tree.Node;
import pal.tree.Tree;

public class Trie {
	
	
	static class TrieNode {
	    private final Map<Character, TrieNode> children = new HashMap<>();
	    Integer v = null;
	    private boolean endOfWord;

	    public Integer getVal() {
	    	if(v!=null) return v;
	    	else return children.values().iterator().next().getVal();
		}
	    
	    TrieNode getChild(char c){
	    	return this.children.get(c);
	    }
	    public TrieNode addChild(char c) {
	    	TrieNode tn = new TrieNode();
			this.children.put(c, tn);
			return tn;
		}
	    Map<Character, TrieNode> getChildren() {
	        return children;
	    }

	    boolean isEndOfWord() {
	        return endOfWord;
	    }

	    void setEndOfWord(boolean endOfWord) {
	        this.endOfWord = endOfWord;
	    }
	    void setInteger(int v){
	    	this.v = v;
	    }
		

		
	}
	int MAX_CNT = Integer.MAX_VALUE;
	public Trie(File in) throws FileNotFoundException, IOException{
		this.root = new TrieNode();
		InputStream is = new FileInputStream(in);
		InputStream is1 = in.getName().endsWith(".gz") ? new GZIPInputStream(is) : is;
		BufferedReader br = new BufferedReader(new InputStreamReader(is1));
		String st = "";
		int cnt=0;
		while(cnt < MAX_CNT && (st = br.readLine())!=null){
			String[] str = st.split("\t");
			insert(str[2], Integer.parseInt(str[0]));
			cnt=cnt+1;
		//	System.err.println(cnt);
		}
		br.close();
		System.err.println("finished");
	}
	
	public Trie(Tree t) throws FileNotFoundException, IOException{
		this.root = new TrieNode();
		int cnt=0;
		int len = t.getExternalNodeCount();
		for(int i=0; i<len; i++){
			Node n = t.getExternalNode(i);
			Identifier id = n.getIdentifier();
			insert(id.getName().replace("+-", "").trim(),(Integer) id.getAttribute("taxon"));
			cnt=cnt+1;
		//	System.err.println(cnt);
		}
		System.err.println("finished");
	}
	
	public char[] slug(String in){
		return in.toLowerCase().replace("\"", "").toCharArray();
	}
	//new File("genomeDB.fna.gz")
	public static Trie getIndexFile(File taxaDir, File fastaFile, String suffix, File outF) throws IOException{
		
		Trie tr = null;
		
			if(!fastaFile.exists()) throw new IOException("File does not exist "+fastaFile.getAbsolutePath());

			 tr = new Trie(new File(taxaDir,"names.dmp"));
			 tr.getIndex(fastaFile, outF);
		
			return tr;
	
	}
	public static int max_depth = 60;
	
	public void getIndex(File refFile, File outF) throws IOException{
			SequenceReader reader = SequenceReader.getReader(refFile.getAbsolutePath());
			Alphabet alphabet = Alphabet.DNA();
		
			PrintWriter pw = new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outF))));
			while (true){
				Sequence genome = reader.nextSequence(alphabet);
				if (genome == null)break;
				Integer sze = genome.length();
				String desc = genome.getDesc();
				int val = this.find(desc);
				int ind = desc.indexOf(' ');
				
			//	System.err.println(desc+" "+val+" "+foundString);
				pw.println(desc+"\t"+genome.getName()+"\t"+val+"\t"+sze); //+"\t"+foundString);
			}
			
		pw.close();
	//	return lenF;
	}
	    


		private final TrieNode root;

	    

	    String foundString = "";
	   public Integer find(String string) {
			char[] c = slug(string);
			  TrieNode current = root;
			  int i=0;
			  for ( i=0; i<c.length; i++) {
		        	TrieNode child = current.getChild(c[i]);
		        	if(child==null){
		        		break;
		        	}
		        	current = child;
		        }
			  foundString = string.substring(0,i);
			 // System.err.println("found: "+foundString);
			  return current.getVal();
		}
	    void insert(String word, int val) {
	        TrieNode current = root;
	        char[] c = slug(word);
	        for (int i=0; i<c.length && i<=max_depth; i++) {
	        	TrieNode child = current.getChild(c[i]);
	        	if(child==null){
	        		child = current.addChild(c[i]);
	        	}
	        	current = child;
	        }
	     //   System.err.println(word);
	        current.setEndOfWord(true);
	        current.setInteger(val);
	    }

	    boolean delete(String word) {
	        return delete(root, word, 0);
	    }

	    boolean containsNode(String word) {
	        TrieNode current = root;

	        for (int i = 0; i < word.length(); i++) {
	            char ch = word.charAt(i);
	            TrieNode node = current.getChildren().get(ch);
	            if (node == null) {
	                return false;
	            }
	            current = node;
	        }
	        return current.isEndOfWord();
	    }

	    boolean isEmpty() {
	        return root == null;
	    }

	    private boolean delete(TrieNode current, String word, int index) {
	        if (index == word.length()) {
	            if (!current.isEndOfWord()) {
	                return false;
	            }
	            current.setEndOfWord(false);
	            return current.getChildren().isEmpty();
	        }
	        char ch = word.charAt(index);
	        TrieNode node = current.getChildren().get(ch);
	        if (node == null) {
	            return false;
	        }
	        boolean shouldDeleteCurrentNode = delete(node, word, index + 1) && !node.isEndOfWord();

	        if (shouldDeleteCurrentNode) {
	            current.getChildren().remove(ch);
	            return current.getChildren().isEmpty();
	        }
	        return false;
	    }
	
}
