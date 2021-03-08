package japsa.bio.phylo;


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/*import japsadev.bio.phylo.AntibioticTree;
#import japsadev.bio.phylo.CommonTree;
import japsadev.bio.phylo.GetTaxonID;
import japsadev.bio.phylo.NCBITree;
*/
import pal.tree.Node;
import pal.tree.Tree;

/** commands for making CSS tree */

/* requires a directory with taxdump/names.dmp  from here
 * ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
 * as well as speciesIndex file
 */
public class CSSProcessCommand {
	  private static final Logger LOG = LoggerFactory.getLogger(CSSProcessCommand.class);

	public static File getTree(File taxdir, File db, File speciesIndex, boolean addExtraNodes){
		 File resistant_gene_list = new File(taxdir,"resistant_genes_list.txt");
		//  File taxdir = new File("taxdump");
		 
		  if(!speciesIndex.exists()) throw new RuntimeException("!!");
	//	  if(!taxdir.exists()) taxdir =  new File("../taxdump");
		  if(!taxdir.exists()){
			  throw new RuntimeException("cannot find taxdump/ directory.  Please make a symbolic link into the "
			  		+ "working directory, or into parent directory");
			  
		  }
		  File taxdump = new File(taxdir,"names.dmp");
		  File nodesdmp= new File(taxdir,"nodes.dmp");
		//  File speciesIndex = new File("speciesIndex");
		  
		  //output files
		  File treein = new File(db,"commontree.txt");  //obtained from Step 3
		 // File taxon_file = new File("taxonid.txt");
		  File resistance_treeout = new File(db,"resistancetree.txt.css");
		  File treeout = new File(db,"commontree.txt.css");
		  File treeout_mod = new File(db,"commontree.txt.css.mod");
	
		/*Step 1 read resistance gene tree and color it */
		  if(resistant_gene_list.exists() && ! resistance_treeout.exists()){
			  makeResistanceTree(resistant_gene_list, resistance_treeout);
		  }
		  // Find the list of taxon to include
		  Set<Integer> taxon_set1 = new HashSet<Integer>();
		  try{
			  int col_ind = 2;
				  BufferedReader br1 = GetTaxonID.getBR(speciesIndex);
				  String st1;
				  while((st1 = br1.readLine())!=null){
					  String[] str = st1.split("\t");
					  taxon_set1.add(Integer.parseInt(str[col_ind]));
				  }
				  br1.close();
			  }catch(Exception exc){
				  exc.printStackTrace();
			  }
		 
		  try{  
		  
			
		  
		  //Step2 get taxon information
			/*
		  if(taxdump.exists() && speciesIndex.exists() && !taxon_file.exists()){
			  LOG.info("getting taxon information");
			 getTaxaForSpecies(gid, speciesIndex, taxon_file);
		  }*/
		//  if(true) return;
		  /*step -3 read tree */
		  if( !treein.exists() && !treeout.exists()){
			  LOG.info("making tree");
			int col_ind = 2;
			  GetTaxonID gid  = new GetTaxonID(taxdump, nodesdmp, taxon_set1);
			NCBITree trees = new NCBITree(gid);
			trees.print(treeout);
			 
		  }
		 
		  /*Step 4 color species tree  */
		  if(treein.exists() && !treeout.exists()){
			  
			  LOG.info("adding CSS to tree");
			  addCSSToTree(treein, treeout);
		  }
		
		  /*Step 5 place the lines from speciesIndex in the tree and add color */
		  if(treeout.exists() && ! treeout_mod.exists() && addExtraNodes){
			  
			  LOG.info("adding extra nodes from speciesIndex");
			 addExtraNodesFromSpeciesIndex( treeout,  speciesIndex, treeout_mod);
		  } 
		  
		  }catch(Exception exc){
			  exc.printStackTrace();
		  }
		//  if(treeout.exists()) return treeout;
		 return treeout;
		  
		
	}
	
	public static void test(String[] args){
		  File treeout_mod = new File("commontree.txt.css.mod");
		  File resistance_treeout = new File("resistancetree.txt.css");

		  /*Step 6 testing */
		  if(treeout_mod.exists()){
			  String totest = "NC_023018.1:Homo sapiens:Capnocytophaga canimorsus:Staphylococcus aureus:NC_023018.1:NC_002645.1";
			  LOG.info("testing");
			  test(treeout_mod, totest.split(":"));
		  }
		  
		  if(resistance_treeout.exists()){
			  LOG.info("testing");
			  String totest = "blaCTX:tetX";
			  test(resistance_treeout, totest.split(":"));
		  }
	  
	}
	
	public static GetTaxonID getTaxaForSpecies(GetTaxonID gid, File speciesIndex, File output){
		//GetTaxonID gid = null;
		  try{
			
			  if(true){
			 // gid.processGenBank(assembly_summary_inputs);
			  gid.process(speciesIndex);
			  gid.print(output);
			  }
		  }catch(Exception exc){
			  exc.printStackTrace();
		  }
		  return gid;
	}
	
	public static NCBITree readTaxaTree(GetTaxonID gid,File treeout, int col){
		 NCBITree trees = null;
		try{ 
		//	  gid.read(taxonFile, col); // this tells which nodes to include in tree
			trees = new NCBITree(gid);
			trees.print(treeout);

		}catch(Exception exc){
			exc.printStackTrace();
			
		}
		return trees;
	}
	
	public static  void color(Tree[] trees){
		try{
			for(int i=0; i<trees.length; i++){
				Tree tree = trees[i];
		ColorTree ct = new ColorTree(tree, true);
		ct.color();
		if(tree.getExternalNodeCount()>1){
		Node actualroot = tree.getRoot();
		Node root = ct.tree.getRoot();
		Object css = root.getIdentifier().getAttribute("css");
		while(root!=actualroot){
			System.err.println("coloring down");
			actualroot.getIdentifier().setAttribute("css", css);
			actualroot = actualroot.getChild(0);
		}
			}
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	public static  void colorEachLevel(Tree[] trees){
		try{
			for(int i=0; i<trees.length; i++){
				Tree tree = trees[i];
		ColorTree ct = new ColorTree(tree, false);
		ct.colorEachLevel();
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	public static  void colorRecursive(Tree[] trees, boolean even){
		try{
			for(int i=0; i<trees.length; i++){
				Tree tree = trees[i];
		ColorTree ct = new ColorTree(tree, false);
		ct.colorRecurisvely(even);
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	public static CommonTree addCSSToTree( File treein, File treeout){
		CommonTree trees = null;
		try{ 
		  trees = //trees1==null ? 
				  new NCBITree(treein, false);  //need to re-read
		  //: trees1;
		 Tree[] tree = trees.getTrees();
			for(int i=0; i<tree.length; i++){
				if(tree[i].getExternalNodeCount()>10){
					
					ColorTree ct = new ColorTree(tree[i], true);
					ct.color();
					Node actualroot = trees.roots.get(i);
					Node root = ct.tree.getRoot();
					Object css = root.getIdentifier().getAttribute("css");
					while(root!=actualroot){
						System.err.println("coloring down");
						actualroot.getIdentifier().setAttribute("css", css);
						actualroot = actualroot.getChild(0);
					}
				}
			}
			trees.print(treeout);

		}catch(Exception exc){
			exc.printStackTrace();
			
		}
		return trees;
	}
	
	public static AntibioticTree makeResistanceTree(File resistant_gene_list, File treeout){
		try{ 
		 AntibioticTree trees = new AntibioticTree(resistant_gene_list);
		 Tree[] tree = trees.getTrees();
			for(int i=0; i<tree.length; i++){
				if(tree[i].getExternalNodeCount()>10){
					ColorTree ct = new ColorTree(tree[i], true);
					ct.color();
				}
			}
			trees.print(treeout);
			return trees;
		}catch(Exception exc){
			exc.printStackTrace();
			
		}
		return null;
	}
	
	public static NCBITree addExtraNodesFromSpeciesIndex( File treein,  File speciesIndex, File treeout){
		NCBITree t  = null;
		try{
		//	System.err.println(gid1.getTaxa("Sclerophthora macrospora virus A"));
			   t = new NCBITree(treein, true) ;
			 //  t.gid  = gid1;//, new File(args[1]));
			   
				t.addSpeciesIndex(speciesIndex);
				if(treeout!=null) t.print(treeout);
		   }catch(Exception exc){
			   exc.printStackTrace();
		   }
		return t;
	   }
	
	
	public static void test(File treein, String[] totest){
		
		try{
		NCBITree t = new NCBITree(treein, false);
		for(int i=0; i<totest.length; i++){
		  String[][] taxa =  t.getTaxonomy(totest[i]);
		  LOG.info(totest[i]);
		  LOG.info(""+Arrays.asList(taxa[0]));
			LOG.info(""+Arrays.asList(taxa[1]));
		}
		 
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
}
