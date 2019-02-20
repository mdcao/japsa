package japsa.bio.phylo;


import java.io.File;
import java.util.Arrays;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import pal.tree.Tree;

/** commands for making CSS tree */

/* requires a directory with taxdump/names.dmp  from here
 * ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
 * as well as speciesIndex file
 */
public class CSSProcessCommand {
	  private static final Logger LOG = LoggerFactory.getLogger(CSSProcessCommand.class);

	public static void main(String[] args){
		
		
		  File resistant_gene_list = new File("resistant_genes_list.txt");
		  File taxdir = new File("taxdump");
		  
		  if(!taxdir.exists()) taxdir =  new File("../taxdump");
		  if(!taxdir.exists()){
			  throw new RuntimeException("cannot find taxdump/ directory.  Please make a symbolic link into the "
			  		+ "working directory, or into parent directory");
			  
		  }
		  File taxdump = new File("taxdump/names.dmp");
		  File nodesdmp= new File("taxdump/nodes.dmp");
		  File speciesIndex = new File("speciesIndex");
		
		 
		
		  
		  //output files
		  File treein = new File("commontree.txt");  //obtained from Step 3
		  File taxon_file = new File("taxonid.txt");
		  File resistance_treeout = new File("resistancetree.txt.css");
		  File treeout = new File("commontree.txt.css");
		  File treeout_mod = new File("commontree.txt.css.mod");
		  
		 
		  //NCBITree trees_species = null ;
		//  AntibioticTree trees_drugs = null;
		  
		/*Step 1 read resistance gene tree and color it */
		  if(resistant_gene_list.exists() && ! resistance_treeout.exists()){
			  makeResistanceTree(resistant_gene_list, resistance_treeout);
		  }
		 
		  
		 
		  try{  
		  
			  GetTaxonID gid  = new GetTaxonID(taxdump, nodesdmp);
		  
		  //Step2 get taxon information
		  if(taxdump.exists() && speciesIndex.exists() && !taxon_file.exists()){
			  LOG.info("getting taxon information");
			 getTaxaForSpecies(gid, speciesIndex, taxon_file);
		  }
		//  if(true) return;
		  /*step -3 read tree */
		  if( taxon_file.exists() && !treein.exists()){
			  LOG.info("making tree");
			
			 readTaxaTree(gid,taxon_file,  treein);
			 
		  }
		 
		  /*Step 4 color species tree  */
		  if(treein.exists() && !treeout.exists()){
			  LOG.info("adding CSS to tree");
			  addCSSToTree(treein, treeout);
		  }
		
		  /*Step 5 place the lines from speciesIndex in the tree and add color  */
		  if(treeout.exists() && ! treeout_mod.exists()){
			  LOG.info("adding extra nodes from speciesIndex");
			 addExtraNodesFromSpeciesIndex( gid, treeout, taxon_file, taxdump, speciesIndex, treeout_mod);
		  }
		  
		  }catch(Exception exc){
			  exc.printStackTrace();
		  }
		 
		  
		
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
	
	public static NCBITree readTaxaTree(GetTaxonID gid, File taxonFile ,File treeout){
		 NCBITree trees = null;
		try{ 
			  gid.read(taxonFile); // this tells which nodes to include in tree
			 boolean cts =  gid.taxon_set.contains("191289");
			//System.err.println(gid.taxon_set.size());
			//System.err.println(gid.taxa2Sci.get("191289"));
			trees = new NCBITree(gid);
			
			trees.print(treeout);

		}catch(Exception exc){
			exc.printStackTrace();
			
		}
		return trees;
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
					ColorTree ct = new ColorTree(tree[i]);
					ct.color();
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
					ColorTree ct = new ColorTree(tree[i]);
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
	
	public static void addExtraNodesFromSpeciesIndex(GetTaxonID gid1, File treein, File taxonid, File taxdump, File speciesIndex, File treeout){
		NCBITree t  = null;
		try{
			System.err.println(gid1.getTaxa("Sclerophthora macrospora virus A"));
			   t = new NCBITree(treein, true) ;
			   t.gid  = gid1;//, new File(args[1]));
			   
				t.addSpeciesIndex(speciesIndex);
				t.print(treeout);
		   }catch(Exception exc){
			   exc.printStackTrace();
		   }
		//return t;
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
