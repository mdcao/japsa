package japsadev.tools.makeCSS;

import static junit.framework.Assert.assertTrue;

import java.io.File;
import java.util.Arrays;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import japsadev.bio.phylo.AntibioticTree;
import japsadev.bio.phylo.CommonTree;
import japsadev.bio.phylo.GetTaxonID;
import japsadev.bio.phylo.NCBITree;
import japsadev.bio.phylo.NCBITreeTest;
import pal.tree.Tree;

/** commands for making CSS tree */

/* requires a directory with taxdump/names.dmp  from here
 * ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
 * as well as speciesIndex file
 */
public class CSSProcessCommand {
	  private static final Logger LOG = LoggerFactory.getLogger(CSSProcessCommand.class);

	public static void main(String[] args){
		//input files
		  File resistant_gene_list = new File("resistant_genes_list.txt");
		  File taxdump = new File("taxdump/names.dmp");
		  File speciesIndex = new File("speciesIndex");
		
		 
		  File treein = new File("commontree.txt");  //obtained from Step 3
		
		  
		  //output files
		  File taxon_file = new File("taxonid.txt");
		  File resistance_treeout = new File("resistancetree.txt.css");
		  File treeout = new File("commontree.txt.css");
		  File treeout_mod = new File("commontree.txt.css.mod");
		  
		/*Step 1 read resistance gene tree and color it */
		  if(resistant_gene_list.exists() && ! resistance_treeout.exists()){
			  makeResistanceTree(resistant_gene_list, resistance_treeout);
		  }
		  
		  //Step2 get taxon information
		  if(taxdump.exists() && speciesIndex.exists() && !taxon_file.exists()){
			  LOG.info("getting taxon information");
			  getTaxaForSpecies(taxdump, speciesIndex, taxon_file);
		  }
		  
		  /*Step -3 is manual step
		   the file taxon_file should be uploaded to 
		   https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi
		   and then download commontree.txt and put it in same directory  (selecting root and including unranked taxa)
		   */
		
		  /*Step 4 color species tree  */
		  if(treein.exists() && !treeout.exists()){
			  LOG.info("adding CSS to tree");
			  addCSSToTree(treein, treeout);
		  }
		
		  /*Step 5 place the lines from speciesIndex in the tree and add color  */
		  if(treeout.exists() && ! treeout_mod.exists()){
			  LOG.info("adding extra nodes from speciesIndex");
			  addExtraNodesFromSpeciesIndex(treeout, taxon_file, taxdump, speciesIndex, treeout_mod);
		  }
		  
		  /*Step 6 testing */
		  if(treeout_mod.exists()){
			  String totest = "Homo sapiens:Capnocytophaga canimorsus:Staphylococcus aureus:NC_023018.1";

			  LOG.info("testing");

			  test(treeout_mod, totest.split(":"));
		  }
		  
		  if(resistance_treeout.exists()){
			  LOG.info("testing");
			  String totest = "blaCTX:tetX";
			  test(resistance_treeout, totest.split(":"));
		  }
	}
	
	public static void getTaxaForSpecies(File taxdump, File speciesIndex, File output){
		  try{
			  File taxaToInclude = null;
			  GetTaxonID gid = new GetTaxonID(taxaToInclude, taxdump);
			  
			  if(true){
			 // gid.processGenBank(assembly_summary_inputs);
			  gid.process(speciesIndex);
			
			  gid.print(output);
			  }
		  }catch(Exception exc){
			  exc.printStackTrace();
		  }
	}
	
	public static void addCSSToTree(File treein, File treeout){
		try{ 
		 CommonTree trees = new NCBITree(treein);
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
	}
	
	public static void makeResistanceTree(File resistant_gene_list, File treeout){
		try{ 
		 CommonTree trees = new AntibioticTree(resistant_gene_list);
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
	}
	
	public static void addExtraNodesFromSpeciesIndex(File treein, File taxonid, File taxdump, File speciesIndex, File treeout){
		   try{
			   NCBITree t = new NCBITree(treein, taxonid, taxdump);//, new File(args[1]));
				t.addSpeciesIndex(speciesIndex);
				t.print(treeout);
		   }catch(Exception exc){
			   exc.printStackTrace();
		   }
	   }
	
	
	public static void test(File treein, String[] totest){
		
		try{
		NCBITree t = new NCBITree(treein);
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
