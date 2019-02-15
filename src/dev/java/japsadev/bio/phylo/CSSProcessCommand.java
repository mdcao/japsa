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
		  File nodesdmp= new File("taxdump/nodes.dmp");
		  File speciesIndex = new File("speciesIndex");
		
		 
		
		  
		  //output files
		  File treein = new File("commontree.txt");  //obtained from Step 3
		  File taxon_file = new File("taxonid.txt");
		  File resistance_treeout = new File("resistancetree.txt.css");
		  File treeout = new File("commontree.txt.css");
		  File treeout_mod = new File("commontree.txt.css.mod");
		  
		  GetTaxonID gid = null;
		  NCBITree trees_species = null ;
		  AntibioticTree trees_drugs = null;
		  
		/*Step 1 read resistance gene tree and color it */
		  if(resistant_gene_list.exists() && ! resistance_treeout.exists()){
			  trees_drugs  = makeResistanceTree(resistant_gene_list, resistance_treeout);
		  }
		  
		  //Step2 get taxon information
		  if(taxdump.exists() && speciesIndex.exists() && !taxon_file.exists()){
			  LOG.info("getting taxon information");
			 gid =  getTaxaForSpecies(taxdump, speciesIndex, taxon_file);
		  }
		//  if(true) return;
		  /*step -3 read tree */
		  if(taxdump.exists() && taxon_file.exists() && !treein.exists()){
			  LOG.info("making tree");
			  trees_species = readTaxaTree(gid, nodesdmp,  taxon_file,  taxdump, treein);
			 
		  }
		 
		  /*Step 4 color species tree  */
		  if(treein.exists() && !treeout.exists()){
			  LOG.info("adding CSS to tree");
			  addCSSToTree(trees_species, treein, treeout);
		  }
		
		  /*Step 5 place the lines from speciesIndex in the tree and add color  */
		  if(treeout.exists() && ! treeout_mod.exists()){
			  LOG.info("adding extra nodes from speciesIndex");
			  trees_species = addExtraNodesFromSpeciesIndex(trees_species, treeout, taxon_file, taxdump, speciesIndex, treeout_mod);
		  }
		  
		  /*Step 6 testing */
		  if(treeout_mod.exists()){
			  String totest = "NC_004355.1:AJ717516.1:NC_023018.1:Homo sapiens:Capnocytophaga canimorsus:Staphylococcus aureus:NC_023018.1:NC_002645.1";
			 
			  LOG.info("testing");

			  test(treeout_mod, totest.split(":"));
		

		  }
		  
		  if(resistance_treeout.exists()){
			  LOG.info("testing");
			  String totest = "blaCTX:tetX";
			  test(resistance_treeout, totest.split(":"));
		  }
	}
	
	public static GetTaxonID getTaxaForSpecies(File taxdump, File speciesIndex, File output){
		GetTaxonID gid = null;
		  try{
			  File taxaToInclude = null;
			gid = new GetTaxonID(taxaToInclude, taxdump);
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
	
	public static NCBITree readTaxaTree(GetTaxonID gid, File nodesdmp, File taxonid, File taxdump, File treeout){
		 NCBITree trees = null;
		try{ 
		 trees= 
				 gid == null ?  new NCBITree(nodesdmp, taxonid, taxdump) : new NCBITree(nodesdmp, gid);
				
		 
			trees.print(treeout);

		}catch(Exception exc){
			exc.printStackTrace();
			
		}
		return trees;
	}
	
	public static CommonTree addCSSToTree(CommonTree trees1, File treein, File treeout){
		CommonTree trees = null;
		try{ 
		  trees = //trees1==null ? 
				  new NCBITree(treein, null, null);  //need to re-read
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
	
	public static NCBITree addExtraNodesFromSpeciesIndex(NCBITree trees1, File treein, File taxonid, File taxdump, File speciesIndex, File treeout){
		NCBITree t  = null;
		try{
			   t = trees1 == null ? new NCBITree(treein, taxonid, taxdump) : trees1;//, new File(args[1]));
				t.addSpeciesIndex(speciesIndex);
				t.print(treeout);
		   }catch(Exception exc){
			   exc.printStackTrace();
		   }
		return t;
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
