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
		
		  File taxdump = new File("taxdump/names.dmp");
		  File speciesIndex = new File("speciesIndex");
		  File taxon_file = new File("taxonid.txt");
		  if(taxdump.exists() && speciesIndex.exists() && !taxon_file.exists()){
			  LOG.info("getting taxon information");
			  getTaxaForSpecies(taxdump, speciesIndex, taxon_file);
		  }
		 
		   /*the file taxonid.txt should be uploaded to 
		   https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi
		   and then download commontree.txt and put it in same directory  (selecting root and including unranked taxa)
		   */
		  File treein = new File("commontree.txt");
		  File treeout = new File("commontree.txt.css");
		  if(treein.exists() && !treeout.exists()){
			  LOG.info("adding CSS to tree");
			  addCSSToTree(treein, treeout);
		  }
		  File treeout_mod = new File("commontree.txt.css.mod");
		  if(treeout.exists() && ! treeout_mod.exists()){
			  LOG.info("adding extra nodes from speciesIndex");
			  addExtraNodesFromSpeciesIndex(treeout, taxon_file, taxdump, speciesIndex, treeout_mod);
		  }
		  String totest = "Homo sapiens:Capnocytophaga canimorsus:Staphylococcus aureus:NC_023018.1";
		  if(treeout_mod.exists()){
			  LOG.info("testing");

			  test(treeout_mod, totest.split(":"));
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
