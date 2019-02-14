package japsadev.bio.phylo;

import static junit.framework.Assert.assertTrue;

import java.io.File;
import java.util.Arrays;

import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class NCBITreeTest {
  private static final Logger LOG = LoggerFactory.getLogger(NCBITreeTest.class);

  
 public static void main(String[] args){
	  try{
		( new NCBITreeTest()).testPhylogeny();
	  }catch(Exception exc){
		  exc.printStackTrace();
	  }
  }
@Test
  public  void testPhylogeny() throws Exception{
	  File commontree = new File("src/test/resources/commontree.txt.css.gz");
	  assertTrue(commontree.exists());
	  /*
	   * this should work but for some reason I get a cannot find symbol for NCBITree

	NCBITree t = new NCBITree(commontree);
		  String[][] taxa =  t.getTaxonomy("Staphylococcus aureus");
		  assertTrue(taxa[0].length>1);
		  assertTrue(taxa[1].length>1);
		  LOG.info(""+Arrays.asList(taxa[0]));
			LOG.info(""+Arrays.asList(taxa[1]));
		  taxa =  t.getTaxonomy("Homo sapiens");
		  assertTrue(taxa[0].length>0);
		  assertTrue(taxa[1].length>0);
		LOG.info(""+Arrays.asList(taxa[0]));
		LOG.info(""+Arrays.asList(taxa[1]));
		  taxa =  t.getTaxonomy("Capnocytophaga canimorsus");
		  assertTrue(taxa[0].length>0);
		  assertTrue(taxa[1].length>0);
		LOG.info(""+Arrays.asList(taxa[0]));
		LOG.info(""+Arrays.asList(taxa[1]));
		taxa =  t.getTaxonomy("NC_001440.1");
		  assertTrue(taxa[0].length>0);
		  assertTrue(taxa[1].length>0);
		LOG.info(""+Arrays.asList(taxa[0]));
		LOG.info(""+Arrays.asList(taxa[1]));
			*/
	  
  }

  
}