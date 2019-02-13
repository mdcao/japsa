package japsadev.bio.phylo;

import japsa.seq.FastaReader;
import japsa.seq.SequenceOutputStream;
import org.apache.commons.io.FileUtils;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.Arrays;

import static junit.framework.TestCase.assertTrue;

public class NCBITreeTest {
  private static final Logger LOG = LoggerFactory.getLogger(NCBITreeTest.class);

  
 /* public static void main(String[] args){
	  try{
		( new NCBITreeTest()).testPhylogeny();
	  }catch(Exception exc){
		  exc.printStackTrace();
	  }
  }*/
  @Test
  public  void testPhylogeny() throws Exception{
	  File commontree = new File("src/test/resources/commontree.txt.css.gz");
	  assertTrue(commontree.exists());
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
	  
  }

  
}