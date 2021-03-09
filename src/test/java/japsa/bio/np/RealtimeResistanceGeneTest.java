package japsa.bio.np;

import japsa.seq.FastaReader;
import japsa.seq.SequenceOutputStream;
import org.apache.commons.io.FileUtils;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReaderFactory;

import java.io.*;

import static junit.framework.TestCase.assertTrue;

public class RealtimeResistanceGeneTest {
  private static final Logger LOG = LoggerFactory.getLogger(RealtimeResistanceGeneTest.class);

  //	public RealtimeResistanceGene(int read, int time, String output, String resDB, String tmp) throws IOException{
  @Test
  public void testTyping() throws Exception {
    int readNumber = 0;
    int timeNumber = 1;
    Double scoreThreshold = 1.99D;
    String resDir  = "src/test/resources/resFinder";
    String recordPrefix = "JUNIT";

    File outFile = File.createTempFile("AnalysisResult_",".json");
    LOG.info("output tmp file = "+outFile.getAbsolutePath());

    File fastaFile = new File("src/test/resources/resFinder/DB.fasta");
    assertTrue(fastaFile.exists());
    InputStream fastaInputStream0 = new FileInputStream(fastaFile);

    File resDBFile = new File("src/test/resources/resFinder/geneList");
    assertTrue(resDBFile.exists());
    InputStream resDBInputStream0 = new FileInputStream(resDBFile);

    File bamFile = new File("src/test/resources/jsa361.sam");
    assertTrue(bamFile.exists());
    InputStream bamInputStream0 = new FileInputStream(bamFile);

    BufferedReader bamReader0 = new BufferedReader(new InputStreamReader(bamInputStream0));
    assertTrue(bamReader0.ready());

    ByteArrayOutputStream outStream = new ByteArrayOutputStream();
    RealtimeResistanceGene rg = new RealtimeResistanceGene(readNumber, timeNumber, outStream, resDBInputStream0, fastaInputStream0, recordPrefix);
    RealtimeResistanceGene.JSON = true;
    rg.setScoreThreshold(scoreThreshold);
    rg.typing(	SamReaderFactory.makeDefault().open(SamInputResource.of(bamInputStream0)).iterator());

    try {
      Thread.sleep(10000);
    } catch (InterruptedException e) {
      e.printStackTrace();
    }

    String jsonLine = null;

    BufferedReader br = new BufferedReader(new InputStreamReader(new ByteArrayInputStream(outStream.toByteArray())));
    jsonLine = br.readLine();
    LOG.info(jsonLine);
    assertTrue(jsonLine.indexOf("lastReadNumber\":0}") > 0);

    jsonLine = br.readLine();
    LOG.info(jsonLine);
    assertTrue(jsonLine.indexOf("lastReadNumber\":1}") > 0);

    jsonLine = br.readLine();
    LOG.info(jsonLine);
    assertTrue(jsonLine.indexOf("JSA_361") > 0);
  }
}