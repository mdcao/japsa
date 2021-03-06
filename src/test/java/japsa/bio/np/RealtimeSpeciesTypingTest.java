package japsa.bio.np;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import junit.framework.TestCase;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;

public class RealtimeSpeciesTypingTest extends TestCase {
  private static final Logger LOG = LoggerFactory.getLogger(RealtimeSpeciesTypingTest.class);

  //public void typing(String bamFile, int readNumber, int timeNumber) throws IOException, InterruptedException
  @Test
  public void testTyping() throws Exception {
    RealtimeSpeciesTyping.JSON = true;
    int readNumber = 0;
    int timeNumber = 1;
    String idxFile = "src/test/resources/Bacterial_speciesIndex";
    String bamFile = "src/test/resources/NC_014923.1.sam";

    LOG.debug("species index file = " + idxFile);
    LOG.debug("bam file = " + bamFile);

    File b = new File(bamFile);
    LOG.info("bam input file exists? "+b.exists());

    File outFile = File.createTempFile("AnalysisResult_",".json");
    LOG.info("output tmp file = "+outFile.getAbsolutePath());

    BufferedReader br = new BufferedReader(new FileReader(new File(idxFile)));
    ByteArrayOutputStream os = new ByteArrayOutputStream();

    RealtimeSpeciesTyping typing;

    GsonBuilder gson_builder = new GsonBuilder();
    Gson gson = gson_builder.create();

    JsonElement element;
    JsonObject object;
    String jsonLine;


    //
    // BufferedReader / OutputStream tests
    //
    typing = new RealtimeSpeciesTyping(br, os);
    typing.setMinQual(1);
    typing.setTwoOnly(false);
    typing.typing(bamFile, readNumber, timeNumber);
    try {
      Thread.sleep(2000);
    } catch (InterruptedException e) {
      e.printStackTrace();
    }

    String output = os.toString("UTF-8");
    //LOG.info("json output = "+output);
    BufferedReader sr = new BufferedReader(new StringReader(output));

    jsonLine = sr.readLine();
    assertTrue(jsonLine.indexOf("others") > 0);

    jsonLine = sr.readLine();
    assertTrue(jsonLine.indexOf("Mesorhizobium ciceri") > 0);

    element = gson.fromJson(jsonLine, JsonElement.class);
    object = element.getAsJsonObject();
    assertTrue(true);

    //
    // Filename String tests
    //
    typing = new RealtimeSpeciesTyping(idxFile, outFile.getAbsolutePath());
    typing.setMinQual(1);
    typing.setTwoOnly(false);
    typing.typing(bamFile, readNumber, timeNumber);
    try {
      LOG.info("start sleep");
      Thread.sleep(5000);
      LOG.info("done sleeping");
    } catch (InterruptedException e) {
      e.printStackTrace();
    }

    BufferedReader outReader = new BufferedReader(new FileReader(outFile));

    assertTrue(outFile.exists());

    jsonLine = outReader.readLine();
    assertTrue(jsonLine.indexOf("others") > 0);

    jsonLine = outReader.readLine();
    LOG.info(jsonLine);
    assertTrue(jsonLine.indexOf("Mesorhizobium ciceri") > 0);

    element = gson.fromJson(jsonLine, JsonElement.class);
    object = element.getAsJsonObject();
    assertTrue(true);

  }
}
