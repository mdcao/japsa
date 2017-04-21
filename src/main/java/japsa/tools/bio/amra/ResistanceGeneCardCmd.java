/*****************************************************************************
 * Copyright (c) Minh Duc Cao, Monash Uni & UQ, All rights reserved.         *
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  * 
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 * 3. Neither the names of the institutions nor the names of the contributors*
 *    may be used to endorse or promote products derived from this software  *
 *    without specific prior written permission.                             *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS   *
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, *
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR    *
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR         *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 ****************************************************************************/

/*****************************************************************************
 *                           Revision History                                
 * 7 Sep 2015 - Minh Duc Cao: Created                                        
 *
 ****************************************************************************/
package japsa.tools.bio.amra;

import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import javax.json.Json;
import javax.json.JsonObject;
import javax.json.JsonReader;
import javax.json.JsonValue;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Collection;
import java.util.Set;


/**
 * @author minhduc
 */
@Deployable(
        scriptName = "jsa.amra.rescard",
        scriptDesc = "Finding resistance genes/classes in a genome using card database"
)
public class ResistanceGeneCardCmd extends CommandLine {

    public ResistanceGeneCardCmd() {
        super();
        Deployable annotation = getClass().getAnnotation(Deployable.class);
        setUsage(annotation.scriptName() + " [options]");
        setDesc(annotation.scriptDesc());

        //addString("input", null, "Name of the genome file",true);
        //addString("output", null, "Name of the output file",true);
        //addString("resDB", null, "Name of the resistance gene database",true);

        //addDouble("identity", 0.85, "Minimum identity");
        //addDouble("coverage", 0.85, "Minimum coverage of gene");

        addStdHelp();
    }

    /**
     * @param args
     */
    public static void main(String[] args) throws FileNotFoundException {
        CommandLine cmdLine = new ResistanceGeneCardCmd();
        args = cmdLine.stdParseLine(args);

        /**********************************************************************/
        //String input =  cmdLine.getStringVal("input");
        //String output = cmdLine.getStringVal("output");
        //String resDBPath = cmdLine.getStringVal("resDB");

        //double identity = cmdLine.getDoubleVal("identity");
        //double coverage = cmdLine.getDoubleVal("coverage");

        //script
        //curl -o card.tar.bz2  https://card.mcmaster.ca/download/0/broadstreet-v1.1.7.tar.gz
        //tar jxvf card.tar.bz2

        JsonReader reader = Json.createReader(new FileReader("card.json"));
        JsonObject sampleObject = reader.readObject();
        reader.close();


        Set<String> keys = sampleObject.keySet();

        for (String key : keys) {
            JsonValue jsonValue = sampleObject.get(key);

            if (jsonValue instanceof JsonObject) {
                try {
                    JsonObject jSeq = (JsonObject) jsonValue;
                    if (jSeq == null)
                        continue;

                    String desc = "";
                    String info;

                    info = jSeq.getString("ARO_accession");
                    if (info == null)
                        continue;
                    desc += "ARO" + info;

                    info = jSeq.getString("ARO_name");
                    if (info == null)
                        continue;
                    desc += " ~~~" + info;

                    info = jSeq.getString("ARO_description");
                    if (info == null)
                        continue;
                    desc += "~~~" + info;

                    //TODO: check this

                    jSeq = (JsonObject) jSeq.get("model_sequences");
                    jSeq = (JsonObject) jSeq.get("sequence");
                    Collection<JsonValue> values = jSeq.values();

                    for (JsonValue value : values) {
                        jSeq = (JsonObject) value;
                        String proteinSeq = jSeq.getJsonObject("protein_sequence").getString("sequence");
                        System.out.println(">" + desc + "\n" + proteinSeq.toString());

                        break;//for
                    }
                } catch (Exception e) {
                    System.err.println(e.getMessage() + '\n' + key);
                }


            }

        }

    }

}
