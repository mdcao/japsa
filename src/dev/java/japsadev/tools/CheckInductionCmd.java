/*
 * Copyright (c) 2017  Minh Duc Cao (minhduc.cao@gmail.com).
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the names of the institutions nor the names of the contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/*                           Revision History                                
 * 11/01/2012 - Minh Duc Cao: Revised 
 * 01/01/2013 - Minh Duc Cao, revised                                       
 ****************************************************************************/

package japsadev.tools;

import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;


/**
 * @author Minh Duc Cao
 *
 */
@Deployable(
        scriptName = "jsa.dev.checkInduction",
        scriptDesc = "Sample script description"
)
public class CheckInductionCmd extends CommandLine{
    private static final Logger LOG = LoggerFactory.getLogger(CheckInductionCmd.class);
    public CheckInductionCmd(){
        super();
        Deployable annotation = getClass().getAnnotation(Deployable.class);
        setUsage(annotation.scriptName() + " [options]");
        setDesc(annotation.scriptDesc());

        addString("path",null, "path to the group");
        addString("list",null, "file contain list of samples");
        //addBoolean("reverse",false,"Reverse sort order");
        addStdHelp();
    }

    public static void main(String[] args) throws IOException {

        /*********************** Setting up script ****************************/
        CommandLine cmdLine = new CheckInductionCmd();
        args = cmdLine.stdParseLine(args);
        /**********************************************************************/

        String path = cmdLine.getStringVal("path");
        String list = cmdLine.getStringVal("list");

        BufferedReader reader = SequenceReader.openFile(list);
        String line = null;

        TreeSet<String> sampleSet = new TreeSet<String>();
        HashMap<String, Set<String>> geneMap = new HashMap<String, Set<String>>();
        HashSet<String> unionSet = new HashSet<String>();

        while ((line = reader.readLine())!=null){
            sampleSet.add(line.trim());
        }
        reader.close();


        Files.walk(Paths.get(path))
                .filter(Files::isRegularFile)
                .filter(p -> p.toString().endsWith("gff"))
                .filter(p -> sampleSet.contains(getSample(p)))
                .forEach(p-> {
                    Set<String> mySet = openGFF(p);
                    unionSet.addAll(mySet);
                    geneMap.put(getSample(p), mySet);
                    LOG.info("Union size " + unionSet.size());
                });
        ArrayList<String> unionList = new ArrayList<String>(unionSet);
        Collections.sort(unionList);

        HashSet<String> coreGene = new HashSet<String>();

        for (String gene:unionSet){
           // System.out.println(gene);
            boolean good = true;
            for (Set<String> mySet:geneMap.values()){
                if (!mySet.contains(gene)){
                    good = false;
                    break;
                }
            }//for
            if (good){
                coreGene.add(gene);
            }
        }



        System.out.print("gene");
        for (String gene:unionList) {
            if (coreGene.contains(gene))
                continue;
            System.out.print("\t" +gene);
        }
        System.out.println();

        for (String sample:sampleSet){
            System.out.print(sample);
            Set<String> mySet = geneMap.get(sample);
            for (String gene:unionList){
                if (coreGene.contains(gene))
                    continue;
                System.out.print("\t" + (mySet.contains(gene)?"Y":"N"));
            }
            System.out.println();
        }

        //System.out.println("====================================================");
        //for (String sample:sampleSet){
        //   System.out.print(sample);
        //    Set<String> mySet = geneMap.get(sample);
        //    for (String gene:unionList){
        //        System.out.print("\t" + (mySet.contains(gene)?"Y":"N"));
        //    }
        //    System.out.println();
       // }
    }

    public static String getSample(Path filePath){
        return filePath.getFileName().toString().replace(".gff","");
    }

    public static HashSet<String>  openGFF(Path fileName){

        HashSet<String> geneSet = new HashSet<String>();

        try {
            FileInputStream in = new FileInputStream(fileName.toFile());
            ArrayList<JapsaAnnotation> annoGFF = JapsaAnnotation.readMGFF(in, 0, 0, "CDS");
            in.close();

            for (JapsaAnnotation anno : annoGFF) {
                for (JapsaFeature f : anno.getFeatureList()) {
                    String desc = f.getDesc();
                    //String [] toks = desc.split(";");
                    int index = desc.indexOf(":UniProtKB:");
                    if (index >= 0) {
                        geneSet.add(desc.substring(index + 1, index + 17));
                        continue;
                    }
                    index = desc.indexOf(":CARD:");
                    if (index >= 0) {
                        geneSet.add(desc.substring(index + 1, index + 16));
                        continue;
                    }
                    index = desc.indexOf(":CLUSTERS:");
                    if (index >= 0) {
                        geneSet.add(desc.substring(index + 1, index + 18));
                        continue;//for
                    }
                    index = desc.indexOf(":Pfam:");
                    if (index >= 0) {
                        geneSet.add(desc.substring(index + 1, index + 15));
                        continue;
                    }
                    index = desc.indexOf(":HAMAP:");
                    if (index >= 0) {
                        geneSet.add(desc.substring(index + 1, index + 15));
                        continue;
                    }
                }//for
            }//for
        }catch (Exception e){
            e.printStackTrace();
        }

        LOG.info("Read " + geneSet.size() + " from " + fileName);
        return geneSet;
    }


}
/*RST*



 
  
  
  
  
  
  
  
  
*RST*/
  
