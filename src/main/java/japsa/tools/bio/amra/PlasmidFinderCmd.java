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

/*                           Revision History                                
 * 18/03/2017 - Minh Duc Cao: Created
 ****************************************************************************/

package japsa.tools.bio.amra;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Calendar;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import static java.nio.file.StandardCopyOption.REPLACE_EXISTING;

/**
 * Identify plasmids from an assembly
 * @author minhduc
 *
 */
@Deployable(
        scriptName = "jsa.amra.plasmidfinder",
        scriptDesc = "Multi-locus strain typing"
)
public class PlasmidFinderCmd extends CommandLine{
    private static final Logger LOG = LoggerFactory.getLogger(PlasmidFinderCmd.class);

    //CommandLine cmdLine;
    public PlasmidFinderCmd(){
        super();
        Deployable annotation = getClass().getAnnotation(Deployable.class);
        setUsage(annotation.scriptName() + " [options]");
        setDesc(annotation.scriptDesc());

        addBoolean("build", false, "To build the latest from plasmidFinder");
        addString("input", null, "Name of the genome file");
        addString("plasmiddb", null, "Folder of the the plasmid database",true);

        addStdHelp();
    }

    public static void main(String [] args) throws IOException, InterruptedException{
        PlasmidFinderCmd cmdLine = new PlasmidFinderCmd ();
        args = cmdLine.stdParseLine(args);

        boolean build = cmdLine.getBooleanVal("build");
        String input = cmdLine.getStringVal("input");
        String plasmidFolder = cmdLine.getStringVal("plasmiddb");

        if (build){
            try {
                File buildFolder = new File(plasmidFolder);
                if (!buildFolder.isDirectory()) {
                    if (!buildFolder.mkdirs()) {
                        LOG.error("Cannot create folder " + plasmidFolder);
                        System.exit(1);
                    }
                }
                ProcessBuilder setup = new ProcessBuilder("curl", "-o", plasmidFolder + File.separator + "data.zip",
                        "--data", "folder=plasmidfinder&filename=plasmidfinder.zip","https://cge.cbs.dtu.dk/cge/download_data.php");
                Process process = setup.inheritIO().start();
                int status = process.waitFor();

                if (status != 0) {
                    LOG.error("Problem downloading the current database from plasmidfilder server");
                    System.exit(1);
                }

                setup = new ProcessBuilder("unzip", "-o","-d", plasmidFolder, plasmidFolder + File.separator + "data.zip");
                process = setup.inheritIO().start();
                status = process.waitFor();
                if (status != 0) {
                    LOG.error("Problem unzip file " + plasmidFolder + File.separator + "data.zip");
                    System.exit(1);
                }
                Long timestamp = Calendar.getInstance().getTimeInMillis();
                Path actualFile = Paths.get(plasmidFolder + File.separator + "ORI" + timestamp + ".fasta");
                actualFile = Files.move(Paths.get(plasmidFolder + File.separator + "plasmid_database.fsa"), actualFile, REPLACE_EXISTING);
                Files.copy(actualFile, Paths.get(plasmidFolder + File.separator + "ORI.fasta"),REPLACE_EXISTING);
            }catch (IOException e){
                e.printStackTrace();
                System.exit(1);
            }
            System.exit(0);
        }

        if (input == null){
            System.err.println("Please specify an input file");
            System.exit(1);
        }

        ProcessBuilder pb = new ProcessBuilder("blastn", "-subject", plasmidFolder + File.separator + "ORI.fasta",
                "-query", input, "-outfmt", "6 qseqid qlen qstart qend sseqid slen sstart send length frames pident nident gaps mismatch score bitscore");
        //qseqid qlen qstart qend sseqid slen sstart send length frames pident nident gaps mismatch score bitscore
        Process process = pb.start();
        BufferedReader br = new BufferedReader(new InputStreamReader(process.getInputStream()));
        String line;
        while ((line = br.readLine()) != null) {
            String [] toks = line.trim().split("\t");
            double oLength = Double.parseDouble(toks[5]);
            double oCov = Math.abs(Double.parseDouble(toks[7])- Double.parseDouble(toks[6])) + 1;

            double ratio = oCov / oLength;
            double identity = Double.parseDouble(toks[10]) / 10;

            if (ratio > 0.9 && identity > 0.9){
                System.out.println(toks[0] + "\t" + toks[4] + "\t" + ratio + "\t" + identity);
            }

        }
        br.close();
        process.waitFor();
    }
}
