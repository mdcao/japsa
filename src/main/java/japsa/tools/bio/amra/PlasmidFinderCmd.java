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
 * 18/01/2017 - Minh Duc Cao: Created                                        
 ****************************************************************************/

package japsa.tools.bio.amra;


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Calendar;

import javax.xml.parsers.ParserConfigurationException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.xml.sax.SAXException;


import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
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

		addString("build", null, "Build the databases to this directory only");
		addString("input", null, "Name of the genome file");
		addString("plasmiddb", null, "Plasmid database");

		addStdHelp();
	}

	public static void main(String [] args) throws IOException, InterruptedException, ParserConfigurationException, SAXException{
		PlasmidFinderCmd cmdLine = new PlasmidFinderCmd ();
		args = cmdLine.stdParseLine(args);

		String build = cmdLine.getStringVal("build");
		String input = cmdLine.getStringVal("input");
		String plasmiddb = cmdLine.getStringVal("plasmiddb");


		if (build != null){
		    try {
                File buildFolder = new File(build);
                if (!buildFolder.isDirectory()) {
                    if (!buildFolder.mkdirs()) {
                        LOG.error("Cannot create folder " + build);
                        System.exit(1);
                    }
                }
                ProcessBuilder setup = new ProcessBuilder("curl", "-o", build + File.separator + "data.zip",
                        "--data", "folder=plasmidfinder&filename=plasmidfinder.zip","https://cge.cbs.dtu.dk/cge/download_data.php");
                Process process = setup.inheritIO().start();
                int status = process.waitFor();
                if (status != 0) {
                    LOG.error("Problem downloading the current database from plasmidfilder server");
                    System.exit(1);
                }

                setup = new ProcessBuilder("unzip", "-d", build, build + File.separator + "data.zip");
                process = setup.inheritIO().start();
                status = process.waitFor();
                if (status != 0) {
                    LOG.error("Problem unzip file " + build + File.separator + "data.zip");
                    System.exit(1);
                }
                Long timestamp = Calendar.getInstance().getTimeInMillis();
                Path actualFile = Paths.get(build + File.separator + "ORI" + timestamp + ".fasta");
                actualFile = Files.move(Paths.get(build + File.separator + "plasmid_database.fsa"), actualFile, REPLACE_EXISTING);
                Files.copy(actualFile, Paths.get(build + File.separator + "ORI.fasta"),REPLACE_EXISTING);


            }catch (IOException e){
		        e.printStackTrace();
                System.exit(1);
            }
		    System.exit(0);

		}
        //String blastn = cmdLine.getStringVal("blastn");
        ArrayList<Sequence> seqs = FastaReader.readAll(input, Alphabet.DNA());


        ProcessBuilder pb = new ProcessBuilder("blastn", "-subject", plasmiddb,
				"-query", input, "-outfmt", "6 qseqid qlen nident gaps mismatch");

		//
		Process process = pb.start();

		//SequenceOutputStream out = new SequenceOutputStream(process.getOutputStream());
		//for (Sequence seq:seqs){
		//	seq.writeFasta(out);
		//}
		//out.close();

		BufferedReader br = new BufferedReader(new InputStreamReader(process.getInputStream()));			

		String line;		
		while ((line = br.readLine()) != null) {
			System.out.println(line);
		}
		br.close();
		process.waitFor();		
	}
}
