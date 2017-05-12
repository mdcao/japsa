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
 * 28/05/2014 - Minh Duc Cao: Created                                        
 ****************************************************************************/

package japsa.tools.bio.amra;


import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Date;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import com.google.common.io.Files;
import com.google.common.io.Resources;

import japsa.bio.amra.MLSTyping;
import japsa.bio.amra.MLSTyping.MLSType;
import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

//import org.slf4j.Logger;
//import org.slf4j.LoggerFactory;

/**
 * @author minhduc
 *
 */
@Deployable(
	scriptName = "jsa.amra.mlst",
	scriptDesc = "Multi-locus strain typing"
	)
public class MLSTCmd extends CommandLine{
    //private static final Logger LOG = LoggerFactory.getLogger(MLSTCmd.class);
	//CommandLine cmdLine;
	public MLSTCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("build", null, "Build the databases to this directory only");		
		addString("input", null, "Name of the genome file");
		addString("mlstScheme", null, "Folder containing the MLST scheme");
		addInt("top", 0, "If > 0, will provide top closest profile");

		addStdHelp();
	}

	public static void main(String [] args) throws IOException, InterruptedException, ParserConfigurationException, SAXException{
		MLSTCmd cmdLine = new MLSTCmd ();
		args = cmdLine.stdParseLine(args);

		String input = cmdLine.getStringVal("input");
		String buildPath = cmdLine.getStringVal("build");
		String mlstDir = cmdLine.getStringVal("mlstScheme");
		int top = cmdLine.getIntVal("top");		

		if (buildPath != null){
			build(buildPath);
			System.exit(0);
		}


		if (input == null){
			System.err.println("ERROR: The required param 'input' is not specified.");
			System.err.println(cmdLine.usageString());
			System.exit(-1);
		}
		if (mlstDir == null){
			System.err.println("ERROR: The required param 'mlstScheme' is not specified.");
			System.err.println(cmdLine.usageString());
			System.exit(-1);
		}



		//String blastn = cmdLine.getStringVal("blastn");		
		ArrayList<Sequence> seqs = FastaReader.readAll(input, Alphabet.DNA());
		if (top <= 0)
			System.out.println(MLSTyping.bestMlst(seqs, mlstDir));
		else{

			MLSTyping t = MLSTyping.topMlst(seqs, mlstDir);
			if (top > t.getProfiles().size())
				top = t.getProfiles().size();
			for (int i = 0; i < top; i++){
				MLSType p = t.getProfiles().get(i);
				System.out.println(p.getST() + " " + p.getScore());
			}

		}

	}

	public static void build(String base) throws ParserConfigurationException, SAXException, IOException {
		String SEP = File.separator;

		if (!base.endsWith(SEP))
			base = base + SEP;

		//Get Document Builder
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		DocumentBuilder builder = factory.newDocumentBuilder();

		//Build Document
		String link = "http://pubmlst.org/data/dbases.xml"; 
		Document document = builder.parse(link);		

		//Normalize the XML Structure; It's just too important !!
		document.getDocumentElement().normalize();

		//Here comes the root node
		Element root = document.getDocumentElement();
		System.out.println(root.getNodeName());

		//Get all employees
		NodeList nList = document.getElementsByTagName("species");		

		for (int temp = 0; temp < nList.getLength(); temp++){			
			//get species name
			Node node = nList.item(temp).getFirstChild();
			String species = node.getNodeValue();
			species = species.trim().replaceAll(" ", "_");
			File baseDir = new File(base + species);
			baseDir.mkdirs();

			SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(baseDir.getAbsolutePath()+ SEP + "logs");
			sos.print("Date " + new Date() + "\n");
			//go to mlst
			node = node.getNextSibling();

			//mlst field (blank)
			node = node.getFirstChild();

			//get to database
			node = node.getNextSibling();

			//first child = blank			
			NodeList childNodes =  node.getChildNodes();

			for (int x = 0; x < childNodes.getLength(); x++){
				node = childNodes.item(x);
				if ("url".equals(node.getNodeName()))				
					sos.print("URL: " + node.getFirstChild().getNodeValue () + "\n");
				else if ("retrieved".equals(node.getNodeName()))
					sos.print("retrieved: " + node.getFirstChild().getNodeValue () + "\n");
				else if ("profiles".equals(node.getNodeName())){
					String profileULR = node.getChildNodes().item(3).getFirstChild().getNodeValue();
					sos.print("Prifle: "+profileULR +"\n");
					Resources.asByteSource(new URL(profileULR)).copyTo(Files.asByteSink(new File(baseDir.getAbsolutePath()+SEP+"profile.dat")));
				}else if ("loci".equals(node.getNodeName())){
					NodeList lociList = node.getChildNodes();
					for (int i = 0; i< lociList.getLength();i++){
						Node locusNode = lociList.item(i);
						if ("locus".equals(locusNode.getNodeName ())){
							String geneName = locusNode.getChildNodes().item(0).getNodeValue().trim();
							String geneULR = locusNode.getChildNodes().item(1).getChildNodes().item(0).getNodeValue();
							Resources.asByteSource(new URL(geneULR)).copyTo(Files.asByteSink(new File(baseDir.getAbsolutePath()+SEP+geneName + ".fas")));							
							sos.print(geneName + " "+ geneULR + "\n" );							
						}
					}//for
				}//if

			}//for
			sos.close();
		}//for
	}
}
