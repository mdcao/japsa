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
 * 10 Dec 2015 - Minh Duc Cao: Created                                        
 * 
 ****************************************************************************/
package japsadev.tools.work;

import japsa.seq.SequenceOutputStream;

import java.io.File;
import java.io.IOException;
import java.net.URL;
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

/**
 * @author minhduc
 *
 */
public class Test {

	/**
	 * @param args
	 * @throws ParserConfigurationException 
	 * @throws IOException 
	 * @throws SAXException 
	 */
	public static void main(String[] args) throws ParserConfigurationException, SAXException, IOException {
		
		String BASE = "/home/minhduc/MLST";
		
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
		System.out.println("============================");

		for (int temp = 0; temp < nList.getLength(); temp++){			
			//get species name
			Node node = nList.item(temp).getFirstChild();
			String species = node.getNodeValue();
			species = species.trim().replaceAll(" ", "_");
			File baseDir = new File(BASE+"/" + species);
			baseDir.mkdirs();

			SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(baseDir.getAbsolutePath()+"/logs");
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
					Resources.asByteSource(new URL(profileULR)).copyTo(Files.asByteSink(new File(baseDir.getAbsolutePath()+"/profile.dat")));
				}else if ("loci".equals(node.getNodeName())){
					NodeList lociList = node.getChildNodes();
					for (int i = 0; i< lociList.getLength();i++){
						Node locusNode = lociList.item(i);
						if ("locus".equals(locusNode.getNodeName ())){
							String geneName = locusNode.getChildNodes().item(0).getNodeValue().trim();
							String geneULR = locusNode.getChildNodes().item(1).getChildNodes().item(0).getNodeValue();
							Resources.asByteSource(new URL(geneULR)).copyTo(Files.asByteSink(new File(baseDir.getAbsolutePath()+"/" + geneName + ".fas")));							
							sos.print(geneName + " "+ geneULR + "\n" );							
						}
					}//for
				}//if

			}//for
			sos.close();
		}//for

	}

}
