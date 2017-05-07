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
 * 11/01/2012 - Minh Duc Cao: Revised 
 * 01/01/2013 - Minh Duc Cao, revised                                       
 ****************************************************************************/

package japsadev.tools.work;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;


/**
 * @author Minh Duc Cao
 * 
 */
@Deployable(
	scriptName = "jsa.dev.cdhit",
	scriptDesc = "Sample script description"
	)
public class GetCDHitCmd extends CommandLine{
	public GetCDHitCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

        addString("sequence",null,"sequence",true);
		addString("input",null,"name",true);
        addString("output",null,"name",true);
		//addBoolean("reverse",false,"Reverse sort order");
		addStdHelp();		
	} 

	public static void main(String[] args) throws IOException {		

		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new GetCDHitCmd();
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/

        String sequence = cmdLine.getStringVal("sequence");
		String input = cmdLine.getStringVal("input");
        String output = cmdLine.getStringVal("output");

		ArrayList<Sequence> seqs = SequenceReader.readAll(sequence, Alphabet.DNA());
		HashMap<String, Sequence> map = new HashMap<String, Sequence>();

		for (Sequence seq:seqs){
			map.put(seq.getName(),seq);
		}

		BufferedReader reader = SequenceReader.openFile(input + ".clstr");
		ArrayList<Group> groups = new ArrayList<Group>();

		String line = "";

		Group group = null;
		while ((line = reader.readLine())!=null){
		    if (line.startsWith(">Cluste")){
		        if (group != null){
		            groups.add(group);

                }
		        group = new Group();
		        continue;
            }

            String [] toks = line.trim().split("\\s");
		    String name = toks[2].substring(1);
		    name = name.substring(0, name.length()-3);
		    group.count ++;
		    group.appendList(name);


		    Sequence s = map.get(name);
		    String desc = s.getDesc();
		    if (desc != null) {
                String[] descs = desc.split("\\s");
                try{
                    int reads = Integer.parseInt(descs[0]);
                    group.countRead +=  reads;
                }catch (Exception e){
                }
            }


		    if (toks[3].equals("*"))
		        group.seq = s;
		}//while
        if (group != null){
            groups.add(group);
        }

		reader.close();
        Collections.sort(groups);


        SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(output + ".fasta");
        SequenceOutputStream sosDat = SequenceOutputStream.makeOutputStream(output + ".dat");
        int index = 0;
        for (Group g: groups){
            index ++;
            Sequence s = g.seq;
            s.setName("group"+index);

            s.setDesc(g.count + " " + g.countRead);
            s.writeFasta(sos);

            sosDat.print("group" + index + " " + g.count + " " + g.countRead + " " + g.list);
            sosDat.println();

        }
        sos.close();
        sosDat.close();

    }

	static class Group implements Comparable<Group>{
	    Sequence seq = null;
	    int count = 0;
	    int countRead = 0;
	    String list = "";

	    public void appendList(String s){
	        if (list.length() > 0) list = list + ",";
	        list = list  + s;

        }

        @Override
        public int compareTo(Group o) {
            return o.count - this.count;
        }
    }
}

