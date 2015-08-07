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
 * 20/09/2012 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.tools.hts.tr;

import japsa.seq.SequenceOutputStream;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;






/**
 * @author Minh Duc Cao (minhduc \dot cao \at gmail \dot com)
 *
 */
public class PEFragment implements Comparable<PEFragment>{
	/**
	 * Store a pair of mapped reads
	 */

	String readID = "";
	int refIndex = -1;

	//String refID = "";

	int start1, start2, len1, len2;
	int iSize;

	int direction;

	int fragmentBegin, fragmentEnd;
	int length;

	public int getStart(){
		return fragmentBegin;
	}

	public int getEnd(){
		return fragmentEnd;
	}
	/**
	 * Note that the length is not neccessarily the distance between the start
	 * and the end.
	 * @return
	 */
	public int getLength(){
		return length;
	}

	public int getISize(){
		return iSize;
	}

	public int getReferenceIndex(){
		return refIndex;
	}

	//public String getSeqID(){
	//	return refID;
	//}

	boolean done = false;


	PEFragment(int  refIdx, String rID) {		
		refIndex = refIdx;
		readID = rID;		
		//if(refID.length()<=1){
		//	refID = '0'+refID;
		//}
	}

	public void println(SequenceOutputStream out) throws IOException{		
		out.print(refIndex + "\t" + readID + "\t" + direction +  "\t" + fragmentBegin + "\t" + fragmentEnd + "\t" + length + "\t"+(iSize>0?iSize:-iSize)+'\n');
	}

	public static PEFragment readLine(String line){
		String [] toks   = line.split("\\t");

		PEFragment pair    = new PEFragment(Integer.parseInt(toks[0]), toks[1]);
		pair.direction   = Integer.parseInt(toks[2]);
		pair.fragmentBegin = Integer.parseInt(toks[3]);
		pair.fragmentEnd = Integer.parseInt(toks[4]);
		pair.length      = Integer.parseInt(toks[5]);
		pair.iSize       = Integer.parseInt(toks[6]);		

		return pair;		
	}

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(PEFragment o) {
		int c = refIndex - o.refIndex;
		if (c == 0)
			if (fragmentBegin == o.fragmentBegin)
				return fragmentEnd - o.fragmentEnd;
			else return fragmentBegin - o.fragmentBegin;			
		else
			return c;
	}

	public String toString(){
		return refIndex+"\t" + fragmentBegin + "\t" + fragmentEnd;
	}

	static class LinkedPEFragment extends PEFragment{

		LinkedPEFragment next, previous;	
		/**
		 * @param ref
		 * @param rID
		 */
		LinkedPEFragment(int ref, String rID) {
			super(ref, rID);
		}

		//public String toString(){
		//	return refIndex + "\t"+this.fragmentBegin + "\t" + this.fragmentEnd;
		//}

		public static LinkedPEFragment readLine(String line){
			String [] toks   = line.split("\\t");

			LinkedPEFragment pair    = new LinkedPEFragment(Integer.parseInt(toks[0]), toks[1]);
			pair.direction   = Integer.parseInt(toks[2]);
			pair.fragmentBegin = Integer.parseInt(toks[3]);
			pair.fragmentEnd = Integer.parseInt(toks[4]);
			pair.length      = Integer.parseInt(toks[5]);
			pair.iSize       = Integer.parseInt(toks[6]);	

			pair.previous = pair.next = null;

			return pair;		
		}
		/**
		 * Read from an input stream a list of fragments, and output them in sorted
		 * order. 
		 * 
		 * @param in
		 * @param out
		 * @param maxLength
		 * @throws IOException
		 */

		public static void read(BufferedReader in, SequenceOutputStream out, int maxLength) throws IOException{
			//Use insertion sort as the list is almost sorted
			LinkedPEFragment head = null, last = null;		
			String line = "";
			int pairIn = 0, pairOut = 0;
			ArrayList<String> desc = new ArrayList<String>();
			//Read the insert file to get all stats


			while ((line = in.readLine()) != null){

				//comments
				if (line.startsWith("#")){
					desc.add(line);				
					continue;
				}

				//headers: copy over for list of reference sequences
				if (line.startsWith("@")){
					out.print(line);
					out.print('\n');				
					continue;
				}

				line = line.trim();
				if (line.length() ==0 )
					continue;

				pairIn ++;			

				LinkedPEFragment pair = readLine(line);

				//Release if neccerasy: fragments from another reference sequence
				while (head != null && head.refIndex != pair.refIndex){
					head.println(out);
					pairOut ++;
					head = head.next;
					if (head != null)
						head.previous = null;
				}

				if (head == null){
					//if the queue is empty
					last = head = pair;				
				}else{
					LinkedPEFragment ptr = last;
					//search for the right place
					while (ptr.compareTo(pair) > 0 && ptr.previous != null){
						ptr = ptr.previous;
					}

					//find the right place				 
					if (ptr.compareTo(pair) > 0){//head
						pair.next = head;
						head.previous = pair;
						head = pair;				
					}else if (ptr.next == null){//the last
						last.next = pair;
						pair.previous = last;
						last = pair;
					}else{
						pair.next = ptr.next;
						ptr.next.previous = pair;
						pair.previous = ptr;
						ptr.next = pair;
					}//end of finding the right place
				}//head is not null


				while (pair.fragmentEnd - head.fragmentBegin > maxLength && head.next != null){
					head.println(out);
					pairOut ++;
					head = head.next;
					head.previous = null;
				}

				if (pairIn % Sam2FragmentSizeCmd.checkPoint == 0){
					Date date = new Date();
					System.err.println("No. of pairIn : " + pairIn + "No. of pairOut : " + pairOut +  " at  " + date.toString());
					System.gc();
				}

			}//while

			while (head != null){
				head.println(out);
				head = head.next;
			}
			for (String li: desc){		
				out.print(li);
				out.print('\n');
			}
		}	
	}
}
