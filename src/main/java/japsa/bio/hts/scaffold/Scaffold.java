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

/**************************     REVISION HISTORY    **************************
 * 20/12/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.bio.hts.scaffold;

import java.util.ArrayList;
import java.util.Collections;


public final class Scaffold extends ArrayList<Contig>{
		/**
		 * 
		 */
		private static final long serialVersionUID = -4310125261868862931L;
		private int start, end;
		public Scaffold(Contig myFContig){
			super();
			add(myFContig);
			start = 1;
			end = myFContig.length();
		}		
		
		public void addContig(Contig contig){
			add(contig);
		}
		
		public int getEnd(){
			return end;
		}
		
		public int getStart(){
			return start;
		}
		
		
		/**
		 * @param start the start to set
		 */
		public void setStart(int start) {
			this.start = start;
		}

		/**
		 * @param end the end to set
		 */
		public void setEnd(int end) {
			this.end = end;
		}

		public void viewNew(){

			System.out.println("========================== START =============================");			 
			int firstIndex = 0;
			int lastIndex = 0;

			for (Contig ctg: this){
				ctg.computePosition();				

				if (firstIndex > ctg.getStartPosition())
					firstIndex = ctg.getStartPosition();

				if (lastIndex > ctg.getEndPosition())
					lastIndex = ctg.getEndPosition();
			}//for
			Collections.sort(this);
			int previousIndex = firstIndex;
			for (int i = 0; i < size();i++){
				Contig ctg = get(i);
				System.out.printf("  contig %3d  ======" + (ctg.getRelDir() > 0?">":"<") + "%6d  %6d %s gaps = %6d\n",ctg.getIndex(), ctg.getStartPosition(),ctg.getEndPosition(), ctg.getName(), ctg.getStartPosition() - previousIndex );
				previousIndex = ctg.getEndPosition();
			}
			System.out.println("============================ END ===========================");
		}

	
/*************************************************
		public Sequence scaffoldSequence(ScaffoldVector [] vectors, ArrayList<Sequence> seqs){			
			//FIXME: this is very very slow
			ArrayList<Contig> ctgs = new ArrayList<Contig>(contigs.size());
			int firstIndex = 0;
			int lastIndex = 0;

			for (int i = 0; i < contigs.size();i++ ){
				Contig ctg = new Contig(contigs.get(i),null);				

				ctg.end = vectors[ctg.index].pos + vectors[ctg.index].dir * seqs.get(ctg.index).length();
				if (vectors[ctg.index].dir > 0){
					ctg.start = vectors[ctg.index].pos;
					ctg.end = vectors[ctg.index].pos + seqs.get(ctg.index).length();
				}else{
					ctg.end   = vectors[ctg.index].pos;
					ctg.start = vectors[ctg.index].pos - seqs.get(ctg.index).length();
				}
				ctgs.add(ctg);
				if (firstIndex > ctg.start)
					firstIndex = ctg.start;

				if (lastIndex > ctg.end)
					lastIndex = ctg.end;
			}//for
			Collections.sort(ctgs);
			//Sequence scaffoldSeq = new Sequence(Alphabet.DNA(),(lastIndex - firstIndex));
			int previousIndex = firstIndex;
			ByteArray bArray = new ByteArray();			

			for (int i = 0; i < ctgs.size();i++){
				Contig ctg = ctgs.get(i);
				Sequence seq = seqs.get(ctg.index);
				if (vectors[ctg.index].dir < 0)
					seq = Alphabet.DNA.complement(seq);
				//gaps
				while (previousIndex < ctg.start){
					bArray.add((byte) Alphabet.DNA.N);
					previousIndex ++;
				}
				for (int x = 0; x< seq.length();x++){
					bArray.add(seq.getBase(x));
				}
				previousIndex += seq.length();
			}
			return new Sequence(Alphabet.DNA(),bArray,"");
		}
		/*************************************************/
	}