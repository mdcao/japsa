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

/****************************************************************************
 *                           Revision History                                
 * 25/06/2014 - Minh Duc Cao: Started
 *  
 ****************************************************************************/
package japsa.bio.alignment;

import java.io.IOException;
import java.util.Arrays;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.util.Logging;

/**
 * Implementation of multiple aligment of (long) short reads to the reference
 * genome. The alignment is implemeted as a linked list of sites, each of which
 * is an array of the corresponding positions from all sequences.
 * 
 * @author minhduc
 * 
 */
public class MultipleAlignment {
	static Alphabet alphabet = Alphabet.DNA5();
	// The head and the tail of the list
	NodeAlignment head = null;// , tail = null;

	// Number of sequences

	Sequence[] seqs;
	// NB: seqHash[0] = reference:

	int seqIndex = 0;// point to the next avaibale one

	public MultipleAlignment(int nSeq, Sequence ref) {
		seqs = new Sequence[nSeq];
		seqs[0] = ref;
		seqIndex = 1;
	}

	public void addRead(SAMRecord rec) {
		if (seqIndex >= seqs.length) {
			Logging.warn("More sequences are added " + seqIndex);
			return;
		}
		Sequence seq = new Sequence(alphabet, rec.getReadLength(),
				rec.getReadName());
		byte[] bases = rec.getReadBases();
		for (int i = 0; i < seq.length(); i++)
			seq.setBase(i, alphabet.byte2index(bases[i]));

		seqs[seqIndex] = seq;
		// seqIndex ++;

		// //////////////////////////////////////////////////////////
		int refBase = rec.getAlignmentStart();
		NodeAlignment current;
		// check if refBase already in the alignment
		if (head == null) {
			head = new NodeAlignment();
			head.site[0] = refBase;
			current = head;
		} else {
			current = head;
			if (current.site[0] < refBase) {
				// go forward to search
				while (true) {
					if (current.next != null)
						current = current.next;
					else {
						current.next = new NodeAlignment();
						current.next.prev = current;
						current = current.next;
						current.site[0] = current.prev.site[0] + 1;
					}

					if (current.site[0] == refBase)
						break;// while
				}// while
				// assert curent.site[0] == refBase
			}// if
			else {
				while (current.site[0] > refBase) {
					head = new NodeAlignment();
					head.next = current;
					current.prev = head;
					head.site[0] = current.site[0] - 1;
					current = head;
				}// while
			}// else
		}// else
		// assert: current.site[0] == refBase;

		int readBase = 1;
		int length;

		for (final CigarElement e : rec.getCigar().getCigarElements()) {
			switch (e.getOperator()) {
			case H:
				break; // ignore hard clips
			case P:
				break; // ignore pads
			case S:
				readBase += e.getLength();
				break; // soft clip read bases
			case N: // N ~ D
			case D:
				length = e.getLength();
				while (length > 0) {
					if (current.site[0] != 0) {
						if (current.site[0] != refBase) {
							Logging.exit("Fatal error " + refBase + " vs "
									+ current.site[0], 1);
						}
						length--;
						refBase++;
					}

					if (current.next == null) {
						current.next = new NodeAlignment();
						current.next.prev = current;
						current.next.site[0] = refBase;
						current = current.next;

					} else
						current = current.next;
				}// while

				break;// case
			case I:
				length = e.getLength();

				while (current.site[0] == 0 && length > 0) {
					length--;
					current.site[seqIndex] = readBase;
					readBase++;

					if (current.next == null) {
						current.next = new NodeAlignment();
						current.next.prev = current;
						current.next.site[0] = refBase;
						current = current.next;

					} else
						current = current.next;
				}// while
				while (length > 0) {
					NodeAlignment newNode = new NodeAlignment();

					newNode.prev = current.prev;
					if (newNode.prev != null)
						newNode.prev.next = newNode;
					else
						head = newNode;

					newNode.next = current;
					current.prev = newNode;

					length--;
					newNode.site[seqIndex] = readBase;
					readBase++;
				}// while
				break;

			case M:
			case EQ:
			case X:
				length = e.getLength();
				while (length > 0) {
					if (current.site[0] != 0) {
						if (current.site[0] != refBase) {
							Logging.exit("Fatal error " + refBase + " vs "
									+ current.site[0], 1);
						}
						length--;
						current.site[seqIndex] = readBase;
						readBase++;
						refBase++;
					}

					if (current.next == null) {
						current.next = new NodeAlignment();
						current.next.prev = current;
						current.next.site[0] = refBase;
						current = current.next;

					} else
						current = current.next;
				}
				break;
			default:
				throw new IllegalStateException(
						"Case statement didn't deal with cigar op: "
								+ e.getOperator());
			}// case
		}// for
		// /////////////////////////////////////////////////////////
		seqIndex++;
	}

	public void printAlignment() {
		for (int s = 0; s < seqIndex; s++) {
			System.out.printf("%20s : ", seqs[s].getName());
			NodeAlignment current = head;
			while (current != null) {
				if (current.site[s] == 0)
					System.out.print('-');
				else
					System.out.print(seqs[s].charAt(current.site[s] - 1));
				current = current.next;
			}// while
			System.out.println();
		}// for
		System.out.println("##########################");
	}

	public MultipleAlignment reduceAlignment(int f, int t) {

		MultipleAlignment reduce = new MultipleAlignment(2, seqs[0]);
		Sequence consensus = new Sequence(alphabet, seqs[1].length(),
				"consensus");
		reduce.seqs[1] = consensus;
		reduce.seqIndex = 2;

		int conIdx = 0;
		NodeAlignment reduceCurrent = null;

		NodeAlignment current = head, fN = null, tN = null;
		while (current != null) {
			if (current.site[0] == f)
				fN = current;

			if (current.site[0] > t) {
				tN = current;
				break;
			}
			current = current.next;
		}

		if (fN == null || tN == null)
			return null;

		// head = fN; head.prev = null;//gabbabe collection
		// tN.next = null;//gabbabe collection

		// so the aligment reduced to fN -> tN
		current = fN;

		int lastInx = alphabet.size();
		int[] votes = new int[lastInx + 1];

		while (current != null && current != tN.next) {
			Arrays.fill(votes, 0);
			// get the votes
			for (int s = 1; s < seqIndex; s++) {
				int loc = current.site[s] - 1;
				if (loc < 0)
					votes[lastInx]++;
				else
					votes[seqs[s].getBase(loc)]++;

			}
			// check the highest
			int best = lastInx;
			for (int i = best - 1; i >= 0; i--)
				if (votes[i] >= votes[best])
					best = i;

			if (best == lastInx) {// is a gap
				if (current.site[0] == 0) { // also a gap
					current = current.next;
					continue;// while
				} else {
					if (reduceCurrent == null) {
						reduce.head = reduceCurrent = new NodeAlignment();
					} else {
						reduceCurrent.next = new NodeAlignment();
						reduceCurrent.next.prev = reduceCurrent;
						reduceCurrent = reduceCurrent.next;
					}
					reduceCurrent.site[0] = current.site[0];
					reduceCurrent.site[1] = 0;// a gap
				}
			} else {// a char
				if (reduceCurrent == null) {
					reduce.head = reduceCurrent = new NodeAlignment();
				} else {
					reduceCurrent.next = new NodeAlignment();
					reduceCurrent.next.prev = reduceCurrent;
					reduceCurrent = reduceCurrent.next;
				}// else

				consensus.setBase(conIdx, (byte) best);
				conIdx++;
				reduceCurrent.site[0] = current.site[0];
				reduceCurrent.site[1] = conIdx;
			}// else
			current = current.next;
		}
		return reduce;
	}

	public void printAlignment(int f, int t) throws IOException {
		// Locate the first node (corresponds to f) and last node (t)
		NodeAlignment current = head, fN = null, tN = null;
		SequenceOutputStream os = SequenceOutputStream.makeOutputStream("-");
		while (current != null) {
			if (current.site[0] == f)
				fN = current;

			if (current.site[0] > t) {
				tN = current;
				break;
			}
			current = current.next;
		}

		if (fN == null || tN == null)
			return;

		current = fN;
		int pos = 0;
		while (current != null && current != tN) {
			if (current.site[0] != 0)
				pos = current.site[0];

			os.print(pos, 7);
			os.print(' ');

			for (int s = 0; s < seqIndex; s++) {
				if (current.site[s] == 0)
					os.print('-');
				else
					os.print(seqs[s].charAt(current.site[s] - 1));
			}
			current = current.next;
			os.print('\n');
		}
	}

	//Number of sequences printed out
	public int printFasta(int f, int t, String fileName) throws IOException {
		SequenceOutputStream os = SequenceOutputStream
				.makeOutputStream(fileName);
		//for (int s = 0; s < seqIndex; s++) {
		//	System.out.println(seqs[s].getName() + "  " + seqs[s].length());
		//}

		StringBuilder[] sbs = new StringBuilder[seqIndex];
		for (int i = 0; i < seqIndex; i++)
			sbs[i] = new StringBuilder(64);

		// Locate the first node (corresponds to f) and last node (t)
		NodeAlignment current = head;
		// Search for start
		while (current != null) {
			if (current.site[0] == f)
				break;
			current = current.next;
		}

		while (current != null && current.site[0] <= t) {
			for (int s = 0; s < seqIndex; s++) {
				if (current.site[s] != 0) {
					sbs[s].append(seqs[s].charAt(current.site[s] - 1));
				}
			}
			current = current.next;
		}

		for (int s = 1; s < seqIndex; s++) {
			if (sbs[s].length() > 6){
				os.print(">" + seqs[s].getName() + "_" + s);
				for (int i = 0; i < sbs[s].length(); i++) {
					if (i % 60 == 0)
						os.print("\n");
					os.print("" + sbs[s].charAt(i));
				}
				os.print("\n");
			}
		}
		os.close();

		return seqIndex;

	}

	/**
	 * A node in the list of
	 * 
	 * @author minhduc
	 * 
	 */
	class NodeAlignment {
		NodeAlignment next = null, prev = null;
		int[] site = new int[seqs.length];//
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
