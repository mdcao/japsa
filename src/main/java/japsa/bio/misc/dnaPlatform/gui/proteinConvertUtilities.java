/******************************************************************************
 * Copyright (C) 2006-2010 Minh Duc Cao                                        *
 *                                                                             *
 * This program is free software; you can redistribute it and/or modify it     *
 * under the terms of the GNU General Public License as published by the Free  *
 * Software Foundation; either version 2 of the License, or (at your option)   *
 * any later version. This program is distributed in the hope that it will be  *
 * useful, but WITHOUT ANY WARRANTY; without even the implied warranty of      *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General    *
 * Public License for more details.                                            *
 *                                                                             *
 * You should have received a copy of the GNU General Public License along with*
 * this program; if not, write to the Free Software  Foundation, Inc.,         *
 * 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.                   *
 ******************************************************************************/

//This class is written by Hoang Anh Nguyen and is subsequently modified and maintained
//by Minh Duc Cao

package japsa.bio.misc.dnaPlatform.gui;

/**
 * 
 * @author hoangnguyen
 */
public class proteinConvertUtilities {

	public static int A = 0;
	public static int C = 1;
	public static int G = 2;
	public static int T = 3;
	public static char STOP_CONDON = '<';
	public static char START_CODON = '>';
	public static char Methionie = 'M';// also can be start codon
	public static int NUMBER_OF_PROTEIN = 64;
	// table contains the
	public static char[] proteinTable = new char[] { 'K', 'N', 'K', 'N', 'T',
			'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H',
			'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L',
			'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G',
			'V', 'V', 'V', 'V', '<', 'Y', '<', 'Y', 'S', 'S', 'S', 'S', '<',
			'C', 'W', 'C', 'L', 'F', 'L', 'F' };

	/**
	 * Method converts array of characters into protein codons
	 * 
	 * @param seqData
	 *            CharSequenceData - sequence of A, C, G,T
	 * @return SequenceData - sequence of protein codons
	 * 
	 */
	public static char[] convertToProtein(char[] seqData)
			throws RuntimeException {

		char[] charSeq = seqData;
		// now run through the sequence, convert 3 contiguous characters into 1
		// protein
		// the number of proteins represented by that sequence will be the
		// length of the sequence - 2
		char[] proteinSeq = new char[charSeq.length - 2];// may it cause
															// exception when
															// charSeq < 2 ?
		for (int i = 0; i < charSeq.length - 2; i++) {
			proteinSeq[i] = getCodon(charSeq[i], charSeq[i + 1], charSeq[i + 2]);
		}
		return proteinSeq;
	}

	/**
	 * method returns the corrensponding codon characters to 3 characters
	 * 
	 * @param char , char, char
	 * @return char
	 */
	public static char getCodon(char c1, char c2, char c3) {
		return proteinTable[4 * 4 * getOrder(c1) + 4 * getOrder(c2)
				+ getOrder(c3)];
	}

	/**
	 * method finds the value of the character
	 * 
	 * @param char c if c='a' or 'A' --> return 0 if c='c' or 'G' --> return 1
	 *        if c='g' or 'G' --> reuturn 2 if c= 't' or 'T'--> return 3
	 * @return int
	 */
	public static int getOrder(char c) {
		char character = Character.toLowerCase(c);
		int returnValue;
		switch (character) {
		case 'a':
			returnValue = A;
			break;
		case 'c':
			returnValue = C;
			break;
		case 'g':
			returnValue = G;
			break;
		case 't':
			returnValue = T;
			break;
		default:
			returnValue = -1;
			break;
		}
		return returnValue;
	}

	/**
	 * this method filters the aminoAcidCodons to aminoAcids in protein
	 * 
	 * @param char[] aminoAcidCodons
	 * @return char[] aminoAcids in protein
	 * 
	 * 
	 */
	public static void proteinFilter(char[] aminoAcidCodons) {
		// int lastStop =0;
		// int lastBegin =0;
		boolean isInProtein = false;
		for (int i = 0; i < aminoAcidCodons.length; i++) {
			// if current char is stop and it is in the protien
			if (aminoAcidCodons[i] == STOP_CONDON && isInProtein == true) {
				aminoAcidCodons[i] = STOP_CONDON;
				isInProtein = false;
			}
			// if current char is start and not in the protein
			else if (aminoAcidCodons[i] == Methionie && isInProtein == false) {
				aminoAcidCodons[i] = START_CODON;
				isInProtein = true;
			} else if (isInProtein == false) {
				aminoAcidCodons[i] = ' ';
			}

		}

	}

}
