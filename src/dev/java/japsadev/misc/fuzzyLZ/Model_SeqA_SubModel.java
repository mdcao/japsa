/*  
 *  Copyright (c) David Powell <david@drp.id.au>
 *
 * 
 * This file is part of FuzzyLZ 
 *
 * FuzzyLZ is a program orginally intended for the
 * compression of DNA sequeces.  It can be viewed as a
 * compression model like Lempel-Ziv 77, but instead of
 * exact matches, allowing matches that contain
 * inserts/deletes/mismatches.
 *
 */

package japsadev.misc.fuzzyLZ;

import japsadev.misc.common.MyMath;
import japsadev.misc.common.Params;
import japsadev.misc.common.Seq_Model;

/**
 * 
 * This is a simple extension to {@link Model_SeqA}. The difference is that this
 * model takes a parameter, a {@link japsadev.misc.common.Seq_Model}, that is used to
 * encode the characters from Sequence A, which <b>do not</b> match characters
 * from Sequence B. This model can be seen as a generalisation of
 * {@link Model_SeqA} which uses a uniform model for characters from unmatched
 * characters from sequence A.
 */
class Model_SeqA_SubModel extends Model_SeqA {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	Seq_Model subMdl;

	public Model_SeqA_SubModel(Seq_Model subMdl, Params p, int alphaSize,
			int countIndex) {
		super(p, alphaSize, countIndex);
		this.subMdl = subMdl;
	}

	public double encA(char a, int i) {
		return subMdl.encodeLen(a, i);
	}

	public double encB(char a, int i) {
		return 0;
	}

	public double encBoth(char a, char b, int i, int j) {
		if (a == b)
			return match_cost;

		// Encode char 'a', but normalise since we know 'a' is different to 'b'
		double costB = subMdl.encodeLen(b, i);
		double norm = MyMath.log2(1 - MyMath.exp2(-costB));
		return change_cost + subMdl.encodeLen(a, i) - norm;
	}

}
