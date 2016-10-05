/*  
 *  Copyright (c) David Powell <david@drp.id.au>
 *
 * 
 * This file is used by both FuzzyLZ and AlignCompress

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


 */

package misc.common;

public interface Two_Seq_Model_Counts extends Two_Seq_Model {
	public abstract void update_count_encA(Counts c, double w, char a, int i);

	public abstract void update_count_encB(Counts c, double w, char a, int i);

	public abstract void update_count_encBoth(Counts c, double w, char a,
			char b, int i, int j);

	public abstract Params counts_to_params(Counts c);
}
