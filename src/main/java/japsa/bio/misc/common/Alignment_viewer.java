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

package japsa.bio.misc.common;

import java.io.*;

import japsa.bio.misc.dnaPlatform.gui.MainFrame;

import java.awt.*;
import java.awt.event.*;

class Alignment_viewer extends Panel {

	static final long serialVersionUID = MainFrame.serialVersionUID;
	Label l1, l2;
	ScrollPane sp;
	Panel p;

	public Alignment_viewer(String a1) {
		int p_ = a1.indexOf("\n");
		init(a1.substring(0, p_), a1.substring(p_ + 1));
	}

	public Alignment_viewer(String a1, String a2) {
		init(a1, a2);
	}

	private void init(String a1, String a2) {
		StringBuffer mat = new StringBuffer();
		for (int i = 0; i < a1.length() && i < a2.length(); i++) {
			char c1 = a1.charAt(i);
			char c2 = a2.charAt(i);
			mat.append((!Character.isLowerCase(c1) && Character.isLetter(c1) && c1 == c2) ? "|"
					: " ");
		}

		setFont(new Font("Fixed", Font.PLAIN, 12));
		setLayout(new GridLayout(1, 1));

		p = new Panel();
		p.setLayout(new GridLayout(3, 1));

		l1 = new Label(a1);
		l2 = new Label(a2);
		p.add(l1);
		p.add(new Label(mat.toString()));
		p.add(l2);

		sp = new ScrollPane();
		sp.add(p);
		add(sp);
	}

	public void init_size() {
		sp.setSize((int) p.getPreferredSize().getWidth(), (int) p
				.getPreferredSize().getHeight() + 80);
		setSize(sp.getPreferredSize());
	}

	public static void main(String args[]) {
		Alignment_viewer t = null;

		if (args.length == 0) {
			// No args, read sequence from stdin.
			String a1 = null, a2 = null;
			try {
				BufferedReader in = new BufferedReader(new InputStreamReader(
						System.in));
				a1 = in.readLine();
				a2 = in.readLine();
			} catch (Exception e) {
				System.err.println("Error stdin: " + e);
			}
			t = new Alignment_viewer(a1, a2);
		} else if (args.length == 1) {
			// One arg, assume it is a filename to read the sequences from
			String a1 = null, a2 = null;
			try {
				BufferedReader in = new BufferedReader(new FileReader(args[0]));
				a1 = in.readLine();
				a2 = in.readLine();
				in.close();
			} catch (Exception e) {
				System.err.println("Error reading '" + args[0] + "': " + e);
			}
			t = new Alignment_viewer(a1, a2);
		} else if (args.length == 2)
			t = new Alignment_viewer(args[0], args[1]);
		else {
			System.err.println("Usage: java Alignment_viewer <alignment>");
			System.exit(-1);
		}

		Frame f = new Frame("Sequence Alignment");
		f.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});

		f.add("Center", t);
		// f.setSize((int)t.getPreferredSize().getWidth(),
		// (int)t.getMinimumSize().getHeight());
		// f.show();
		f.setVisible(true);
		t.init_size();
		f.setSize(t.getPreferredSize().getWidth() < 800 ? (int) t
				.getPreferredSize().getWidth() : 800, (int) t
				.getPreferredSize().getHeight());
	}
}
