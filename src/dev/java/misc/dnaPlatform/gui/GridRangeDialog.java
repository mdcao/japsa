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

//This class is written by Julie Bernal and subsequently modified and maintained
//by Minh Duc Cao

package misc.dnaPlatform.gui;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;

/**
 * <p>
 * Title: GridRangeDialog
 * </p>
 * 
 * <p>
 * Description: This is a dialog used to enter minimum and maximum values for
 * both x and y axis
 * </p>
 * 
 * 
 * @author Julie Bernal
 * @version 1.0
 */
public class GridRangeDialog extends JDialog {
	public static final long serialVersionUID = MainFrame.serialVersionUID;

	JPanel contentPane;
	JPanel titlePane = new JPanel();
	JPanel coordPane = new JPanel();
	JPanel buttonPane = new JPanel();

	JLabel titleLabel = new JLabel("  Plot Title: ");
	JTextField titleField = new JTextField();

	JLabel yMinLabel = new JLabel("Min y: ");
	JTextField yMinField = new DoubleField();
	JLabel yMaxLabel = new JLabel("    Max y: ");
	JTextField yMaxField = new DoubleField();

	JLabel xMinLabel = new JLabel("Min x: ");
	JTextField xMinField = new DoubleField();
	JLabel xMaxLabel = new JLabel("    Max x: ");
	JTextField xMaxField = new DoubleField();

	JButton okButton = new JButton("OK");
	JButton cancelButton = new JButton("Cancel");

	// to hold new grid range entered by user
	double yMin;
	double yMax;
	double xMin;
	double xMax;

	// to hold title entered by user
	String title;

	/**
	 * Constructor takes as parameters current range values for a plot
	 * 
	 * @param xMinVal
	 *            double
	 * @param xMaxVal
	 *            double
	 * @param yMinVal
	 *            double
	 * @param yMaxVal
	 *            double
	 */
	public GridRangeDialog(String title, double xMinVal, double xMaxVal,
			double yMinVal, double yMaxVal) {

		this.title = title;
		xMin = xMinVal;
		xMax = xMaxVal;
		yMin = yMinVal;
		yMax = yMaxVal;

		enableEvents(AWTEvent.WINDOW_EVENT_MASK);
		try {
			jbInit();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		pack();
	}

	/**
	 * Creating the GUI of GridRangeDialog
	 * 
	 * @throws Exception
	 */
	private void jbInit() throws Exception {
		this.setTitle("Grid Range");
		this.setResizable(false);

		/* creating coordPane */
		yMinField.setText(yMin + "");
		yMaxField.setText(yMax + "");
		xMinField.setText(xMin + "");
		xMaxField.setText(xMax + "");

		coordPane.setLayout(new GridLayout(2, 4));
		coordPane.setBorder(BorderFactory
				.createTitledBorder("Enter Grid Range"));
		titleLabel.setHorizontalAlignment(SwingConstants.RIGHT);
		coordPane.add(xMinLabel);
		coordPane.add(xMinField);
		coordPane.add(xMaxLabel);
		coordPane.add(xMaxField);
		coordPane.add(yMinLabel);
		coordPane.add(yMinField);
		coordPane.add(yMaxLabel);
		coordPane.add(yMaxField);

		/* creating buttonPane */
		buttonPane.add(okButton);
		buttonPane.add(cancelButton);

		cancelButton
				.addMouseListener(new GridRangeDialog_cancelButton_mouseAdapter(
						this));
		okButton.addMouseListener(new GridRangeDialog_okButton_mouseAdapter(
				this));

		/* creating titlePane */
		titleField.setText(title);
		titlePane.setLayout(new BorderLayout());
		titlePane.add(titleLabel, BorderLayout.WEST);
		titlePane.add(titleField, BorderLayout.CENTER);

		/* adding all panels to contentPane */
		contentPane = (JPanel) this.getContentPane();
		contentPane.setLayout(new BorderLayout(2, 10));

		contentPane.add(titlePane, BorderLayout.NORTH);
		contentPane.add(coordPane, BorderLayout.CENTER);
		contentPane.add(buttonPane, BorderLayout.SOUTH);

	}

	/**
	 * Returns title of plot
	 * 
	 * @return String
	 */
	public String getTitle() {
		return title;
	}

	/**
	 * Returns minimum y value entered
	 * 
	 * @return double
	 */
	public double getYminValue() {
		return yMin;
	}

	/**
	 * Returns maximum y value entered
	 * 
	 * @return double
	 */
	public double getYmaxValue() {
		return yMax;
	}

	/**
	 * Returns minimum x value entered
	 * 
	 * @return double
	 */
	public double getXminValue() {
		return xMin;
	}

	/**
	 * Returns maximum x value entered
	 * 
	 * @return double
	 */
	public double getXmaxValue() {
		return xMax;
	}

	/**
	 * Gets minumum and maximum values for x and y axis
	 * 
	 * @param e
	 *            MouseEvent
	 */
	public void okButton_mouseClicked(MouseEvent e) {
		xMin = ((DoubleField) xMinField).getDouble();
		xMax = ((DoubleField) xMaxField).getDouble();
		yMin = ((DoubleField) yMinField).getDouble();
		yMax = ((DoubleField) yMaxField).getDouble();

		title = titleField.getText();
	}

	/**
	 * DoubleField is a JTextField for doubles
	 * */
	private class DoubleField extends JTextField {
		// public static final long serialVersionUID =
		// MainFrame.serialVersionUID;

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		DoubleField() {
			this.addKeyListener(new KeyAdapter() {
				public void keyTyped(KeyEvent e) {
					char c = e.getKeyChar();
					if (!(Character.isDigit(c) || c == KeyEvent.VK_BACK_SPACE
							|| c == KeyEvent.VK_DELETE
							|| (c == '.' && getText().indexOf('.') == -1) || (c == '-' && getText()
							.length() == 0))) {
						getToolkit().beep();
						e.consume();
					}
				}
			});
		}

		/**
		 * Returns double in text box, if the DoubleField is empty then 0 is
		 * returned
		 * 
		 * @return double
		 */
		double getDouble() {
			if (super.getText().length() < 0)
				return 0;

			return Double.parseDouble(super.getText());
		}
	}

}

class GridRangeDialog_cancelButton_mouseAdapter extends MouseAdapter {
	private GridRangeDialog adaptee;

	GridRangeDialog_cancelButton_mouseAdapter(GridRangeDialog adaptee) {
		this.adaptee = adaptee;
	}

	public void mouseClicked(MouseEvent e) {
		adaptee.dispose();
	}
}

class GridRangeDialog_okButton_mouseAdapter extends MouseAdapter {
	private GridRangeDialog adaptee;

	GridRangeDialog_okButton_mouseAdapter(GridRangeDialog adaptee) {
		this.adaptee = adaptee;
	}

	public void mouseClicked(MouseEvent e) {
		adaptee.okButton_mouseClicked(e);
		adaptee.dispose();
	}
}
