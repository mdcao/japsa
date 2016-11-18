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

package japsadev.misc.dnaPlatform.gui;

import javax.swing.*;

import japsadev.misc.dnaPlatform.OptionsHandle;
import japsadev.misc.dnaPlatform.sequence.*;

import java.awt.*;
import java.awt.event.KeyEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.util.*;

/**
 * <p>
 * Title: OptionsHandleDialog
 * </p>
 * 
 * <p>
 * Description: This is a dialog that displays options in a OptionsHandle object
 * and lets user set options
 * </p>
 * 
 * @author Julie Bernal
 * @version 1.0
 */
@SuppressWarnings("rawtypes")
class OptionsHandleDialog extends JDialog {
	/**
	 * 
	 */
	private static final long serialVersionUID = -6467708920193257046L;

	OptionsHandle myHandle;

	private JLabel[] options;
	private Object[] values;

	JPanel contentPane;
	JPanel northPane = new JPanel();
	JPanel southPane = new JPanel();
	JScrollPane centerPane = new JScrollPane();
	JPanel optionPane = new JPanel();

	JLabel title_lbl = new JLabel();
	JLabel subtitle_lbl = new JLabel("Options: ");
	JButton cancel_btn = new JButton("Cancel");
	JButton set_btn = new JButton("Run");

	Vector sequenceData;

	// used to indicate whether user clicked on ok or cancel
	private boolean accepted = false;

	public OptionsHandleDialog(Frame parent, OptionsHandle options, Vector data) {
		super(parent);
		enableEvents(AWTEvent.WINDOW_EVENT_MASK);
		myHandle = options;
		sequenceData = data;

		if (myHandle != null) {
			try {
				jbInit();
			} catch (Exception ex) {
				ex.printStackTrace();
			}
		}
		pack();
	}

	public OptionsHandleDialog(JDialog parent, OptionsHandle options,
			Vector data) {
		super(parent);
		enableEvents(AWTEvent.WINDOW_EVENT_MASK);
		myHandle = options;
		sequenceData = data;

		if (myHandle != null) {
			try {
				jbInit();
			} catch (Exception ex) {
				ex.printStackTrace();
			}
		}
		pack();
	}

	/*
	 * jbInit: in this function we create a ParamPanel depending on its
	 * Parameters, four types of parameters are recognized: boolean, integer,
	 * double and string
	 */
	private void jbInit() throws Exception {
		int numberParams = myHandle.getNumberOptions();

		this.setTitle("Set Options");
		this.setSize(500, 300);
		contentPane = (JPanel) this.getContentPane();
		contentPane.setLayout(new BorderLayout());

		title_lbl.setFont(new java.awt.Font("Dialog", Font.BOLD, 14));

		title_lbl.setText(myHandle.getOwner().toString());
		northPane.add(title_lbl);
		northPane.add(subtitle_lbl);

		cancel_btn
				.addActionListener(new OptionsHandleDialog_cancel_btn_actionAdapter(
						this));
		set_btn.addActionListener(new OptionsHandleDialog_set_btn_actionAdapter(
				this));
		southPane.add(cancel_btn);
		southPane.add(set_btn);

		// Creating the option panel:

		optionPane.setLayout(new GridLayout(numberParams, 2, 4, 4));
		createMyOptionsValues(numberParams);

		// adding options and values into this panel
		for (int i = 0; i < numberParams; i++) {
			String o = myHandle.getOptionAt(i);
			Object val = myHandle.getOptionValue(o);

			optionPane.add(options[i]);

			if (val instanceof Boolean)
				optionPane.add((JCheckBox) values[i]);

			else if (val instanceof Integer)
				optionPane.add((IntField) values[i]);

			else if (val instanceof Double)
				optionPane.add((DoubleField) values[i]);

			else if (val instanceof String)
				optionPane.add((JTextField) values[i]);

			else if (val instanceof SequenceData)
				optionPane.add((JComboBox) values[i]);

			else if (val instanceof OptionsHandle)
				optionPane.add((HandleButton) values[i]);

		}

		centerPane.getViewport().add(optionPane);
		contentPane.add(northPane, BorderLayout.NORTH);
		contentPane.add(southPane, BorderLayout.SOUTH);
		contentPane.add(centerPane, BorderLayout.CENTER);
	}

	@SuppressWarnings("unchecked")
	private void createMyOptionsValues(int numberParams) {
		options = new JLabel[numberParams];
		values = new Object[numberParams];

		// creating options[] and values[] from Options
		for (int i = 0; i < numberParams; i++) {
			String option = myHandle.getOptionAt(i);

			options[i] = new JLabel(option);
			Object val = myHandle.getOptionValue(option);

			if (val instanceof Boolean) {
				JCheckBox checkBox = new JCheckBox();
				checkBox.setSelected(((Boolean) val).booleanValue());
				values[i] = checkBox;
			}

			else if (val instanceof Integer) {
				IntField intField = new IntField();
				intField.setText(((Integer) val).intValue());
				values[i] = intField;
			}

			else if (val instanceof Double) {
				DoubleField doubleField = new DoubleField();
				doubleField.setText(((Double) val).doubleValue());
				values[i] = doubleField;

			}

			else if (val instanceof String) {
				values[i] = new JTextField((String) val);
			}

			else if (val instanceof SequenceData) {
				JComboBox seqBox = new JComboBox(sequenceData);
				values[i] = seqBox;
			}

			else if (val instanceof OptionsHandle) {
				OptionsHandle subHandle = (OptionsHandle) val;
				HandleButton handleButton = new HandleButton(
						option + "_Handle", subHandle);
				handleButton.addActionListener(new HandleButton_actionAdapter(
						this));
				values[i] = handleButton;
			}

		}
	}

	public void openSubHandleDialog(ActionEvent e) {

		// get option of button
		HandleButton button = (HandleButton) e.getSource();
		OptionsHandle handle = button.getOptionsHandle();

		OptionsHandleDialog dlg = new OptionsHandleDialog(this, handle,
				sequenceData);
		Dimension dlgSize = dlg.getPreferredSize();
		Dimension frmSize = getSize();
		Point loc = getLocation();
		dlg.setLocation((frmSize.width - dlgSize.width) / 2 + loc.x,
				(frmSize.height - dlgSize.height) / 2 + loc.y);
		dlg.setModal(true);

		dlg.setVisible(true);

	}

	/*
	 * setOptionsHandle: this function changes the option values in
	 * OptionsHandle
	 */
	public void setOptionsHandle() {
		int numberParams = myHandle.getNumberOptions();

		for (int i = 0; i < numberParams; i++) {
			String option = myHandle.getOptionAt(i);
			Object val = myHandle.getOptionValue(option);

			if (val instanceof Boolean) {
				JCheckBox checkBox = (JCheckBox) values[i];
				myHandle.setOptionValue(option,
						new Boolean(checkBox.isSelected()));
			}

			else if (val instanceof Integer) {
				IntField intField = (IntField) values[i];
				if (intField.getText().length() > 0)
					myHandle.setOptionValue(option,
							new Integer(intField.getText()));
			}

			else if (val instanceof Double) {
				DoubleField doubleField = (DoubleField) values[i];
				if (doubleField.getText().length() > 0)
					myHandle.setOptionValue(option,
							new Double(doubleField.getText()));

			}

			else if (val instanceof String) {
				JTextField txt = (JTextField) values[i];
				if (txt.getText().length() > 0)
					myHandle.setOptionValue(option, txt.getText());
			}

			else if (val instanceof SequenceData) {
				JComboBox seqBox = (JComboBox) values[i];
				SequenceData data = (SequenceData) seqBox.getSelectedItem();
				if (data instanceof SequenceData)
					myHandle.setSequenceDataValue(option, data);
			}

			else if (val instanceof OptionsHandle) {
				HandleButton handleButton = (HandleButton) values[i];
				OptionsHandle subhandle = handleButton.getOptionsHandle();
				if (subhandle instanceof OptionsHandle)
					myHandle.setOptionsHandleValue(option, subhandle);

			}
		}

		accepted = true;

	}

	/* IntField is a JTextField that only allows user to enter integers */
	private class IntField extends JTextField {
		/**
		 * 
		 */
		private static final long serialVersionUID = -4791170332146987317L;

		IntField() {
			this.addKeyListener(new KeyAdapter() {
				public void keyTyped(KeyEvent e) {
					char c = e.getKeyChar();
					if (!(Character.isDigit(c) || c == KeyEvent.VK_BACK_SPACE
							|| c == KeyEvent.VK_DELETE || (c == '-' && getText()
							.length() == 0))) {
						getToolkit().beep();
						e.consume();
					}
				}
			});
		}

		void setText(int number) {
			super.setText(String.valueOf(number));
		}

	}

	/* DoubleField is a JTextField that only allows user to enter doubles */
	private class DoubleField extends JTextField {
		/**
		 * 
		 */
		private static final long serialVersionUID = -7249595281088991083L;

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

		void setText(double number) {
			super.setText(String.valueOf(number));
		}
	}

	/* HandleButton is a JButton that contains a OptionsHandle object */
	private class HandleButton extends JButton {
		/**
		 * 
		 */
		private static final long serialVersionUID = 9207613471442424795L;
		private OptionsHandle myHandle;

		HandleButton(String text, OptionsHandle handle) {
			super(text);
			myHandle = handle;
		}

		public OptionsHandle getOptionsHandle() {
			return myHandle;
		}

	}

	/**
	 * This function returns a boolean indicating whether user clicked on button
	 * to set options
	 * 
	 * @return boolean
	 */
	public boolean isOptionsHandleSet() {
		return accepted;
	}

	protected void processWindowEvent(WindowEvent e) {
		if (e.getID() == WindowEvent.WINDOW_CLOSING) {
			dispose();
		}
		super.processWindowEvent(e);
	}

}

/*--------------*/
class HandleButton_actionAdapter implements ActionListener {
	private OptionsHandleDialog adaptee;

	HandleButton_actionAdapter(OptionsHandleDialog adaptee) {
		this.adaptee = adaptee;
	}

	public void actionPerformed(ActionEvent e) {
		adaptee.openSubHandleDialog(e);
	}
}

/*--------------*/
class OptionsHandleDialog_cancel_btn_actionAdapter implements ActionListener {
	private OptionsHandleDialog adaptee;

	OptionsHandleDialog_cancel_btn_actionAdapter(OptionsHandleDialog adaptee) {
		this.adaptee = adaptee;
	}

	public void actionPerformed(ActionEvent e) {
		adaptee.dispose();
	}
}

/*--------------*/
class OptionsHandleDialog_set_btn_actionAdapter implements ActionListener {
	private OptionsHandleDialog adaptee;

	OptionsHandleDialog_set_btn_actionAdapter(OptionsHandleDialog adaptee) {
		this.adaptee = adaptee;
	}

	public void actionPerformed(ActionEvent e) {
		adaptee.setOptionsHandle();
		adaptee.dispose();
	}
}
