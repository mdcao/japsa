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
import java.awt.*;
import javax.swing.tree.*;
import java.util.*;

/**
 * <p>
 * Title: DataTreePanel
 * </p>
 * 
 * <p>
 * Description: This is a Panel that displays a JTree given a group of objects
 * inside a Vector. toString method is called for every object in the tree and
 * if this method returns a string with new line characters each line is added
 * as a new node in the tree. If the string has curly braces, all lines between
 * the curly braces are added as child nodes
 * </p>
 * 
 * @author Julie Bernal
 * @version 1.0
 */
@SuppressWarnings("rawtypes")
public class DataTreePanel extends JPanel {
	public static final long serialVersionUID = MainFrame.serialVersionUID;

	private JTree tree;
	JScrollPane treeView;

	private DefaultMutableTreeNode rootNode;
	private DefaultTreeModel treeModel;

	public DataTreePanel() {

		try {
			jbInit();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	private void jbInit() throws Exception {
		setLayout(new BorderLayout());
		setBorder(BorderFactory.createTitledBorder("Data"));

		createTree();
		treeView = new JScrollPane(tree);
		this.add(treeView, BorderLayout.CENTER);
	}

	private void createTree() {
		rootNode = new DefaultMutableTreeNode("Data");
		treeModel = new DefaultTreeModel(rootNode);

		tree = new JTree(treeModel);
		tree.getSelectionModel().setSelectionMode(
				TreeSelectionModel.SINGLE_TREE_SELECTION);
		tree.setRootVisible(false);
		tree.setShowsRootHandles(true);
		tree.putClientProperty("JTree.lineStyle", "Angled");

		DefaultTreeCellRenderer renderer2 = new DefaultTreeCellRenderer();
		renderer2.setOpenIcon(null);
		renderer2.setClosedIcon(null);
		renderer2.setLeafIcon(null);
		tree.setCellRenderer(renderer2);

	}

	/**
	 * Returns the JTree displayed in this panel
	 * 
	 * @return JTree
	 */
	public JTree getTree() {
		return tree;
	}

	private Vector<String> createVectorString(String description) {
		Vector<String> v = new Vector<String>();
		String[] children = description.split("[\\n]");

		for (int i = 0; i < children.length; i++) {
			if (!children[i].matches("^\\s*$"))
				v.add(children[i]);
		}
		return v;
	}

	/**
	 * Given a vector of vectors and strings describing the history of a
	 * SequenceData object it adds this history to the rootNode of tree to be
	 * displayed
	 * 
	 * @param history
	 *            Vector
	 */
	public void addDataHistory(Vector history) {
		DefaultMutableTreeNode parent = rootNode;
		DefaultMutableTreeNode current = parent;

		Iterator it = history.iterator();
		while (it.hasNext()) {
			Object o = it.next();
			if (o instanceof Vector) {
				// if element is a vector then call function to add elements
				// of vector as children to current node
				addDataHistory((Vector) o, current);
			} else {
				String description = o.toString();

				// if toString() returns a string with curly braces
				// add all lines in between as child nodes
				if ((description.indexOf('{') > -1)
						&& (description.indexOf('}') > -1)) {

					// add anything before brakets into tree
					String desStart = description.substring(0,
							description.indexOf('{'));
					if (desStart.matches("\\S+"))
						current = addObject(parent, desStart, true);

					// add children to a vector and call function to add
					// children to tree
					String des = description.substring(
							(description.indexOf('{') + 1),
							(description.lastIndexOf('}')));

					addDataHistory(createVectorString(des), current);

					// add anything after children
					String end = description.substring(
							description.lastIndexOf('}'),
							description.length() - 1);
					if (end.matches("\\S+"))
						current = addObject(parent, end, true);
				} else
					current = addObject(parent, o.toString(), true);

			}

		}
	}

	/**
	 * Given a vector of vectors and strings describing the history of a
	 * SequenceData object it adds this history to a parent
	 * DefaultMutableTreeNode to be displayed.
	 * 
	 * @param history
	 *            Vector
	 * @param parent
	 *            DefaultMutableTreeNode
	 */
	public void addDataHistory(Vector history, DefaultMutableTreeNode parent) {
		DefaultMutableTreeNode current = parent;
		Iterator it = history.iterator();
		while (it.hasNext()) {
			Object o = it.next();
			if (o instanceof Vector) {
				// if element is a vector then call function to add elements
				// of vector as children to current node
				addDataHistory((Vector) o, current);
			} else {
				String description = o.toString();

				// if toString() returns a string with curly braces
				// add all lines in between as child nodes
				if ((description.indexOf('{') > -1)
						&& (description.indexOf('}') > -1)) {

					// add anything before brakets into tree
					String desStart = description.substring(0,
							description.indexOf('{'));
					if (desStart.matches("\\S+"))
						current = addObject(parent, desStart, false);

					// add children to a vector and call function to add
					// children to tree
					String des = description.substring(
							(description.indexOf('{') + 1),
							(description.lastIndexOf('}')));

					addDataHistory(createVectorString(des), current);

					// add anything after children
					String end = description.substring(
							description.lastIndexOf('}'),
							description.length() - 1);
					if (end.matches("\\S+"))
						current = addObject(parent, end, false);
				} else
					current = addObject(parent, o.toString(), false);
			}
		}
	}

	/**
	 * Adds child object to given DefaultMutableTreeNode parent.
	 * 
	 * @param parent
	 *            DefaultMutableTreeNode
	 * @param child
	 *            Object
	 * @param shouldBeVisible
	 *            boolean
	 * @return DefaultMutableTreeNode
	 */
	public DefaultMutableTreeNode addObject(DefaultMutableTreeNode parent,
			Object child, boolean shouldBeVisible) {
		DefaultMutableTreeNode childNode = new DefaultMutableTreeNode(child);

		treeModel.insertNodeInto(childNode, parent, parent.getChildCount());

		// Make sure the user can see the lovely new node.
		if (shouldBeVisible) {
			tree.scrollPathToVisible(new TreePath(childNode.getPath()));
		}

		return childNode;
	}

	/**
	 * Returns selected sequence data name if a path has been selected.
	 * 
	 * @return String
	 */
	public String getFirstSelectedPath() {
		if (!tree.isSelectionEmpty()) {
			Object[] path = tree.getSelectionPath().getPath();
			return path[1].toString();
		} else
			return null;
	}

	/**
	 * Removes selected node
	 * 
	 * @return boolean
	 */
	public boolean removeSelectedPath() {
		TreePath currentSelection = tree.getSelectionPath();
		if (currentSelection != null) {
			DefaultMutableTreeNode currentNode = (DefaultMutableTreeNode) (currentSelection
					.getPathComponent(1));
			MutableTreeNode parent = (MutableTreeNode) (currentNode.getParent());
			if (parent != null) {
				treeModel.removeNodeFromParent(currentNode);
				return true;
			}
		}
		return false;
	}

	/*
	 * // Adding child node to selected path in tree. // Might need this later.
	 * 
	 * public DefaultMutableTreeNode addObject(Object child) {
	 * DefaultMutableTreeNode parentNode = null; TreePath parentPath =
	 * tree.getSelectionPath();
	 * 
	 * if (parentPath == null) { //There's no selection. Default to the root
	 * node. parentNode = rootNode; } else { parentNode =
	 * (DefaultMutableTreeNode) (parentPath.getLastPathComponent()); }
	 * 
	 * return addObject(parentNode, child, true); }
	 */

}
