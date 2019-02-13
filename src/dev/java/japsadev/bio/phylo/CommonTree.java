package japsadev.bio.phylo;

import java.io.File;
import java.io.IOException;

import pal.tree.Tree;

public interface CommonTree {

	/**Returns String[][] res where res[0] is an array of taxonomy and res[1] is an array of css colors 
	  * includes current node all the way up to the root
	  * */
	String[][] getTaxonomy(String in);

	void print(File out) throws IOException;

	Tree[] getTrees();

}