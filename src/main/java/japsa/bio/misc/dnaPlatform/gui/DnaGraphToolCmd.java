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

package japsa.bio.misc.dnaPlatform.gui;

import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import javax.swing.*;

/**
 * <p>
 * Title: DNA Graph Tool
 * </p>
 * 
 * <p>
 * Description: This is the main class of the DNAPlatform
 * </p>
 * 
 * <p>
 * Copyright: Copyright (c) 2005
 * </p>
 * 
 * <p>
 * Company: Monash
 * </p>
 * 
 * @author Julie Bernal
 * @version 1.0
 */

@Deployable(scriptName = "jsa.dnaGraph", scriptDesc = "Visualisation")
public class DnaGraphToolCmd extends CommandLine{	
	public DnaGraphToolCmd() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addStdHelp();		
	}

	public void show(){
		@SuppressWarnings("unused")
		MainFrame frame = new MainFrame();

	}


	public static void main(String[] args){ 
		DnaGraphToolCmd cmdLine = new DnaGraphToolCmd();
		cmdLine.stdParseLine(args);

		/**********************************************************************/
		System.setProperty("java.awt.headless", "false");

		try {
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		}catch (Exception e) {
			e.printStackTrace();
		}
		cmdLine.show();
	}

}
