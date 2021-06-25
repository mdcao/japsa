package japsadev.bio.hts.clustering;

import java.util.ArrayList;
import java.io.*;

import japsa.seq.SequenceOutputStream;
import japsadev.bio.hts.clustering.PairDistance;
import japsadev.bio.hts.clustering.GettingTreadsFromFasta;
import japsa.seq.*;

/**
 * @author buvan.suji
 *
 */

public class WriteClusterResultOnFile {	
	
	public static void writeOnFile(ArrayList<ArrayList<String>> 
	list, Sequence cons1, Sequence cons2, String FileName ) throws Exception{
	
	File file = new File("ClusterResult_"+FileName+".fasta");
	FileWriter fw = new FileWriter(file.getAbsoluteFile());
	BufferedWriter bw = new BufferedWriter(fw);	
				
	
	ArrayList<String> list1 = new ArrayList<String>();
	list1 = list.get(0);
	
	bw.write("Minimum Read Length: "+ list1.get(0));
	bw.newLine();
	bw.write("Maximum Read Length: "+ list1.get(1));
	bw.newLine();
	bw.write("Estimated Time: "+ list1.get(2));
	bw.newLine();
	bw.newLine();
	bw.write("The consensus sequence of cluster1: ");
	bw.newLine();
	bw.write(""+cons1);
	bw.newLine();
	bw.newLine();
	bw.write("The C1 members are: ");
	bw.newLine();
	for(int x=1;x<list.size();x++){
		String s = "C"+(x)+"{ ";
		bw.write(s);
		ArrayList<String> tempList = new ArrayList<String>();
		tempList = list.get(x);
				
		for(int y=0;y<tempList.size();y++){
			bw.write(tempList.get(y));
			bw.newLine();		
			
		}
		bw.write("}");
		bw.newLine();
		if(x==1){
			bw.newLine();
			bw.write("The consensus sequence of cluster2: ");
			bw.newLine();
			bw.write(""+cons2);
			bw.newLine();
			bw.newLine();
			bw.write("The C2 members are: ");
			bw.newLine();
		}
		
		
	}
		
	bw.close();
	System.out.println("Cluster job has been done");
}
}
