package clustering;

import java.io.*;
import java.util.*;



/**
 * @author buvan.suji
 *
 */

public class GettingTreadsFromFasta {

	private String[] description;
	private String[] sequence;

	static ArrayList<String> Rname = new ArrayList<String>();
	static ArrayList<String> TReads = new ArrayList<String>();
	static String FnamePath;
	static int sequenceLength;

	public GettingTreadsFromFasta(String filename) {
		readSequenceFromFile(filename);
	}

	@SuppressWarnings("resource")
	void readSequenceFromFile(String file) {

		ArrayList<String> desc = new ArrayList<String>();
		ArrayList<String> seq = new ArrayList<String>();
		BufferedReader in = null;

		try {
			in = new BufferedReader(new FileReader(file));
			StringBuffer buffer = new StringBuffer();
			String line = in.readLine();

			if (line == null)
				throw new IOException(file + " is an empty file");

			if (line.charAt(0) != '>')
				throw new IOException("First line of " + file
						+ " should start with '>'");
			else
				desc.add(line);

			for (line = in.readLine().trim(); line != null; line = in
					.readLine()) {
				if (line.length() > 0 && line.charAt(0) == '>') {
					seq.add(buffer.toString());
					buffer = new StringBuffer();
					desc.add(line);
				} else
					buffer.append(line.trim());
			}
			if (buffer.length() != 0)
				seq.add(buffer.toString());
		} catch (IOException e) {
			System.out.println("Error when reading " + file);
			e.printStackTrace();
		}

		description = new String[desc.size()];
		sequence = new String[seq.size()];
		for (int i = 0; i < seq.size(); i++) {
			description[i] = (String) desc.get(i);
			sequence[i] = (String) seq.get(i);
		}
		sequenceLength = seq.size();

	}

	// return first sequence as a String
	public String getSequence() {
		return sequence[0];
	}

	// return first xdescription as String
	public String getDescription() {
		return description[0];
	}

	// return sequence as a String
	public String getSequence(int i) {
		return sequence[i];
	}

	// return description as String
	public String getDescription(int i) {
		return description[i];
	}

	public int size() {
		return sequence.length;
	}

	public static void FileReading(String filename) {
		GettingTreadsFromFasta fsf = new GettingTreadsFromFasta(
				filename);
		StringBuffer buffer1 = new StringBuffer();

		for (int i = 0; i < fsf.size(); i++) {
			String temp = ((buffer1.append(fsf.getDescription(i)))
					.deleteCharAt(0)).toString();
			Rname.add(temp);
			TReads.add(fsf.getSequence(i));
			buffer1 = new StringBuffer();
		}
	}

	// call this method-1
	public static void DestReads() throws Exception {
		String fname = "";
		System.out.print("Enter the name of the FastaFile:");
		fname = (new BufferedReader(new InputStreamReader(System.in)))
				.readLine();
		FnamePath = new File(fname).getName();
		// System.out.println(FnamePath);
		FileReading(fname);
	}

	public static void main(String[] args) throws Exception {
		DestReads();
		ViewList();
	}

	// call this method-2
	public static ArrayList<String> GetRname() {
		return Rname;
	}

	// call this method-3
	public static ArrayList<String> GetTReads() {
		return TReads;
	}

	public static int NumberReads() {
		return TReads.size();
	}

	// call this method-4
	public static String GetFileName() {
		return FnamePath;
	}

	public static int SeqLength() {
		return sequenceLength;
	}

	public static int NelementsClustering() {
		return sequenceLength * (sequenceLength - 1) / 2;
	}

	public static void ViewList() {
		System.out.println(GetRname());
		// System.out.println(GetTReads());
		System.out.println(SeqLength());
		System.out.println(NelementsClustering());
	}



}
