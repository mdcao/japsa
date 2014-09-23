/*****************************************************************************
 * Copyright (c) Minh Duc Cao, Monash Uni & UQ, All rights reserved.         *
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  * 
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 * 3. Neither the names of the institutions nor the names of the contributors*
 *    may be used to endorse or promote products derived from this software  *
 *    without specific prior written permission.                             *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS   *
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, *
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR    *
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR         *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 ****************************************************************************/

/*                           Revision History                                
 * 04/01/2012 - Minh Duc Cao:                                         
 *  
 ****************************************************************************/

package japsa.xm;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import japsa.xm.expert.*;

import java.io.File;
import java.util.Random;

/**
 * @author Minh Duc Cao
 * 
 */

@Deployable(scriptName = "jsa.xm.compress", scriptDesc = "Compression of DNA/protein sequences")
public class ExpertModelDriver {


	public static void main(String[] args) throws Exception {
		Deployable annotation = ExpertModelDriver.class.getAnnotation(Deployable.class);
		CommandLine cmdLine = new CommandLine("\nUsage: "
				+ annotation.scriptName() + " [options] file1 file2 ...",
				annotation.scriptDesc());
		
			cmdLine.addInt("hashSize", 11, "Hash size");
			cmdLine.addInt("context", 15, "Length of the context");
			cmdLine.addInt("limit", 200, "Expert Limit");
			cmdLine.addDouble("threshold", 0.15, "Listen threshold");
			cmdLine.addInt("chance", 20, "Chances");
			cmdLine.addBoolean("binaryHash", false, "Use binary hash or not");
			cmdLine.addString("offsetType", "counts",
							"Way of update offset/palindrome expert: possible value count, subs");

			cmdLine.addString("real", null, "File name of the real compression");
			cmdLine.addString("decode", null, "File name of the encoded");
			cmdLine.addString("output", "decoded",
					"The output file of decoded file");
			cmdLine.addString("info", null, "File name of the infomation content");
			cmdLine.addString("markov", null,
					"File name of the markov infomation content");
			cmdLine.addBoolean("optimise", false,
							"Running in optimise mode, just report the entropy,recommended for long sequence");

			cmdLine.addInt("checkPoint", 1000000, "Frequency of check point");
			cmdLine.addString("hashType", "hash",
							"Type of Hash table: hash=hashtable, sft=SuffixTree,sfa = SuffixArray");
			cmdLine.addBoolean("selfRep", true,
					"Propose experts from the sequence to compressed?");
			

			//args = cmdLine.parseLine(args);
			args = cmdLine.stdParseLine(args);
			
			
			ExpertModel eModel = getExpertModel(cmdLine);

			System.out.println(ExpertModel.version());

			if (cmdLine.getStringVal("decode") == null
					&& (args == null || args.length <= 0)) {
				System.err
						.println(cmdLine.usage() + "\n");
				System.exit(1);
			}

			// Print out all params
			eModel.printParams();

			// cmdLine.printOptions();

			// #TIME_BEGIN
			long timeStart = System.currentTimeMillis();
			System.out.println(" #Reading file(s)");
			// #TIME_END

			Sequence[] dnaArray;
			if (cmdLine.getStringVal("decode") != null) {
				dnaArray = new Sequence[args.length + 1];
			} else {// Read per normal
				dnaArray = new Sequence[args.length];
			}
			Random rnd = new Random();
			for (int i = 0; i < args.length; i++) {
				// System.out.print("Read " + args[i] + "...");
				dnaArray[i] = SequenceReader.getReader(args[i]).nextSequence(null);// IOTools.read(args[i]);
				if(dnaArray[i].alphabet() != Alphabet.DNA4()){
					Sequence seq = new Sequence(Alphabet.DNA4(), dnaArray[i].toBytes());
					for (int xx = 0; xx < seq.length(); xx++){
						if (seq.getBase(xx) >= 4)
							seq.setBase(xx, (byte) rnd.nextInt(4) );
					}
					dnaArray[i] = seq;
				}
				// System.out.println(" done");
			}
			
			// #TIME_BEGIN
			long timeEnd = System.currentTimeMillis();
			System.out.println(" #Read file(s) in " + (timeEnd - timeStart)
					+ "ms ");
			// #TIME_END

			/*********************************************
			 * if (cmdLine.getStringVal("hashType").equals("sft")){ //Augment
			 * all the background sequences BioCompDNA[] augdnaArray = new
			 * BioCompDNA[2]; augdnaArray[1] = dnaArray[dnaArray.length - 1];
			 * augdnaArray[0] = new BioCompDNA(new byte[0],"Combine");
			 * 
			 * for (int x = 0; x < dnaArray.length - 1; x++){
			 * augdnaArray[0].concatenate(dnaArray[x]);
			 * augdnaArray[0].concatenate(dnaArray[x].reverseComplement()); }
			 * dnaArray = augdnaArray; }else if
			 * (cmdLine.getStringVal("hashType").equals("sfa") &&
			 * dnaArray.length > 1){ //Augment all the background sequences for
			 * suffix array BioCompDNA[] augdnaArray = new BioCompDNA[2];
			 * augdnaArray[1] = dnaArray[dnaArray.length - 1];
			 * 
			 * int combileLength = 0; for (int x = 0; x < dnaArray.length - 1;
			 * x++){ combileLength += dnaArray[x].length(); }
			 * 
			 * //Combine all sequences, including the rev comp into 1/ //in the
			 * reverse order. the correct order is retrieved later byte []
			 * combineByte = new byte[combileLength * 2];///Backward and forward
			 * int ind = 0;
			 * 
			 * for (int x = 0; x < dnaArray.length - 1; x++){ byte[] seqX =
			 * dnaArray[x].toBytes(); for (int y = 0; y < seqX.length; y++){ //
			 * combineByte[ind ++] = seqX[y];//(byte)(3 - seqX[y]);//compliment
			 * seqX[y]; combineByte[combineByte.length - ind] = (byte)(3 -
			 * seqX[y]);//seqX[y];//(byte)(3 - seqX[y]);//compliment } }
			 * 
			 * augdnaArray[0] = new BioCompDNA(combineByte,"Combine"); dnaArray
			 * = augdnaArray; } /
			 *********************************************/

			for (int i = 0; i < dnaArray.length; i++) {
				if (i == dnaArray.length - 1) {
					if (cmdLine.getStringVal("decode") == null) {
						System.out.println("Encode  : " + dnaArray[i].getName() + "\t"
								+ dnaArray[i].length() + "");
					}
				} else
					System.out.println("Context : " + dnaArray[i].getName() + "\t"
							+ dnaArray[i].length() + "");
			}
			System.out
					.println("----------------------------------------------------------------------");

			long start;
			/*************************************************************************/
			if (cmdLine.getBooleanVal("optimise")) {
				System.gc();
				start = System.currentTimeMillis();

				double cost = eModel.encode_optimise(dnaArray);
				long time = (System.currentTimeMillis() - start); //
				System.out.printf("%f bps in %d ms\n", cost, time);
				System.out
						.println("=============================================================================");

			} else if (cmdLine.getStringVal("decode") != null) {
				String file = cmdLine.getStringVal("decode");
				String output = cmdLine.getStringVal("output");

				// dnaArray[dnaArray.length - 1] = new BioCompDNA();

				System.out.println("Decoding");
				start = System.currentTimeMillis();

				eModel.decode(dnaArray, new File(file));

				System.out.println(" Time decode "
						+ (System.currentTimeMillis() - start) / 1000.0
						+ "seconds\n");
				dnaArray[dnaArray.length - 1].writeFasta(output);

			} else if (cmdLine.getStringVal("real") != null) {
				System.gc();
				System.out.println("Real  encoding");
				start = System.currentTimeMillis();
				File outputFile = eModel.realEncode(dnaArray, cmdLine
						.getStringVal("real"));

				System.out.println(" Time encode "
						+ (System.currentTimeMillis() - start) / 1000.0
						+ "seconds");
				System.out.println(" Encoding cost = "
						+ (outputFile.length() * 8.0)
						/ dnaArray[dnaArray.length - 1].length() + "bps");

			} else if (cmdLine.getStringVal("info") != null) {
				System.gc();
				start = System.currentTimeMillis();
				eModel.encode(dnaArray, cmdLine.getStringVal("info"), cmdLine
						.getStringVal("markov"));// ,args[1]);
				System.out.println("Total time "
						+ (System.currentTimeMillis() - start) + "ms");
			} else {// Normal
				/*************************************************************************/
				System.gc();
				start = System.currentTimeMillis();
				// seqHash[1] = japsa.seq;
				double costs = eModel.encode1(dnaArray);// ,args[1]);

				long time = (System.currentTimeMillis() - start);
				// System.out.printf(" Comp rate : #%f#\n",total /
				// costs.length);
				System.out.printf("%f bps in %d ms\n", costs, time);
				System.out
						.println("=============================================================================");
			}
			/*************************************************************************/
		
	}
	
	
	
	public static ExpertModel getExpertModel(CommandLine cmdLine)
			throws Exception {
		
		
		

		// Params for expert model		
		int hashSize = cmdLine.getIntVal("hashSize");
		int context = cmdLine.getIntVal("context");
		int expertsLimit = cmdLine.getIntVal("limit");
		double listenThreshold = cmdLine.getDoubleVal("threshold");
		int chances = cmdLine.getIntVal("chance");

		ExpertModel eModel = new ExpertModel(hashSize, Alphabet.DNA4(), context,
				expertsLimit, listenThreshold, chances, cmdLine
						.getBooleanVal("binaryHash"));
		eModel.setHashType(cmdLine.getStringVal("hashType"));

		eModel.setSelfRep(cmdLine.getBooleanVal("selfRep"));

		if ("subs".equals(cmdLine.getStringVal("offsetType"))) {
			eModel.offSetSeed = new RepeatSubsExpert(new Sequence(null,0), 0, null,
					RepeatExpert.COPY_TYPE);
			eModel.palinSeed = new RepeatSubsExpert(new Sequence(null,0), 0, null,
					RepeatExpert.PALIN_TYPE);
			/**************************************************************/
			// }else if
			// ("viaprotein".equals(cmdLine.getStringVal("offsetType"))){
			// eModel.offSetSeed = new OffsetViaProteinExpert(new
			// byte[0],0,null);
			// eModel.palinSeed = new PalindromeViaProteinExpert(new
			// byte[0],0,null);
			// //}else if ("static".equals(cmdLine.getStringVal("offsetType"))){
			// // eModel.offSetSeed = new OffsetStaticExpert(new
			// byte[0],0,null);
			// // eModel.palinSeed = new PalindromeCountExpert(new
			// byte[0],0,null);
		} else {
			eModel.offSetSeed = new RepeatCountExpert(new Sequence(null,0), 0, null,
					RepeatExpert.COPY_TYPE);
			eModel.palinSeed = new RepeatCountExpert(new Sequence(null,0), 0, null,
					RepeatExpert.PALIN_TYPE);
		}
		/**************************************************************/
		// #CHECKPOINT_BEGIN
		eModel.checkPoint = cmdLine.getIntVal("checkPoint");
		// #CHECKPOINT_END

		return eModel;
	}
	

}
