package japsadev.tools;

import java.io.IOException;

import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import japsadev.bio.hts.barcode.*;


@Deployable(
		scriptName = "jsa.dev.barcode", 
		scriptDesc = "Clustering nanopore sequences based on barcode"
		)
public class BarCodeAnalysisCmd extends CommandLine{
	public BarCodeAnalysisCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc()); 

		addString("bcFile", null, "Barcode file",true);		
		addString("seqFile", null, "Nanopore sequences file",true);
		addString("scriptRun", null, "Invoke command script to run npScarf",true);
		addDouble("threshold", .7, "Minimum identity for barcode alignment");
		addBoolean("twoends", false, "Whether a read must contain barcode sequence from both ends or just one end (default)");
		addBoolean("print", false, "Print out demultiplexed reads to corresponding FASTA file or not.");
		addStdHelp();
	}
	public static void main(String[] args) throws IOException, InterruptedException{
		CommandLine cmdLine = new BarCodeAnalysisCmd ();
		args = cmdLine.stdParseLine(args);

		String bcFile = cmdLine.getStringVal("bcFile");
		String script = cmdLine.getStringVal("scriptRun");
		String seqFile = cmdLine.getStringVal("seqFile");
		Double threshold = cmdLine.getDoubleVal("threshold");

		BarCodeAnalysis.print = cmdLine.getBooleanVal("print");
		BarCodeAnalysis.strict = cmdLine.getBooleanVal("twoends");

		
		BarCodeAnalysis bc = new BarCodeAnalysis(bcFile,script);
		bc.setThreshold(threshold);
		bc.clustering(seqFile);

		
	}
}
/*RST*
---------------------------------------------------------------------------
*barcode*: real-time de-multiplexing Nanopore reads from barcode sequencing
---------------------------------------------------------------------------

*barcode* (jsa.np.barcode) is a program that demultiplex the nanopore reads from 
Nanopore barcode sequencing. Downstream analysis can be invoked concurrently by an input script.

*barcode* is included in the `Japsa package <http://mdcao.github.io/japsa/>`_.

<usage>

~~~~~~~~~~~~~~
Usage examples
~~~~~~~~~~~~~~

A summary of *barcode* usage can be obtained by invoking the --help option::

    jsa.np.barcode --help

Input
=====
 *barcode* takes three files as required input::

	jsa.np.barcode -seq <nanopore reads> -bc <barcode.fasta> -script <analysis_script>

<*nanopore reads*> is either the long reads in FASTA/FASTQ file (after MinION sequencing is 
finished) or standard input ( specified by "-", for real-time analysis). 
	
<*barcode.fasta*> is the FASTA file of barcode sequences (given by ONT) with name correspond to the assigned sample id.

<*analysis_script*> is the script call for further action on the de-multiplexed reads. It always take one argument and be
executable by invoking::

	./analysis_script <id>
	
in which <*id*> is the identifier of a sample as given in the <barcode.fasta>. The script should read the standard input
of long-read streams to do further analysis.
	
	Missing any file would break down the whole pipeline.
	
	*barcode* allows user to set the minimum score of a hit with barcode reference to be considered valid. The default value
is 24, the length of a single barcode sequence from Nanopore Native Kit. Decreasing the threshold will lead to more reads being
clustered but with higher risk of false positive while increasing will generate less but more confident of demultiplexed
reads.

Output
======
*barcode* output depends on the <*analysis script*> because the de-multiplexed reads are streamed directly to its dedicated process.
If ones only interest in de-multiplexing alone, then the script should be as simple as to write stream to file. For example:

.. code-block:: bash
	:linenos:
	
	#!/bin/bash
	while read line
	do
	  echo "$line"
	done >> ${1}_script.fasta

This is equivalent to enable the *-p* option::

	jsa.np.barcode -seq <nanopore reads> -bc <barcode.fasta> -script <analysis_script> -p
	
that would print out de-multiplexed FASTA sequences <id>\_clustered.fasta	

Real-time scaffolding for barcode sequencing
============================================

One use-case for barcode sequencing is to run *npscarf* on the resulted de-multiplexed reads. This could be done by calling a script 
that can take an output folder of long reads from a sample to scaffold its corresponding short-reads (e.g. SPAdes) assembly.
E.g.

.. code-block:: bash
	:linenos:
	
	#!/bin/bash
	dirname=`find /coin/barcode/ -maxdepth 1 -type d -name "*${1}*" -print -quit`
	
	bwa index ${dirname}/contigs.fasta
	
	bwa mem -t 16 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y -K 10000 ${dirname}/contigs.fasta - 2> /dev/null | \
	jsa.np.npscarf -realtime -read 100 -time 1 -b - -seq ${dirname}/contigs.fasta -spadesDir ${dirname} -prefix ${1} > ${1}.log 2>&1
	
In this scenario, we assume the output SPAdes folders locate in one directory and the folder names contain the ID of the corresponding samples.
 *RST*/