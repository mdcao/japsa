package japsa.tools.bio.np;

import java.io.IOException;

import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import japsa.bio.np.barcode.BarCodeAnalysis;


@Deployable(
		scriptName = "jsa.np.barcode", 
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
		addString("scriptRun", null, "Invoke command script to run npScarf");
		addDouble("threshold", 70, "Minimum identity(%) for barcode alignment");
		addDouble("distance", 4, "Minimum identity(%) distance between the best alignment to others");

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
		Double threshold = cmdLine.getDoubleVal("threshold"),
				distance = cmdLine.getDoubleVal("distance");

		BarCodeAnalysis.print = cmdLine.getBooleanVal("print");
		BarCodeAnalysis.twoends = cmdLine.getBooleanVal("twoends");

		
		BarCodeAnalysis bc = new BarCodeAnalysis(bcFile,script);
		bc.setThreshold(threshold);
		bc.setDistance(distance);
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

=====
Input
=====
 *barcode* takes 2 files as required input::

	jsa.np.barcode -seq <nanopore reads> -bc <barcode.fasta>

<*nanopore reads*> is either the long reads in FASTA/FASTQ file (after MinION sequencing is 
finished) or standard input ( specified by "-", for real-time analysis). 
	
<*barcode.fasta*> is the FASTA file of barcode sequences (given by ONT) with name correspond to the assigned sample id.

Missing any file would break down the whole pipeline.
	
In addition, one can provide <*analysis_script*> which is the script call for further action on the de-multiplexed reads. It always take one argument and be
executable by invoking::

	./analysis_script <id>
	
in which <*id*> is the identifier of a sample as given in the <barcode.fasta>. The script should read the standard input
of long-read streams to do further analysis.
	
*barcode* allows user to set the minimum criteria of a hit with barcode reference to be considered valid. The default value
is 70% for minimum identity. At the same time, 4% distance between the best hit and the second best is necessary for differentiation.
Decreasing the thresholds will lead to more reads being clustered but with higher risk of false positive while more stringent parameters
will generate less but more confident of demultiplexed reads.

User can also have control on the matching condition for barcode detection, either one-end match or both-end match. For the first case (default), only the 
a legal maximal hit from one end of a read is enough to label it while in the later case, we take into account a pair from both 5' and 3'terminus. 
Thus the input for each use case should be different. The one-end option can take the simple FASTA file of Nanopore barcodes while the two-end need pairs of
barcode to be specified (e.g. with _F and _R suffix). One of a typical use case for two-end matching is when we want to detect the super-barcode which includes
also tail- and primer-sequences in pre-defined orientation.

======
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

============================================
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