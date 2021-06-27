#!/bin/bash
#SBATCH --job-name=blast_consensus
#SBATCH --output=b_consensus.out
#SBATCH --error=b_consensus.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64000 #mb
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=16

#find output consensus fasta
export JSA_MEM=62800m
files=$(find jST/ -name consensus_output.fa -size +1b)
echo found `echo ${files} | wc -w` files
blastdb="/DataOnline/Data/dbs/Blast/nt/nt"
#run blast in loop
for found in ${files}; do
echo running blast for ${found};	
#change into direcotry of found file
cd $(dirname ${found});	
	#check it contains output
	if [ -s consensus_output.fa ]; then
	/sw/blast/current/bin/blastn -num_threads 16 -outfmt 17 -query consensus_output.fa -db ${blastdb} -out blastn_sam.out;
	fi
#check blast output contains data
if [ -s blastn_sam.out ];
then

#these commented out lines just produce a summary
#bash /DataOnline/Data/Projects/AQIP/get_details_cmd.sh;
#head blastn.out | cut -f 2 | tr "\n" "," | esummary -db nucleotide -id $(cat -) | xtract -pattern DocumentSummary -element Title > top_blast_hits.txt

#need index for java part
if [ -e hits.index.txt ]; then mv hits.index.txt blastn_sam.out.fa.index; fi
bash ~/github/japsa_coverage/scripts/get_details_cmd.sh 
#run java
echo creating results in $found
bash ~/github/japsa_coverage/scripts/run_blastspecies.sh --bamFile=blastn_sam.out
fi

#change back to root directory
cd -;
done
