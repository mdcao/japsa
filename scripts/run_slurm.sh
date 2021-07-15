#!/bin/sh

#SBATCH --job-name=jST
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64000 # mb
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output=jst_%A_%a.out
#SBATCH --error=jst_%A_%a.err


## Load required modules
#module load gcc/8.3.0
##module load samtools/1.9
#module load minimap2/2.17
#module load java

##script for running 
##tip - use symbolic link to put this in the directory with bam files
#run as sbatch --array=1:10 run_slurm.sh species
#todo.txt has paths to fastq

export JSA_MEM=61800m

export japsa_coverage="${HOME}/github/japsa_coverage"
echo ${japsa_coverage}
if [ ! -s 'todo.txt' ];then
 echo " need todo.txt"
 exit 
fi
dat=$(date +%Y%m%d%H%M)

resdir="results"


n=$SLURM_ARRAY_TASK_ID
        fastq=$(cat todo.txt | tail -n +$n | head -n 1)

fastqname=$(echo $fastq | rev | cut -d '/' -f 1 | rev) 
species=$1
bamfiles="--fastqFile=${fastq}"


if [ 'resistance' == $1 ]; then
	mainclass="japsa.tools.bio.np.RealtimeResistanceGeneCmd"
	optsfile="opts_resistance.txt"
else
	mainclass="japsa.tools.bio.np.RealtimeSpeciesTypingCmd"
	optsfile="opts_species.txt"
fi
resdir="${resdir}/${fastqname}"
mkdir -p ${resdir}
cp ${optsfile} ${resdir}
echo $mainclass

if [ ! -f $optsfile ]; then
  echo "need to copy ${optsfile}  into your working directory from ${japsa_coverage}/scripts/"
	exit;
fi

if [ ! $bamfiles ];then
	echo "need to define bamfile"
	exit;
fi

opts=$(grep -v '^#' ${optsfile})
echo $bamfiles
echo $typ
bash ${japsa_coverage}/scripts/run.sh ${mainclass} ${bamfiles} ${opts} --resdir=${resdir}
#echo $cmd
cd ${resdir}
rm consensus_output.fa
rm commontree.txt.css
 find . -name consensus_output.fa | xargs cat >> consensus_output.fa
 files="consensus_output.fa"
found=$files
 blastdb="/DataOnline/Data/dbs/Blast/nt/nt"
#run blast in loop
#echo running blast for ${found} in ${pwd};
#change into direcotry of found file
#cd $(dirname ${found});
rm commontree.txt.css
        #check it contains output and blastn not already run
        if [ ! -s blastn_sam.out ]; then
                if [ -s consensus_output.fa ]; then
                        /sw/blast/current/bin/blastn -num_threads 16 -outfmt 17 -query consensus_output.fa -db ${blastdb} -out blastn_sam.out;
                fi
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
      
