#!/bin/sh

#SBATCH --job-name=jRT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=15800 # mb
#SBATCH --time=100:00:00
#SBATCH --output=jst.stdout
#SBATCH --error=jst.stderr
#SBATCH --cpus-per-task=16


## Load required modules
#module load gcc/8.3.0
##module load samtools/1.9
#module load minimap2/2.17
#module load java

##script for running 
##tip - use symbolic link to put this in the directory with bam files
#run as sbatch run_slurm.sh species --bamFile=file.bam 
#  sbatch run_slurm_combined.sh human combined --RNA=false
export JSA_MEM=15200m

export japsa_coverage="${HOME}/github/japsa_coverage"
echo ${japsa_coverage}

bamfiles=$1

	mainclass="japsa.tools.bio.np.RealtimeResistanceGeneCmd"
	optsfile="opts_resistance.txt"
echo $mainclass

if [ ! -f $optsfile ]; then
  echo "need to copy ${optsfile}  into your working directory from ${japsa_coverage}/scripts/"
	exit;
fi

if [ ! $bamfiles ];then
	echo "need to define bamfile"
	exit;
fi
dat=$(date +%Y%m%d%H%M%S)
resdir="results_${dat}"

opts=$(grep -v '^#' ${optsfile})
echo $bamfiles
echo $typ
bash ${japsa_coverage}/scripts/run.sh ${mainclass} ${bamfiles} ${opts}

