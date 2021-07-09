#!/bin/sh

#SBATCH --job-name=jST
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8000 # mb
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
export JSA_MEM=7800m

export japsa_coverage="${HOME}/github/japsa_coverage"
echo ${japsa_coverage}


mainclass="japsa.bio.phylo.MergeKrakenCmd"
echo $mainclass
dat=$(date +%Y%m%d%H%M%S)
resdir="results_${dat}"


opts="--todo=todo.txt --output0=cumul.css --output1=sep.css   --thresh=0"


if [ $1 ]; then
		find . -name results.krkn | grep blast  | grep $1 > todo.txt
else
	        find . -name results.krkn | grep blast > todo.txt
fi
if [ ! -s cumul.css ]; then 
bash ${japsa_coverage}/scripts/run.sh ${mainclass}  ${opts}
fi
echo "R CMD BATCH ${japsa_coverage}/R/visColors.R"
R CMD BATCH  ${japsa_coverage}/R/visColors.R
#echo $cmd
#pwd=$(pwd)
#for i in resdir/*; do
#	cd $i
#	nme=$(pwd | rev | cut -f 1 -d / | rev)
#	ls | xargs -I {} mv {} ${nme}.{}
#	cd $pwd
#done

