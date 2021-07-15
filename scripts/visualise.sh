#!/bin/sh

#SBATCH --job-name=jST
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8000 # mb
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output=jst_%A_%a.out
#SBATCH --error=jst_%A_%a.err
n=$SLURM_ARRAY_TASK_ID
        fastq=$(cat todo.txt | tail -n +$n | head -n 1)
mkdir -p res
java -Xmx8000m -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8 -classpath ~/github/japsa_coverage/target/japsacov-1.9.5e.jar japsa.bio.phylo.MergeKrakenCmd --input=$fastq --output=res/${fastq}
