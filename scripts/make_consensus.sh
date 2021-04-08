rm output_combined.fasta
for j in *; do
	cd $j
	cnt=$(ls *.fa_* | wc -l)
	echo $j	
	echo $cnt
	 if [ $cnt -gt 0 ]; then
		if [ ! -s output.fasta ]; then
		for i in *.fa_*; do 
			~/abPOA-v1.0.1/bin/abpoa $i  | sed "s/Consensus_sequence/${i} ${j}/g" >> output.fasta
		done
		fi
		cat output.fasta >> ../output_combined.fasta
	fi
	cd ..
done
