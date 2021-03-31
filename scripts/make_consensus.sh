rm output.fasta
for i in *.fa_*; do 
~/abPOA-v1.0.1/bin/abpoa $i  | sed "s/Consensus_sequence/${i}/g" >> output.fasta
done
