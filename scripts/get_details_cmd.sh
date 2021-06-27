blastdb="/DataOnline/Data/dbs/Blast/nt/nt"
/sw/blast/current/bin/blastdbcmd -db ${blastdb} -entry_batch <(samtools view blastn_sam.out | cut -f 1) -outfmt "%t	%a	%l	%T" -out blastn_sam.out.fa.index
