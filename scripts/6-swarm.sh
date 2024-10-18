swarm -d 1 -i ./spades_output/contigs.fasta -o ./swarm/swarm_output.txt -s ./swarm/swarm_output.stats
vsearch --fastx_getseqs ./spades_output/contigs.fasta --labels $(awk '{print $1}' ./swarm/swarm_output.txt) --output ./swarm/swarm_representatives.fasta
