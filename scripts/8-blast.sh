#Convertir la referencia (k-mers) en una base de datos compatible con blast
makeblastdb -in reference.fasta -dbtype nucl -out reference_db

#Blast
blastn -query ./swarm/swarm_representatives.fasta -db reference_db -out blast_output.txt -outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"