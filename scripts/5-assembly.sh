#Illumina
spades.py --careful -1 reads_R1.fastq -2 reads_R2.fastq -o spades_output --reference reference.fasta

#Nanopore (no incluye la opción de genoma referencia, pero siendo secuencias largas no hace falta)
spades.py --nanopore -s reads.fastq -o spades_output