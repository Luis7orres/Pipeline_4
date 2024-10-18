#Mapear los reads recortados contra las secuencias consenso
minimap2 -ax map-ont ./swarm/swarm_representatives.fasta recortadas.fastq > mapped_reads.sam

#Convertir el archivo SAM a BAM y ordenarlo
samtools view -bS mapped_reads.sam | samtools sort -o mapped_reads_sorted.bam

#Generar un archivo de índice BAM
samtools index mapped_reads_sorted.bam

#Generar una tabla de cuentas
samtools idxstats mapped_reads_sorted.bam > counts_table.txt

#Generar un report de mapeo
samtools flagstat mapped_reads_sorted.bam > report.txt