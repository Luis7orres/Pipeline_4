# Renombrar las columnas de idxstats si no tienen encabezado
echo -e "sseqid\tmapped_reads" > counts_table_renamed.txt
cat counts_table.txt >> counts_table_renamed.txt

#Formatear salida blast
cut -f 2,3,4,5,11 blast_output.txt > blast_output_processed.txt

#Combinar ambos archivos:
join -t $'\t' -1 1 -2 1 counts_table_renamed.txt blast_output_processed.txt > combined_output.txt