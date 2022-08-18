cd ../phylogeny/data/fasta_files
for file in *fasta
do
	../../../external/mafft-mac/mafft.bat --auto "$file" > ../../temp/aln_"$file" 
done

cd ../../temp
for file in aln*.fasta
do
	java -jar ../../external/BMGE/BMGE.jar -i "$file" -t DNA -of filtered_"$file" -oh filtered_"$file".html -g 0.1 -on filtered_"$file".nexus
done
