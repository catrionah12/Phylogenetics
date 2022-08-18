cd ../genomes/evolution_rates/ids

for file in OG*txt
do
	while read id; do
		for fastafile in ../../data/cds/*.dna.fasta
		do
			grep -A 1 "$id" "$fastafile"  >> ../evoltree/"$file".fasta
		done
		for fastafile in ../../data/proteins/sequences/*prot.fasta
		do
			grep -A 1 "$id" "$fastafile" >> ../evoltree/"$file"_prot.fasta
		done		
	done<"$file"
done

cd ../evoltree

for file in *prot.fasta
do
	../../../external/mafft-mac/mafft.bat --auto "$file" > aln_"$file"
done

cd ../ids

for file in *.txt
do
	../../../external/pal2nal.pl ../evoltree/aln_"$file"_prot.fasta ../evoltree/"$file".fasta -output paml > ../evoltree/codon_aln/aln_"$file".paml
	rm ../evoltree/"$file"_prot.fasta
	rm ../evoltree/aln_"$file"_prot.fasta
	rm ../evoltree/"$file".fasta
done

cd ../evoltree/codon_aln
for file in aln*.paml
do
	sed -e 's/PSC.*\.1/m_conductrix_genome_proteins/' -e 's/PRW.*\.1/c_sorokiniana_genome_proteins/' -e 's/XP.*\.1/c_variabilis_genome_proteins/' -e 's/sca.*1/A99_genome_proteins/' -e 's/KAI.*\.1/c_vulgaris_genome_proteins/' "$file" > sp_name_"$file"
	rm "$file"
done

cd ../../ids
for file in OG*txt
do
	python3 ../../../scripts/evoltree.py "$file" "$1"
done
