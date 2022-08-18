#ENA accession numbers

while read species rbcl ssu ITS 
do if [ "$ITS" != nd ]
then wget -O - $"https://www.ebi.ac.uk/ena/browser/api/fasta/$ITS" | sed "s/>ENA.*$/>$species /" >> ../phylogeny/data/fasta_files/ITS"$1".fasta; fi
if [ "$rbcl" != nd ]
then wget -O - $"https://www.ebi.ac.uk/ena/browser/api/fasta/$rbcl" | sed "s/>ENA.*$/>$species /" >> ../phylogeny/data/fasta_files/rbcl"$1".fasta; fi
if [ "$ssu" != nd ]
then wget -O - $"https://www.ebi.ac.uk/ena/browser/api/fasta/$ssu" | sed "s/>ENA.*$/>$species /" >> ../phylogeny/data/fasta_files/18s"$1".fasta; fi

done<../phylogeny/data/accession_numbers/"$1"accession_numbers.txt

#NCBI accession numbers

while read species ITS
do
esearch -db nuccore -query "$ITS" | efetch -format fasta | sed "s/>.*$/>$species /" >> ../phylogeny/data/fasta_files/ITS"$1".fasta; 
done<../phylogeny/data/accession_numbers/ncbi_accession_numbers.txt

#Extra sequences

if [ "$1" = 'sym' ]
then
	cat ../phylogeny/data/fasta_files/ITSsymbionts.fas  >> ../phylogeny/data/fasta_files/ITS"$1".fasta
else
	cat ../phylogeny/data/fasta_files/a99_its.fas >> ../phylogeny/data/fasta_files/ITS.fasta
fi
