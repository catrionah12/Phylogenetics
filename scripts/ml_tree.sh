for file in ../phylogeny/temp/filtered*.nexus
do
	../external/iqtree2 -s "$file" -B 1000 -redo
done

cd ../phylogeny/temp
for file in *treefile
do
	mv "$file" ../trees/"$file"
done
