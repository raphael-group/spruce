#!/bin/bash
build_dir="../../build/"

if [ ! -d "A22" ]
then
	echo "Generating A22 results directory"
	mkdir A22
fi

cd A22
f="A22_70_0.cliques"
if [ ! -e "$f" ]
then
	echo "Enumerating cliques for A22 (PPFIA1, S_1)..."
	$build_dir/cliques -s 26 -f 70,0 ../../data/real/A22.tsv > $f
fi

f="A22_70_1.cliques"
if [ ! -e "$f" ]
then
	echo "Enumerating cliques for A22 (PPFIA1, S_2)..."
	$build_dir/cliques -s 26 -f 70,1 ../../data/real/A22.tsv > $f
fi

if [ ! -e "A22_offset_70_0_666.merged.res.gz" ]
then
	echo "Running A22 (PPFIA1, S_1)..."
	$build_dir/enumerate ../../data/real/A22.tsv ../../data/real/A22.intervals -clique A22_70_0.cliques -ll 1200 -purity "0.8 0.675 0.772 0.89 0.808 0.704 0.763 0.672 0.798 0.817" -m -t 2 -lb 0 -s 1 -w 70 -o 69 -r 666 -v 2 | gzip > A22_offset_70_0_666.res.gz
	gzcat A22_offset_70_0_666.res.gz | $build_dir/rank - > A22_offset_70_0_666.merged.res
fi

if [ ! -e "A22_offset_70_1_1136.merged.res.gz" ]
then
	echo "Running A22 (PPFIA1, S_2)..."
	$build_dir/enumerate ../../data/real/A22.tsv ../../data/real/A22.intervals -clique A22_70_1.cliques -ll 1200 -purity "0.8 0.675 0.772 0.89 0.808 0.704 0.763 0.672 0.798 0.817" -m -t 2 -lb 0 -s 1 -w 70 -o 769 -r 1136 -v 2 | gzip > A22_offset_70_1_1136.res.gz
	gzcat A22_offset_70_1_1136.res.gz | $build_dir/rank - > A22_offset_70_1_1136.merged.res
fi

gzcat A22_offset_70_0_666.merged.res.gz | $build_dir/visualize -c ../../data/real/gundem_colors.txt -i 3249 -simple - > A22_offset_70_0_666_3249.dot
dot -Tpdf A22_offset_70_0_666_3249.dot -o A22_offset_70_0_666_3249.pdf

gzcat A22_offset_70_1_1136.merged.res.gz | $build_dir/visualize -c ../../data/real/gundem_colors.txt -i 95 -simple - > A22_offset_70_1_1136_95.dot
dot -Tpdf A22_offset_70_1_1136_95.dot -o A22_offset_70_1_1136_95.pdf

cd ..
