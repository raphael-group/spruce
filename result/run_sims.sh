#!/bin/bash
#MEK, 06/28/2016
build_dir="../../build/"

echo "Generating simulation results directories..."
if [ ! -d "n_5_perfect" ]; then mkdir n_5_perfect; fi
if [ ! -d "n_5_noisy" ]; then mkdir n_5_noisy; fi
if [ ! -d "n_5_inf_alleles_violations" ]; then mkdir n_5_inf_alleles_violations; fi
if [ ! -d "n_15_noisy" ]; then mkdir n_15_noisy; fi

echo "Running instances in n_5_perfect"
cd n_5_perfect
n=5
for r in {0..19}
do
	for m in {2,5,10}
	do
		f=sims_r${r}_m${m}_n${n}
		echo "Running $f..."
		$build_dir/enumerate -p -v 0 ../../data/sims/n_5_perfect/${f}.data > ${f}.res
		$build_dir/rank $f.res | ($build_dir/recall - ../../data/sims/n_5_perfect/${f}.true 2> /dev/null) | head -n 1 > ${f}.recall
		$build_dir/recall $f.res ../../data/sims/n_5_perfect/${f}.true 2> ${f}.median.recall > /dev/null
	done
done
python ../process_perfect.py > ../n_5_perfect.csv
cd ..

echo "Running instances in n_5_noisy"
cd n_5_noisy
n=5
for r in {0..19}
do
	for m in {2,5,10}
	do
		for c in {50,100,500,1000,10000}
		do
			f=sims_r${r}_m${m}_n${n}_c${c}
			echo "Running $f..."
			$build_dir/enumerate -p -v 0 ../../data/sims/n_5_noisy/${f}.data > ${f}.res
			$build_dir/rank $f.res | ($build_dir/recall - ../../data/sims/n_5_noisy/${f}.true 2> /dev/null) | head -n 1 > ${f}.recall
		$build_dir/recall $f.res ../../data/sims/n_5_noisy/${f}.true 2> ${f}.median.recall > /dev/null
		done
	done
done
python ../process_noisy.py > ../n_5_noisy.csv

cd ..

echo "Running instances in n_5_inf_alleles_violations"
cd n_5_inf_alleles_violations
m=10
n=5
c=1000
for v in {0..5}
do
	for r in {0..19}
	do
		f=sims_r${r}_m${m}_n${n}_c${c}_v${v}
		echo "Running $f..."
		$build_dir/enumerate -p -v 0 ../../data/sims/n_5_inf_alleles_violations/${f}.data > ${f}.res
	done
done
cd ..

echo "Running instances in n_15_noisy"
cd n_15_noisy
m=10
n=15
c=1000
for r in {0..19}
do
	f=sims_r${r}_m${m}_n${n}_c${c}
	echo "Running $f..."
	$build_dir/enumerate -p -v 0 ../../data/sims/n_15_noisy/${f}.data > ${f}.res
	$build_dir/rank $f.res | ($build_dir/recall - ../../data/sims/n_15_noisy/${f}.true 2> /dev/null) | head -n 1 > ${f}.recall
	$build_dir/recall $f.res ../../data/sims/n_15_noisy/${f}.true 2> ${f}.median.recall > /dev/null
done
python ../process_noisy.py > ../n_15_noisy.csv
cd ..
