#!/bin/bash
build_dir="../../build2/"

echo "Generating simulation results directories..."
mkdir n_5_perfect
mkdir n_5_noisy
mkdir n_5_inf_alleles_violations
mkdir n_15_noisy

echo "Running instances in n_5_perfect"
cd n_5_perfect
n=5
for r in {0..19}
do
	for m in {2,5,10}
	do
		f=sims_r${r}_m${m}_n${n}
		echo "Running $f"
		$build_dir/enumerate -v 0 ../../data/sims/n_5_perfect/${f}.data > ${f}.res
	done
done
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
			echo "Running $f"
			$build_dir/enumerate -v 0 ../../data/sims/n_5_noisy/${f}.data > ${f}.res
		done
	done
done
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
		$build_dir/enumerate -v 0 ../../data/sims/n_5_inf_alleles_violations/${f}.data > ${f}.res
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
	$build_dir/enumerate -v 0 ../../data/sims/n_15_noisy/${f}.data > ${f}.res
done
cd ..
