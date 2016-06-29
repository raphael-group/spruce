# Results

In the following we describe how to obtain the results presented in the paper.

## Simulations

We consider the following simulated instances:

* [`data/sims/n_5_perfect`](../data/sims/n_5_perfect)
* [`data/sims/n_5_noisy`](../data/sims/n_5_noisy)
* [`data/sims/n_5_inf_alleles_violations`](../data/sims/n_5_inf_alleles_violations)
* [`data/sims/n_15_noisy`](../data/sims/n_15_noisy)

Run the following commands from the repository root directory to obtain the simulation results:

	cd build
	make recall
    cd ../result
    ./run_sims.sh
    
The simulation results (including #solutions and recall values) will be in:

* `result/n_5_perfect.csv`
* `result/n_5_noisy.csv`
* `result/n_15_noisy.csv`

## A22

Run the following commands from the repository root directory to obtain the A22 results:

    cd result
    ./run_sims.sh
    ./run _A22.sh
    
The A22 results will be in:

* `result/A22/A22_offset_70_0_666_3249.pdf` is a visualization of the phylogenetic tree corresponding to state tree S_1 for *PPFIA1*
* `result/A22/A22_offset_70_0_666_summary.pdf` is a visualization of the enumerated solution space of using state tree S_1 for *PPFIA1*
* `result/A22/A22_offset_70_1_1136_95.pdf` is a visualization of the phylogenetic tree corresponding to state tree S_2 for *PPFIA1*
* `result/A22/A22_offset_70_1_1136_summary.pdf` is a visualization of the enumerated solution space of using state tree S_2 for *PPFIA1*
