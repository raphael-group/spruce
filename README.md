# SPRUCE
SPRUCE (Somatic Phylogeny Reconstruction using Combinatorial
Enumeration) is an algorithm for inferring the clonal evolution of
single-nucleotide and copy-number variants given multi-sample bulk tumor
sequencing data.

## Support

For support using SPRUCE, please visit the [SPRUCE Google Group](https://groups.google.com/forum/#!forum/sprucealgorithm).

## Dependencies

SPRUCE is written C++. In addition to a recent C++ compiler (that supports C++11), it has the following dependencies:

* [CMake](http://www.cmake.org/) (>= 2.8)
* [Boost](http://www.boost.org) (>= 1.38)
* [LEMON](http://lemon.cs.elte.hu/trac/lemon) graph library (>= 1.3)

[Graphviz](http://www.graphviz.org) is required to visualize the resulting DOT files, but is not required for compilation.

## Compilation instructions

To compile SPRUCE, execute the following commands from the root of the repository:

    mkdir build
    cd build
    cmake ..
    make
    
In case CMake fails to detect LEMON, run the following command with adjusted paths:

	cmake -DLIBLEMON_ROOT=~/lemon ..
	
The compilation results in the following files in the `build` directory:

EXECUTABLE | DESCRIPTION
-----------|-------------
`cliques`  | Enumerates cliques of the compatibility graph (given a size and a filter)
`enumerate`| Enumerates perfect phylogeny trees
`merge`    | Merges multiple solution files into one
`rank`     | Sorts solution trees by the fraction of common edges (solution with rank 0 is the most representative tree)
`visualize`| Visualizes one solution or the entire solution space


For example usage see [`result/run_A22.sh`](result/run_A22.sh) and corresponding [instructions](result/README.md).
For a description of the input file format see [`data/README.md`](data/README.md).