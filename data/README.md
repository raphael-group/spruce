# Input format

SPRUCE takes two files as input: a data file and an (optional) interval file.

## Data file format
The data file contains the VAF confidence intervals and copy-number mixing proportions for each SNV in each sample. 
The first line of this file contains the number of samples; the next line contains the number of characters (SNVs). 
The third line is a header that describes the (fixed) format of the subsequent lines.

	10 #m
	82 #n
	#sample_index	sample_label	character_label	character_index	vaf_lb	vaf_mean	vaf_ub	x	y	mu	x	y	mu
	0	A22-A	0	ENSG00000214793	0	0	0	1	1	1	2	1	0
	1	A22-C	0	ENSG00000214793	0	0	0	1	1	1	2	1	0
	2	A22-D	0	ENSG00000214793	0.139700581	0.190730838	0.250219406	1	1	0.718981296	2	1	0.281018704

Requirements:

* `sample_index` and `character_index` are numbered starting from 0
* `sample_label` and `character_label` do not contain whitespaces or tabs
* `vaf_mean`, `vaf_lb` and `vaf_ub` are between 0 and 1 such that `vaf_lb <= vaf_mean <= vaf_ub`
* `mu` values (the proportion of cells with copy-number state `(x,y)`) are non-negative and sum up to 1 (for each row)

See [`real/A22.tsv`](real/A22.tsv) for an example.
Please see [`convert.py`](convert.py) for a Python script for obtaining confidence intervals from read count data.

## Interval file format

The optional interval file relates SNVs that are affected by the same copy-number aberration in a sample. 
Each line of this file corresponds to a sample (ordered by `sample_index`). The SNV clusters are separated by tab characters. Each SNV in a cluster is separated by a single whitespace character. An SNV cluster is denoted by its corresponding `character_index`.

See [`real/A22.intervals`](real/A22.intervals) for an example.
