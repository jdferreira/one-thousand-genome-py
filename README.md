# VCF-machines

This project is a first attempt to create a set of `VCF-machines`, programs that read from a VCF file and output computations based on the data in those files. 

A VCF-machine can be a comparer of individuals, outputting a numeric value for each pair of individuals in a VCF file (or a subset of them), or a population predictor, outputting a population label for each individual based on how it compares with the rest of the individuals.

For now, there is code to download and massage the 1kone-sthousand genome project, which contains some data that we use to test our methods.
