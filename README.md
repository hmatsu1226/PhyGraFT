# PhyGraFT
PhyGraFT: a network-based method for phylogenetic trait analysis

## Reference

## Requirements
PhyGraft is written with R and uses RSpectra library for eigenvalues and vectors calculation.

## Functions
PhyGraft includes following functions

- sym_normalized_graph_laplacian(A)
	- Return symmetrically normalized graph Laplacian matrix from the adjacency matrix of A.
- graph_laplacian_regularizer(x, U, eigenvalues)
	- Return aGLR value of phlygenetic signal x using eigenvector matrix U and eigenvalues.
- graph_fourier_transform(L,X,m)
	- Return the top-m eigenvectors and eigenvalues of L, graph Fourier coefficients hF and tF (haf_f and tilde_f) for trait data matrix X.
- barplot_gfdomain(hF,tF,i)
	- Plot barplot for GF coefficients for the i-th eigenvector.
- plot_aGLR(X,U,shuffle_num)
	- Plot aGLR values of trait data matrix X. shuffle_num is the number of negative control data generated from shuffling X.

## Dataset
We have four datasets, and each dataset have the ajdacency matrix of k-NNG (A.txt) and binary trait data matrix (X.txt).

- dataset_HA
	- The dataset based on the HA sequences of influenza type A virus.
- dataset_NA
	- The dataset based on the NA sequences of influenza type A virus.
- dataset_PB2
	- The dataset based on the PB2 sequences of influenza type A virus.
- dataset_virome
	- The dataset based on the virome genome gene-sharing data.

## Example Jupyter Notebook
The analysis_*.ipynb are the jupyter notebook for each dataset analysis.
