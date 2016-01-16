# OrthoClust

OrthoClust is a clustering algorithm built on a multilayer network framework. It concatenates networks from individual species by their orthology relationships, arriving at a multiplex network. By optimizing the a generalized modularity function, OrthoClust returns a set of modules that could be either conserved or species-specific.

INPUT FILES:

There are two kinds of input files: the network files and a coupling information file.
A network file is a tab delimited file which consists of two columns of integers. They are simply indices of the nodes, running from 1 to the size of the network (see data/net_fly as an example). Depending on the number of species, there should be the same number of network files.
A coupling information file is a tab delimited file which consists of five columns. The 1st and 2nd columns are indices of the species (running from 1 to the number of species), the 4th and 5th columns are the indices of genes (the same as those specified in the network files). The 3rd column is the coupling constant between networks specified by the 1st and 2nd columns. It is a parameter which could adjusted. In principle, if there are N species, there could be N(N-1)/2 free coupling constants. e.g. If the coupling between species "1" and "2" is 1, and if gene "17" and gene "2016" of species "1" and species "2" are orthologs, the numbers 1,2,1,17,2016 will show up in column 1,2,3,4,5 respectively (see data/ortho_info as an example). 

OUTPUT FILE:

The output file is a tab delimited file which consists of three columns. The 1st and 2nd column are the species id and the gene id given by the input files. The 3rd column is a module id. Nodes with the same id are assigned to the same module.

USAGE:

OrthoClust is written in Julia. It has been tested in Julia v0.4.0, it should work for earlier versions. If Julia and the required packages are installed (see the first few lines in lib_orthoclust.jl), one could simply run in the command prompt

> julia run_orthoclust.jl net1 net2 ortho_info

The minimal number of network files is 2. If there are more than 2 networks, simply use, for instance,

> julia run_orthoclust.jl net1 net2 net3 net4 ortho_info

If you prefer to run OrthoClust in a Julia interactive shell, see demo.jl and lib_orthoclust.jl

The current implementation is quite time consuming. The test data, which consists of ~35000 genes, take about 3 hour to finish in a single node.

FURTHER DETAILS:

The framework and the objective function of OrthoClust are explained in Yan et al. Genome Biology 2014 (see reference below). The optimization was implemented by simulated annealing in the paper. The Julia version employs a heuristic (the Louvain algorithm) to speed up the process. The algorithm has a stochastic procedure, one could repeat the algorithm for more robust results. We may include auxillary functions like visualization tools later.

While the concept was developed for integrating networks across species by ortholgy relationships, it could be employed in other contexts if the coupling information file is properly defined. 

REFERENCE:

OrthoClust: an orthology-based network framework for clustering data across multiple species.
Yan et al. Genome Biology 2014, 15:R100
