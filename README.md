# gremlin3
Experimental GREMLIN version with residue pair preferences from [1] and RRCE statistical potentials [2]

## TODO

* OpenMP parallelization
* (+) remove all GSL-related code
* generate Rosetta constraints
* PROB8 - default correction
* remove lpair/lskew options
* hide masking/unmasking options
* separate program for Neff calculation
* map_align input
* clean ./data folder, add relevant stuff, introduce GREMLINDAT env. variable

## Installation

### Download and compilation
    git clone https://github.com/gjoni/gremlin3
    cd ./gremlin3
    make

## Issues

* The program is not likely to work on 32bit machines; sse2 support is required
* OpenMP support not tested

## Acknowledgements

This package uses L-BFGS minimizer by J.Nocedal and N.Okazaki available at https://github.com/chokkan/liblbfgs.

Original GREMLIN2 protocol (aka plmDCA) is substantially adopted from the CCMpred package by J.Soeding group https://github.com/soedinglab/CCMpred.

## References

[1] I Anishchenko, S Ovchinnikov, H Kamisetty, D Baker. Origins of coevolution between residues distant in protein 3D structures. (2017) 114(34):9122-7

[2] I Anishchenko, PJ Kundrotas, IA Vakser. Contact potential for structure prediction of proteins and protein complexes from Potts model. (2018) Submitted

