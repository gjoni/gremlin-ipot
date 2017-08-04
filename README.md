# gremlin3
Experimental GREMLIN version with RRCE statistical potentials [1]

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
[1] I Anishchenko, PJ Kundrotas, IA Vakser. Contact energies in proteins and protein complexes inferred from the Potts model. (2017) In preparation

