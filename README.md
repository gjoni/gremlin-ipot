# gremlin3
Experimental GREMLIN version with RRCE contact energies [2]

## TODO

* generate Rosetta constraints - separate program
* add ./example
* ./scripts folder content is irrelevant
* (+) OpenMP parallelization
* (+) remove all GSL-related code
* (+) PROB8 - default correction
* (+) remove lpair/lskew options
* (+) hide masking/unmasking options
* (+) separate program for Neff calculation
* (+) clean ./data folder, add relevant stuff, introduce GREMLINDAT env. variable
* (+) do we really need FN correction? hide it? HIDE!
* (+) matrix output after APC (for bbcontacts)

## Installation

### Download and compilation
```
git clone https://github.com/gjoni/gremlin3
cd ./gremlin3
make
```

### Setup
```
export GREMLINDAT=${INSTALL_DIR}/data
```

## Usage
```
Usage:   ./gremlin3 [-option] [argument]

Options:  -i alignment.a3m               - input, required
          -o matrix.txt                  - output, optional
          -b apcmatrix.txt               - output, optional
          -n number of iterations          (50)
          -r max gaps per row [0;1)        (0.25)
          -c max gaps per column [0;1)     (0.25)
          -R contact matrix correction
             {FN,APC,PROB5,PROB8}          (PROB8)
          -t number of threads             (1)

```

## Acknowledgements

 - [L-BFGS minimizer](https://github.com/chokkan/liblbfgs) by J.Nocedal and N.Okazaki
 - [ccmpred](https://github.com/soedinglab/CCMpred) by J.Soeding group

## References

[1] H Kamisetty, S Ovchinnikov, D Baker. Assessing the utility of coevolution-based residue-residue contact predictions in a sequence- and structure-rich era. PNAS (2013). 110:15674–9

[2] I Anishchenko, PJ Kundrotas, IA Vakser. Contact potential for structure prediction of proteins and protein complexes from Potts model. (2018) Submitted

