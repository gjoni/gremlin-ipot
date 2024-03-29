# gremlin-ipot
GREMLIN version with residue-residue contact energies from the [iPot repo](https://github.com/gjoni/iPot).

## Installation

### Download and compilation
```
git clone https://github.com/gjoni/gremlin3
cd ./gremlin3
make
```

## Programs

* `gremlin` - predict protein contact map from an MSA
* `neff` - calculate effective number of sequences for an MSA
* `rstgen` - generate restraints for Rosetta


## Usage
```
Usage:   ./gremlin3 [-option] [argument]

Options:  -i alignment.a3m               - input, required
          -o matrix.txt                  - output, optional
          -b apcmatrix.txt               - output, optional
          -n number of iterations          50
          -r max gaps per row [0;1)        0.25
          -c max gaps per column [0;1)     0.25
          -R contact matrix correction
             {FN,APC,PROB5,PROB8}          PROB8
          -t number of threads             1

```

```
Usage:   ./rstgen [-option] [argument]

Options:  -i alignment.a3m               - input, required
          -m matrix.txt                  - input, required
          -o restraints.txt              - output, required
          -t restraint type {SIG, BND}     SIG
          -f fraction of top contacts      1.50 * Len
          -p probability cutoff            0.95
          -k sequence separation           3
```

## Acknowledgements

 - [L-BFGS minimizer](https://github.com/chokkan/liblbfgs) by J.Nocedal and N.Okazaki
 - [ccmpred](https://github.com/soedinglab/CCMpred) by J.Soeding group
 - [iPot](https://github.com/gjoni/iPot) statistical potential

## Links
- [GREMLIN](http://gremlin.bakerlab.org/)
- [Baker Lab](http://www.bakerlab.org/)
- [Ovchinnikov Lab](http://site.solab.org/home)


## References

[1] H Kamisetty, S Ovchinnikov, D Baker. Assessing the utility of coevolution-based 
residue-residue contact predictions in a sequence- and structure-rich era. 
[PNAS (2013). 110:15674–9](https://doi.org/10.1073/pnas.1314045110)

[2] S Ovchinnikov, L Kinch, H Park, Y Liao, J Pei, DE Kim, H Kamisetty, NV Grishin, D Baker. 
Large-scale determination of previously unsolved protein structures using evolutionary information. 
[eLife (2015). 4:e09248](https://doi.org/10.7554/eLife.09248)

[3] I Anishchenko, PJ Kundrotas, IA Vakser. Contact potential for structure prediction 
of proteins and protein complexes from Potts model. (2018) Biophys J. 
[DOI: 10.1016/j.bpj.2018.07.035](https://doi.org/10.1016/j.bpj.2018.07.035)

