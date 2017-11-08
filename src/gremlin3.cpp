#include <cstdio>
#include <unistd.h>
#include <ctime>
#include <cstring>
#include <cmath>

#include <string>
#include <map>
#include <vector>
#include <utility>

#include "MSAclass.h"
#include "ProblemFull.h"
//#include "ProblemRRCE.h"
//#include "ProblemPNAS.h"
#include "Minimizer.h"
#include "MRFprocessor.h"
#include "ContactList.h"
#include "RRCE.h"

/* TODO: separate programs
 * 1) gremlin3 - produces an MRF
 *    1a. gremlin3c - for complexes: block-APC
 * 2) score_patchdock
 * 3) score_gramm
 * 4) score_decoy */

/* TODO:
 *   + TEST - convert 'shift' into cleaned MSA numbering
 *   + TEST - add block-APC (BAPC) correction & parameter '-b'
 *   + TEST - output per group scores (for constrained problem only)
 *   - GREMLIN2 symmetric minimizer
 *   - don't store masked edges (simplified problem with less variables)
 *   - GREMLIN1 routine with L1 penalty (ADMM solver needed)
 */

struct OPTS {
	char *a3m; /* A3M file */
	char *mtx; /* file with the computed contact matrix */
	char *mrf; /* file to save MRF */
	size_t niter; /* number of iterations */
	double grow; /* gaps per row */
	double gcol; /* gaps per col */
	char *mask;
	char *umask;
	int rmode; /* regularization mode */
//	size_t shift; /* block A for block-APC */
};

bool GetOpts(int argc, char *argv[], OPTS &opts);
void PrintOpts(const OPTS &opts);

double PairEnergies(const MSAclass &MSA, double **mtx);

int main(int argc, char *argv[]) {

	/*
	 * (0) process input parameters
	 */
	OPTS opts = { NULL, NULL, NULL, 25, 0.25, 0.25, NULL, NULL, 2 };
	if (!GetOpts(argc, argv, opts)) {
		PrintOpts(opts);
		return 1;
	}

	/*
	 * (1) read & clean MSA
	 */
	MSAclass MSA(opts.a3m);
	MSA.CleanMsa(opts.grow, opts.gcol);
	MSA.Reweight();

	/*
	 * (2) calculate basic MSA statistics
	 *     and save it
	 */
	size_t ncol = MSA.GetNcol();
	ContactList Contacts(ncol);
	Contacts.AddTerm( { "log(Neff)", log(MSA.GetNeff()) });
	Contacts.AddTerm( { "log(Ncol)", log(ncol) });
	double **mtx = (double**) malloc(ncol * sizeof(double*));
	for (size_t i = 0; i < ncol; i++) {
		mtx[i] = (double*) malloc(ncol * sizeof(double));
	}

	/* mutual information */
	MSA.MI(mtx);
	Contacts.AddFeature("MI", mtx);

	/* joint entropy */
	MSA.Hxy(mtx);
	Contacts.AddFeature("Hxy", mtx);

	/* statistical potential */
	printf("# E(RRCE)= %.5f\n", PairEnergies(MSA, mtx));
	Contacts.AddFeature("RRCE", mtx);

	/*
	 * (3) read edge constraints
	 */
	EDGE_LIST L;
	if (opts.mask != NULL) {
		L = ReadEdges(opts.mask);
	} else if (opts.umask != NULL) {
		L = ReadEdges(opts.umask);
	}
	for (EDGE_LIST::iterator it = L.begin(); it != L.end(); it++) {
		printf("# %lu edges in group '%s'\n", it->second.size(),
				it->first.c_str());
	}

	/*
	 * (3) set up the problem
	 */
	/*
	 ProblemPNAS P(MSA, 2);
	 ProblemRRCE P(MSA, 7);
	 */
	ProblemFull P(MSA);
	if (opts.mask != NULL) {
		P.UnmaskAllEdges();
		for (EDGE_LIST::iterator it = L.begin(); it != L.end(); it++) {
			P.MaskEdges(it->second);
		}
	} else if (opts.umask != NULL) {
		P.MaskAllEdges();
		for (EDGE_LIST::iterator it = L.begin(); it != L.end(); it++) {
			P.UnmaskEdges(it->second);
		}
	}

	/*
	 * (4) solve P
	 */
	MRFclass MRF = Minimizer::MinimizeLBFGS(P, opts.niter);

	/*
	 * (5) save MRF
	 */
	if (opts.mrf != NULL) {
		MRF.Save(opts.mrf);
	}

	/*
	 * (6) do MRF post-processing
	 */
	if (opts.mtx != NULL) {
		MRFprocessor::MTX result;
		switch (opts.rmode) {
		case 1:
			printf("# Contact matrix correction: FN\n");
			MRFprocessor::FN(MRF, result);
			break;
		case 2:
			printf("# Contact matrix correction: APC\n");
			MRFprocessor::APC(MRF, result);
			break;
		case 3:
			printf("# Contact matrix correction: PROB\n");
			MRFprocessor::APC(MRF, result);
			MRFprocessor::APC(MRF, mtx);
			MRFprocessor::Zscore(ncol, mtx);
			Contacts.AddFeature("Z(APC)", mtx);
			break;
		default:
			printf("!!! ACHTUNG !!! (this should never happen)\n");
			return 1;
		}
		MRFprocessor::SaveMTX(result, opts.mtx);

		/*
		 * (7) calculate per group scores
		 */
		if (opts.umask != NULL) {
			for (EDGE_LIST::iterator it = L.begin(); it != L.end(); it++) {
				std::vector<std::pair<size_t, size_t> > vec = MSA.CastToMsa(
						it->second);
				double score = MRFprocessor::GetScore(result, vec);
				double energy = MRF.GetPairEnergies(MSA, vec);
				printf("# Score(%s)= %.6e, Energy(%s)= %.6e\n",
						it->first.c_str(), score, it->first.c_str(), energy);
			}
		}
	}

	/*
	 *
	 */
	Contacts.Print();

	/*
	 * free
	 */
	for (size_t i = 0; i < ncol; i++) {
		free(mtx[i]);
	}
	free(mtx);

	return 0;

}

bool GetOpts(int argc, char *argv[], OPTS &opts) {

	char tmp;
	while ((tmp = getopt(argc, argv, "hi:o:f:n:r:c:m:u:R:")) != -1) {
		switch (tmp) {
		case 'h': /* help */
			printf("!!! HELP !!!\n");
			return false;
			break;
		case 'i': /* A3M file (in) */
			opts.a3m = optarg;
			break;
		case 'o': /* MTX file (out) */
			opts.mtx = optarg;
			break;
		case 'f': /* MRF file (out) */
			opts.mrf = optarg;
			break;
		case 'n': /* number of iterations */
			opts.niter = atoi(optarg);
			if (!opts.niter) {
				printf("Error: niter shoud be > 0\n");
				return false;
			}
			break;
		case 'r': /* gaps per row */
			opts.grow = atof(optarg);
			if (opts.grow < 0.0 || opts.grow >= 1.0) {
				printf("Error: 0.0 <= gaps_per_row < 1.0 (-r)\n");
				return false;
			}
			break;
		case 'c': /* gaps per column */
			opts.gcol = atof(optarg);
			if (opts.gcol < 0.0 || opts.gcol >= 1.0) {
				printf("Error: 0.0 <= gaps_per_col < 1.0 (-c)\n");
				return false;
			}
			break;
		case 'm': /* a list of residue pairs to be masked */
			opts.mask = optarg;
			break;
		case 'u': /* a list of residue pairs to be unmasked */
			opts.umask = optarg;
			break;
//		case 'b':
//			opts.shift = atoi(optarg);
//			if (!opts.shift) {
//				printf("Error: block size should be > 0 (-b)\n");
//				return false;
//			}
//			break;
		case 'R': /* regularization mode */
			if (strcmp(optarg, "FN") == 0) {
				opts.rmode = 1;
			} else if (strcmp(optarg, "APC") == 0) {
				opts.rmode = 2;
			} else if (strcmp(optarg, "PROB") == 0) {
				opts.rmode = 3;
			} else {
				printf("Error: wrong matrix correction mode '%s'\n", optarg);
				return false;
			}
			break;
		default:
			return false;
			break;
		}
	}

	if (opts.a3m == NULL) {
		printf("Error: A3M file not specified\n");
		return false;
	}

//	if (opts.rmode == 3 && !opts.shift) {
//		printf(
//				"Error: block-APC correction requires setting block size '-b'\n");
//		return false;
//	}
//
//	if (opts.shift && opts.rmode == 2) {
//		printf(
//				"Error: setting block size requires FN or BAPC correction '-R'\n");
//		return false;
//	}

	return true;

}

void PrintOpts(const OPTS &opts) {

	printf("Usage:   ./gremlin3 [-option] [argument]\n");
	printf("Options:  -i alignment.a3m (input, required)\n");
	printf("          -o matrix.txt (output)\n");
	printf("          -f mrf.txt (output)\n");
	printf("          -n number of iterations (%ld)\n", opts.niter);
	printf("          -r gaps per row [0;1) (%.2lf)\n", opts.grow);
	printf("          -c gaps per column [0;1) (%.2lf)\n", opts.gcol);
	printf("          -m list1.txt - residue pairs to be masked\n");
	printf("          -u list2.txt - residue pairs to be unmasked "
			"(all others are masked)\n");
	printf("          -R contact matrix correction {FN,APC,PROB} (APC)\n");
//	printf("          -b block size for the block-APC correction\n");

}

double PairEnergies(const MSAclass &MSA, double **mtx) {

	double E = 0.0;

	RRCE RRCE_(RRCE::RRCE20RC, 7.8, 5);

	size_t nrow = MSA.GetNrow();
	size_t ncol = MSA.GetNcol();

	unsigned char * msa = MSA.GetMsa();
	MSAclass::aatoi(msa, nrow * ncol);

	for (size_t i = 0; i < ncol; i++) {
		memset(mtx[i], 0, ncol * sizeof(double));
	}

	for (size_t i = 0; i < nrow; i++) {

		unsigned char *seq = msa + i * ncol;
		double w = MSA.GetWeight(i);

		for (size_t p = 0; p < ncol; p++) {
			unsigned char a = seq[p];
			if (a < 0 || a >= 20) {
				continue;
			}

			for (size_t q = p + 1; q < ncol; q++) {

				unsigned char b = seq[q];
				if (b < 0 || b >= 20) {
					continue;
				}

				double j = RRCE_.GetJij(a, b);
				mtx[p][q] += j * w;
				E += j;

			}

		}

	}

	for (size_t p = 0; p < ncol; p++) {
		for (size_t q = p + 1; q < ncol; q++) {
			mtx[p][q] /= nrow;
			mtx[q][p] = mtx[p][q];
		}
	}

	return 2.0 * E / (ncol - 1) / ncol;

}
