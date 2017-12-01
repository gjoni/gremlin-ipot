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
#include "ProblemFullAsym.h"
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

/*
 * training set:      1,972 proteins from the Exceptions set
 * distance cut-off:  8A
 * minimal Neff:      10
 * datapoints used:   every 1, every 5, every 5
 * slogreg param A:   -9
 */

// Length (0:150)
const LogRegCoeff coef_range1 = { -1.67473, { 1.17148, -0.54004, -4.6156,
		-1.24735, 0.87817, 0.09511 }, { { 0.00000, 0.04282, 2.11214, 0.23992,
		0.07935, -0.27628 }, { 0.04282, 0.00000, -1.43714, -0.00709, 0.02930,
		0.00367 }, { 2.11214, -1.43714, 0.00000, -0.16130, 0.32952, 0.56220 }, {
		0.23992, -0.00709, -0.16130, 0.00000, 0.24220, 0.08363 }, { 0.07935,
		0.02930, 0.32952, 0.24220, 0.00000, -0.17670 }, { -0.27628, 0.00367,
		0.56220, 0.08363, -0.17670, 0.00000 } } };

// Length [150:300)
const LogRegCoeff coef_range2 = { 2.52312, { -0.87493, -0.6806, -2.8686,
		-1.38143, 0.41752, -0.78191 }, { { 0.00000, 0.11966, 2.18651, 0.19529,
		0.05010, 0.14904 }, { 0.11966, 0.00000, -1.43999, 0.02014, 0.03217,
		0.01518 }, { 2.18651, -1.43999, 0.00000, -0.21606, 0.34207, 0.22288 }, {
		0.19529, 0.02014, -0.21606, 0.00000, 0.26241, 0.07216 }, { 0.05010,
		0.03217, 0.34207, 0.26241, 0.00000, -0.07533 }, { 0.14904, 0.01518,
		0.22288, 0.07216, -0.07533, 0.00000 } } };

// Length [300:inf)
const LogRegCoeff coef_range3 = { 1.74862, { -0.94721, -0.801, -2.08548,
		-0.75452, 0.46976, -0.59321 }, { { 0.00000, 0.19966, 2.06230, 0.10662,
		0.02228, 0.15400 }, { 0.19966, 0.00000, -1.44987, 0.02691, 0.04843,
		0.01527 }, { 2.06230, -1.44987, 0.00000, -0.19650, 0.36617, 0.06135 }, {
		0.10662, 0.02691, -0.19650, 0.00000, 0.23712, -0.00996 }, { 0.02228,
		0.04843, 0.36617, 0.23712, 0.00000, -0.09214 }, { 0.15400, 0.01527,
		0.06135, -0.00996, -0.09214, 0.00000 } } };

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
	double lpair;
	double lskew;
};

bool GetOpts(int argc, char *argv[], OPTS &opts);
void PrintOpts(const OPTS &opts);

double PairEnergies(const MSAclass &MSA, double **mtx);
double PairEnergiesDI(const MSAclass &MSA, double **mtx);

int main(int argc, char *argv[]) {

	/*
	 * (0) process input parameters
	 */
	OPTS opts = { NULL, NULL, NULL, 25, 0.25, 0.25, NULL, NULL, 2, 0.2, 0.2 };
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

	/* temp storage for variour scores */
	size_t ncol = MSA.GetNcol();
	double **mtx = (double**) malloc(ncol * sizeof(double*));
	for (size_t i = 0; i < ncol; i++) {
		mtx[i] = (double*) malloc(ncol * sizeof(double));
	}

	/*
	 * (2) calculate basic MSA statistics
	 *     and save it
	 */
	ContactList Contacts(ncol);
	Contacts.AddTerm( { "log(Neff)", log(MSA.GetNeff()) });
	Contacts.AddTerm( { "log(Ncol)", log(ncol) });

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

//	ProblemFullAsym P(MSA);
//	P.SetL2(opts.lpair, opts.lskew);

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
	MRFclass MRF;
	if (opts.rmode > 0) {
		MRF = Minimizer::MinimizeLBFGS(P, opts.niter);
	}

	/*
	 * (5) save MRF
	 */
	if (opts.mrf != NULL && opts.rmode > 0) {
		MRF.Save(opts.mrf);
	}

	/*
	 * (6) do MRF post-processing
	 */
	switch (opts.rmode) {
	case 0:
		printf("# Contact matrix correction: ZILCH\n");
		break;
	case 1:
		printf("# Contact matrix correction: FN\n");
		MRFprocessor::FN(MRF, mtx);
		Contacts.AddFeature("FN", mtx);
		break;
	case 2:
		printf("# Contact matrix correction: APC\n");
		MRFprocessor::APC(MRF, mtx);
		Contacts.AddFeature("APC", mtx);
		break;
	case 3:
		printf("# Contact matrix correction: PROB\n");
		MRFprocessor::APC(MRF, mtx);
		MRFprocessor::Zscore(ncol, mtx);
		Contacts.AddFeature("Z(APC)", mtx);
		if (ncol < 150) {
			printf("# Coefficients set: 1\n");
			Contacts.RescoreLogistic(coef_range1);
		} else if (ncol >= 150 && ncol < 300) {
			printf("# Coefficients set: 2\n");
			Contacts.RescoreLogistic(coef_range2);
		} else {
			printf("# Coefficients set: 3\n");
			Contacts.RescoreLogistic(coef_range3);
		}
		MRFprocessor::CH(MRF, mtx);
		Contacts.AddFeature("Cor(h)", mtx);
		break;
	default:
		printf("!!! ACHTUNG !!! (this should never happen)\n");
		return 1;
	}

	/*
	 * (7) save MTX
	 */
	if (opts.mtx != NULL) {
		Contacts.SaveMTX(opts.mtx, MSA);
	}

	/*
	 * (7) calculate per group scores
	 */
//	if (opts.umask != NULL) {
//		for (EDGE_LIST::iterator it = L.begin(); it != L.end(); it++) {
//			std::vector<std::pair<size_t, size_t> > vec = MSA.CastToMsa(
//					it->second);
//			double score = MRFprocessor::GetScore(result, vec);
//			double energy = MRF.GetPairEnergies(MSA, vec);
//			printf("# Score(%s)= %.6e, Energy(%s)= %.6e\n", it->first.c_str(),
//					score, it->first.c_str(), energy);
//		}
//	}
//	}
	/*
	 *
	 */
	Contacts.Print(MSA);

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
	while ((tmp = getopt(argc, argv, "hi:o:f:n:r:c:m:u:R:p:s:")) != -1) {
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
		case 'R': /* regularization mode */
			if (strcmp(optarg, "FN") == 0) {
				opts.rmode = 1;
			} else if (strcmp(optarg, "APC") == 0) {
				opts.rmode = 2;
			} else if (strcmp(optarg, "PROB") == 0) {
				opts.rmode = 3;
			} else if (strcmp(optarg, "ZILCH") == 0) {
				opts.rmode = 0;
			} else {
				printf("Error: wrong matrix correction mode '%s'\n", optarg);
				return false;
			}
			break;
		case 'p': /* lpair - pair penalty (symmetric) */
			opts.lpair = atof(optarg);
			if (opts.lpair <= 0.0) {
				printf("Error: 'lpair' should be positive (-p)\n");
				return false;
			}
			break;
		case 's': /* lskew - asymmetric penalty for couplings */
			opts.lskew = atof(optarg);
			if (opts.lskew <= 0.0) {
				printf("Error: 'lskew' should be positive (-s)\n");
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
	printf("          -p lpair\n");
	printf("          -s lskew\n");

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

	free(msa);

	for (size_t p = 0; p < ncol; p++) {
		for (size_t q = p + 1; q < ncol; q++) {
			mtx[p][q] /= nrow;
			mtx[q][p] = mtx[p][q];
		}
	}

	return 2.0 * E / (ncol - 1) / ncol;

}

double PairEnergiesDI(const MSAclass &MSA, double **mtx) {

	double DI = 0.0;

	RRCE RRCE_(RRCE::RRCE20RC, 7.8, 5);

	size_t nrow = MSA.GetNrow();
	size_t ncol = MSA.GetNcol();

	unsigned char * msa = MSA.GetMsa();
	MSAclass::aatoi(msa, nrow * ncol);

	double Neff = MSA.GetNeff();

	for (size_t i = 0; i < ncol; i++) {
		memset(mtx[i], 0, ncol * sizeof(double));
	}

	/*
	 * precompute 1-site probabilities
	 */
	double *Px = (double*) malloc(20 * sizeof(double));
	double Z = 0.0;
	for (int i = 0; i < 20; i++) {
		Px[i] = exp(-RRCE_.GetHi(i));
		Z += Px[i];
	}
	for (int i = 0; i < 20; i++) {
		Px[i] /= Z;
	}

	/*
	 * precomputed 2-site delta-log-likelihoods
	 */
	double **Pxy = (double**) malloc(20 * sizeof(double));
	Z = 0.0;
	for (int i = 0; i < 20; i++) {
		Pxy[i] = (double*) malloc(20 * sizeof(double));
		for (int j = 0; j < 20; j++) {
			Pxy[i][j] = exp(
					-RRCE_.GetHi(i) - RRCE_.GetHi(j) - RRCE_.GetJij(i, j));
			Z += Pxy[i][j];
		}
	}

	for (int i = 0; i < 20; i++) {
		for (int j = 0; j < 20; j++) {
			Pxy[i][j] /= Z;
		}
	}

	/*
	 * evaluate MSA
	 */
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

				mtx[p][q] += w * Pxy[a][b] * log(Pxy[a][b] / Px[a] / Px[b]);

			}
		}
	}

	/* symmetrize */
	for (size_t p = 0; p < ncol; p++) {
		for (size_t q = p + 1; q < ncol; q++) {
			mtx[p][q] /= Neff;
			mtx[q][p] = mtx[p][q];
		}
	}

	/*
	 * free
	 */
	free(msa);
	for (int i = 0; i < 20; i++) {
		free(Pxy[i]);
	}
	free(Pxy);
	free(Px);

	return DI;

}
