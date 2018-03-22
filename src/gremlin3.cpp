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
#include "ProblemFullOMP.h"
#include "Minimizer.h"
#include "MRFprocessor.h"
#include "ContactList.h"
#include "RRCE.h"
#include "LogRegCoeff.h"

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
};

bool GetOpts(int argc, char *argv[], OPTS &opts);
void PrintOpts(const OPTS &opts);

double PairEnergies(const MSAclass &MSA, double **mtx);
double PairEnergiesDI(const MSAclass &MSA, double **mtx);

int main(int argc, char *argv[]) {

	/*
	 * (0) process input parameters
	 */
	OPTS opts = { NULL, NULL, NULL, 50, 0.25, 0.25, NULL, NULL, 4 };
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

	/* temp storage for various scores */
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
	MSA.HxHy(mtx);
	Contacts.AddFeature("Hx+Hy", mtx);

	/* gaps */
//	MSA.Gxy(mtx);
//	Contacts.AddFeature("Gxy", mtx);
//
//	MSA.GxGy(mtx);
//	Contacts.AddFeature("Gx+Gy", mtx);
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
//	ProblemFull P(MSA);
	ProblemFullOMP P(MSA);

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
		printf("# Contact matrix correction: PROB5\n");
		MRFprocessor::APC(MRF, mtx);
		MRFprocessor::Zscore(ncol, mtx);
		Contacts.AddFeature("Z(APC)", mtx);
		if (ncol < 100) {
			printf("#          Coefficients set: 0..100\n");
			Contacts.RescoreLogistic(coef5A_range1);
		} else if (ncol >= 100 && ncol < 150) {
			printf("#          Coefficients set: 100..150\n");
			Contacts.RescoreLogistic(coef5A_range2);
		} else if (ncol >= 150 && ncol < 200) {
			printf("#          Coefficients set: 150..200\n");
			Contacts.RescoreLogistic(coef5A_range3);
		} else if (ncol >= 200 && ncol < 250) {
			printf("#          Coefficients set: 200..250\n");
			Contacts.RescoreLogistic(coef5A_range4);
		} else if (ncol >= 250 && ncol < 300) {
			printf("#          Coefficients set: 250..300\n");
			Contacts.RescoreLogistic(coef5A_range5);
		} else if (ncol >= 300 && ncol < 400) {
			printf("#          Coefficients set: 300..400\n");
			Contacts.RescoreLogistic(coef5A_range6);
		} else {
			printf("#          Coefficients set: 400..inf\n");
			Contacts.RescoreLogistic(coef5A_range7);
		}
		break;
	case 4:
		printf("# Contact matrix correction: PROB8\n");
		MRFprocessor::APC(MRF, mtx);
		MRFprocessor::Zscore(ncol, mtx);
		Contacts.AddFeature("Z(APC)", mtx);
		if (ncol < 100) {
			printf("#          Coefficients set: 0..100\n");
			Contacts.RescoreLogistic(coef8A_range1);
		} else if (ncol >= 100 && ncol < 150) {
			printf("#          Coefficients set: 100..150\n");
			Contacts.RescoreLogistic(coef8A_range2);
		} else if (ncol >= 150 && ncol < 200) {
			printf("#          Coefficients set: 150..200\n");
			Contacts.RescoreLogistic(coef8A_range3);
		} else if (ncol >= 200 && ncol < 250) {
			printf("#          Coefficients set: 200..250\n");
			Contacts.RescoreLogistic(coef8A_range4);
		} else if (ncol >= 250 && ncol < 300) {
			printf("#          Coefficients set: 250..300\n");
			Contacts.RescoreLogistic(coef8A_range5);
		} else if (ncol >= 300 && ncol < 400) {
			printf("#          Coefficients set: 300..400\n");
			Contacts.RescoreLogistic(coef8A_range6);
		} else {
			printf("#          Coefficients set: 400..inf\n");
			Contacts.RescoreLogistic(coef8A_range7);
		}
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
			} else if (strcmp(optarg, "PROB5") == 0) {
				opts.rmode = 3;
			} else if (strcmp(optarg, "PROB8") == 0) {
				opts.rmode = 4;
			} else if (strcmp(optarg, "ZILCH") == 0) {
				opts.rmode = 0;
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
//	printf("          -m list1.txt - residue pairs to be masked\n");
//	printf("          -u list2.txt - residue pairs to be unmasked "
//			"(all others are masked)\n");
	printf(
			"          -R contact matrix correction {FN,APC,PROB5,PROB8} (PROB8)\n");

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
