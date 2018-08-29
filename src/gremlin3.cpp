#include <cstdio>
#include <ctime>
#include <cstring>
#include <cmath>

#include <string>

#include <omp.h>

#include "MSAclass.h"
#include "MRFclass.h"

#include "ProblemL2.h"
#include "ProblemL2_1b.h"

#include "Minimizer.h"
#include "ContactList.h"
#include "RRCE.h"
#include "LogRegCoeff.h"
#include "Options.h"

/* TODO: move this function somewhere else */
double PairEnergies(const MSAclass &MSA, double **mtx);

int main(int argc, char *argv[]) {

	/*
	 * (0) process input parameters
	 */
	OPTS opts = { NULL, NULL, NULL, 25, 0.25, 0.25, 4, 1, NULL };
	if (!GetOpts(argc, argv, opts)) {
		PrintOpts(opts);
		return 1;
	}

	PrintCap(opts);

#if defined(_OPENMP)
	omp_set_num_threads(opts.nthreads);
#endif

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

	if (opts.rmode == 5 || opts.rmode == 6) {
		if (opts.rmode == 6) {
			MRFclass::APC(ncol, mtx);
			Contacts.AddFeature("MIAPC", mtx);
		}
		if (opts.mtx != NULL) {
			Contacts.SaveMTX(opts.mtx, MSA);
		}
		return 0;
	}

	/* joint entropy */
	MSA.HxHy(mtx);
	Contacts.AddFeature("Hx+Hy", mtx);

	/* statistical potential */
	printf("# E(RRCE)= %.5f\n", PairEnergies(MSA, mtx));
	Contacts.AddFeature("RRCE", mtx);

	/*
	 * (3) set up the problem
	 */
	MRFclass MRF(MSA);
	if (opts.rmode > 0) {

		/* train 1-body */
		printf("# %s\n", std::string(70, '-').c_str());
		printf("# step 1: solve for local fields\n");
		printf("# %s\n", std::string(70, '-').c_str());
		ProblemL2_1b P1(MSA);
		Minimizer::MinimizeLBFGS(P1, opts.niter, MRF);

		/* train 1-body & 2-body */
		printf("# %s\n", std::string(70, '-').c_str());
		printf("# step 2: solve for local fields and couplings\n");
		printf("# %s\n", std::string(70, '-').c_str());
		ProblemL2 P(MSA);
		Minimizer::MinimizeLBFGS(P, opts.niter, MRF);

	}

	/*
	 * (5) save MRF
	 */
	if (opts.mrf != NULL && opts.rmode > 0) {
		MRF.Save(opts.mrf);
	}

	/* save APC-corrected contact map
	 * (for bbcontacts) */
	if (opts.apc != NULL) {
		ContactList ContactsTmp(ncol);
		MRF.APC(mtx);
		ContactsTmp.AddFeature("APC", mtx);
		ContactsTmp.SaveMTX(opts.apc, MSA);
	}

	/*
	 * (6) do MRF post-processing
	 */
	switch (opts.rmode) {
	case 0:
		printf("# Contact matrix correction: ZILCH\n");
		break;
	case 1:
		printf("# Contact matrix correction : FN\n");
		MRF.FN(mtx);
		Contacts.AddFeature("FN", mtx);
		break;
	case 2:
		printf("# Contact matrix correction : APC\n");
		MRF.APC(mtx);
		Contacts.AddFeature("APC", mtx);
		break;
	case 3:
		printf("# Contact matrix correction : PROB5\n");
		MRF.APC(mtx);
		MRF.Zscore(ncol, mtx);
		Contacts.AddFeature("Z(APC)", mtx);
		if (ncol < 100) {
			printf("#          Coefficients set : 0..100\n");
			Contacts.RescoreLogistic(coef5A_range1);
		} else if (ncol >= 100 && ncol < 150) {
			printf("#          Coefficients set : 100..150\n");
			Contacts.RescoreLogistic(coef5A_range2);
		} else if (ncol >= 150 && ncol < 200) {
			printf("#          Coefficients set : 150..200\n");
			Contacts.RescoreLogistic(coef5A_range3);
		} else if (ncol >= 200 && ncol < 250) {
			printf("#          Coefficients set : 200..250\n");
			Contacts.RescoreLogistic(coef5A_range4);
		} else if (ncol >= 250 && ncol < 300) {
			printf("#          Coefficients set : 250..300\n");
			Contacts.RescoreLogistic(coef5A_range5);
		} else if (ncol >= 300 && ncol < 400) {
			printf("#          Coefficients set : 300..400\n");
			Contacts.RescoreLogistic(coef5A_range6);
		} else {
			printf("#          Coefficients set : 400..inf\n");
			Contacts.RescoreLogistic(coef5A_range7);
		}
		break;
	case 4:
		printf("# Contact matrix correction : PROB8\n");
		MRF.APC(mtx);
		MRF.Zscore(ncol, mtx);
		Contacts.AddFeature("Z(APC)", mtx);
		if (ncol < 100) {
			printf("#          Coefficients set : 0..100\n");
			Contacts.RescoreLogistic(coef8A_range1);
		} else if (ncol >= 100 && ncol < 150) {
			printf("#          Coefficients set : 100..150\n");
			Contacts.RescoreLogistic(coef8A_range2);
		} else if (ncol >= 150 && ncol < 200) {
			printf("#          Coefficients set : 150..200\n");
			Contacts.RescoreLogistic(coef8A_range3);
		} else if (ncol >= 200 && ncol < 250) {
			printf("#          Coefficients set : 200..250\n");
			Contacts.RescoreLogistic(coef8A_range4);
		} else if (ncol >= 250 && ncol < 300) {
			printf("#          Coefficients set : 250..300\n");
			Contacts.RescoreLogistic(coef8A_range5);
		} else if (ncol >= 300 && ncol < 400) {
			printf("#          Coefficients set : 300..400\n");
			Contacts.RescoreLogistic(coef8A_range6);
		} else {
			printf("#          Coefficients set : 400..inf\n");
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

//	Contacts.Print(MSA);

	/*
	 * (8) finish date/time
	 */
	time_t timer;
	time(&timer);
	struct tm* tm_info = localtime(&timer);
	char buf[100];
	strftime(buf, 26, "%Y:%m:%d / %H:%M:%S", tm_info);
	printf("# %s\n", std::string(70, '-').c_str());
	printf("# %25s : %s\n", "end date/time", buf);
	printf("# %s\n", std::string(70, '-').c_str());

	/*
	 * free
	 */
	for (size_t i = 0; i < ncol; i++) {
		free(mtx[i]);
	}
	free(mtx);

	return 0;

}

double PairEnergies(const MSAclass &MSA, double **mtx) {

	double E = 0.0;

	RRCE RRCE_;

	size_t nrow = MSA.GetNrow();
	size_t ncol = MSA.GetNcol();

	unsigned char * msa = MSA.GetMsa();
	MSAclass::aatoi(msa, nrow * ncol);

	for (size_t i = 0; i < ncol; i++) {
		memset(mtx[i], 0, ncol * sizeof(double));
	}

#if defined(_OPENMP)
#pragma omp parallel for reduction (+:E)
#endif
	for (size_t i = 0; i < nrow; i++) {

		unsigned char *seq = msa + i * ncol;
		double w = MSA.GetWeight(i);

		double e = 0.0;
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
				e += j * w;

			}

		}

		E += e;

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

