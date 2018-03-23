/*
 * neff.cpp
 *
 *  Created on: Mar 23, 2018
 *      Author: aivan
 */

#include <cmath>
#include <omp.h>

#include "MSAclass.h"

int main(int argc, char *argv[]) {

	if (argc != 5) {
		printf("\nUsage:   ./neff alignment.a3m ROW_GAPS COL_GAPS NCPU\n\n");
		return 0;
	}

	double grow = atof(argv[2]);
	double gcol = atof(argv[3]);
	int nthreads = atoi(argv[4]);

#if defined(_OPENMP)
	omp_set_num_threads(nthreads);
#endif

	MSAclass MSA(argv[1]);
	MSA.CleanMsa(grow, gcol);
	MSA.Reweight();

	double Neff = MSA.GetNeff();

	printf("Neff= %.2lf   Nf= %.2f\n", Neff, Neff / sqrt(MSA.GetNcol()));

	return 0;

}
