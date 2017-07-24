#include <cstdio>

#include "RRCE.h"
#include "EigenRRCE.h"
#include "MSAclass.h"
#include "ProblemFull.h"
#include "Minimizer.h"
#include "MRFprocessor.h"

/* TODO: program parameters
 * 1) a3m file (in)
 * 2) mtx file (out)
 * 3) gaps per row (def = 0.25)
 * 4) gaps per col (def = 0.25)
 * 5) number of iterations (def = 25)
 * 6) constraints: mask (-m) OR unmask (-u) edges
 * 7) mrf file (out) */

int main(int argc, char *argv[]) {

	MSAclass MSA(argv[1]);

	MSA.CleanMsa(0.25, 0.25);

	printf("# %ld x %ld\n", MSA.GetNcol(), MSA.GetNrow());

	ProblemFull P(MSA);

	MRFclass MRF = Minimizer::Minimize(P, 100);
	MRF.Save("mrf.txt");

	size_t dim = MRF.GetDim();
	double **mtx = (double**) malloc(dim * sizeof(double*));
	for (size_t i = 0; i < dim; i++) {
		mtx[i] = (double*) malloc(dim * sizeof(double));
	}

	MRFprocessor::APC(MRF, mtx);

	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			printf("%.5e ", mtx[i][j]);
		}
		printf("\n");
	}

	for (size_t i = 0; i < dim; i++) {
		free(mtx[i]);
	}
	free(mtx);

	return 0;

}
