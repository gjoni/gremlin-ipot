#include <cstdio>

#include "RRCE.h"
#include "EigenRRCE.h"
#include "MSAclass.h"
#include "ProblemFull.h"
#include "Minimizer.h"

int main(int argc, char *argv[]) {

	MSAclass MSA(argv[1]);
	MSA.CleanMsa(0.1, 0.1);

	MSAclass MSA2(MSA);

	MSA2.CastToIdx();

	unsigned char const *msa = MSA2.GetMsa();
	size_t nrow = MSA2.GetNrow();
	size_t ncol = MSA2.GetNcol();

	printf("# %ld x %ld\n", nrow, ncol);

	for (size_t i = 0; i < nrow; i++) {
		for (size_t j = 0; j < ncol; j++) {
			printf(" %d", msa[i * ncol + j]);
		}
		printf("\n");
	}

	ProblemFull P(MSA2);

	size_t NAA = MSAclass::NAA;
	gsl_vector*x = gsl_vector_alloc(ncol * NAA * (1 + ncol * NAA));
	double f = P.f(x);
	printf("# f(x)= %f\n", f);

	printf("start\n");
	f = Minimizer::Minimize(P);
	printf("stop\n");

	/*
	 RRCE RRCE_(argv[1]);

	 const int N = 20;
	 double **J = (double**) malloc(N * sizeof(double*));
	 for (int i = 0; i < N; i++) {
	 J[i] = (double*) malloc(N * sizeof(double));
	 }

	 RRCE_.GetCouplings(J);
	 EigenRRCE Eigen(J);

	 printf("%f\n", Eigen.GetReconstructionCorrel(atof(argv[2])));

	 for (int i = 0; i < N; i++) {
	 free(J[i]);
	 }
	 free(J);
	 */

	return 0;

}
