#include <cstdio>

#include "RRCE.h"
#include "EigenRRCE.h"
#include "MSAclass.h"

int main(int argc, char *argv[]) {

	MSAclass MSA(argv[1]);
	MSA.CleanMsa(0.1,0.1);

	MSAclass MSA2(MSA);

	unsigned char const *msa = MSA2.GetMsa();
	size_t nrow = MSA2.GetNrow();
	size_t ncol = MSA2.GetNcol();

	printf("# %ld x %ld\n", nrow, ncol);

	for (size_t i = 0; i < nrow; i++) {
		for (size_t j = 0; j < ncol; j++) {
			printf("%c", msa[i * ncol + j]);
		}
		printf("\n");
	}

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
