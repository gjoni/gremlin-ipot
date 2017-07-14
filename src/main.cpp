
#include <cstdio>

#include "RRCE.h"
#include "EigenRRCE.h"

int main(int argc, char *argv[]) {

	RRCE RRCE_(argv[1]);

	const int N = 20;
	double **J = (double**) malloc(N * sizeof(double*));
	for (int i = 0; i < N; i++) {
		J[i] = (double*) malloc(N * sizeof(double));
	}

	RRCE_.GetCouplings(J);
	EigenRRCE Eigen(J);

	printf("%f\n", Eigen.GetReconstructionCorrel(atof(argv[2])));

	/*
	 * free
	 */
	for (int i = 0; i < N; i++) {
		free(J[i]);
	}
	free(J);

	return 0;

}