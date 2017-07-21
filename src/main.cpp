#include <cstdio>

#include "RRCE.h"
#include "EigenRRCE.h"
#include "MSAclass.h"
#include "ProblemFull.h"
#include "Minimizer.h"

int main(int argc, char *argv[]) {

	MSAclass MSA(argv[1]);

	MSA.CleanMsa(0.25, 0.25);
//	MSA.SaveMSA("SPC19.msa");
//	MSA.CastToIdx();

	printf("# %ld x %ld\n", MSA.GetNcol(), MSA.GetNrow());

	/*
	 for (size_t i = 0; i < nrow; i++) {
	 for (size_t j = 0; j < ncol; j++) {
	 printf(" %d", msa[i * ncol + j]);
	 }
	 printf("\n");
	 }
	 */

	ProblemFull P(MSA);

	MRFclass MRF = Minimizer::Minimize(P);

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
