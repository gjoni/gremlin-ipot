#include <cstdio>
#include <unistd.h>

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

/* TODO: separate programs
 * 1) gremlin3 - produces an MRF
 * 2) score_patchdock
 * 3) score_gramm
 * 4) score_decoy */

struct OPTS {
	char *a3m; /* A3M file */
	char *mtx; /* file with the computed contact matrix */
	char *mrf; /* file to save MRF */
	size_t niter; /* number of iterations */
	double grow; /* gaps per row */
	double gcol; /* gaps per col */
	char *mask;
	char *umask;
};

int main(int argc, char *argv[]) {

	MSAclass MSA(argv[1]);

	MSA.CleanMsa(0.25, 0.25);

	printf("# %ld x %ld\n", MSA.GetNcol(), MSA.GetNrow());

	ProblemFull P(MSA);

	MRFclass MRF = Minimizer::Minimize(P, 25);
	MRF.Save("mrf.txt");

	MRFprocessor::MTX result;
	MRFprocessor::APC(MRF, result);
	MRFprocessor::SaveMTX(result, "mtx.txt");

	return 0;

}
