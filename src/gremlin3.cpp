#include <cstdio>
#include <unistd.h>
#include <ctime>
#include <cstring>

#include "RRCE.h"
#include "EigenRRCE.h"
#include "MSAclass.h"
#include "ProblemFull.h"
#include "Minimizer.h"
#include "MRFprocessor.h"

/* TODO: separate programs
 * 1) gremlin3 - produces an MRF
 * 2) score_patchdock
 * 3) score_gramm
 * 4) score_decoy */

/* TODO: issues
 *   - multiple copies of h,J (huge memory overhead) ???
 *   - minimizer requires a lot of memory
 */

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

int main(int argc, char *argv[]) {

	OPTS opts = { NULL, NULL, NULL, 25, 0.25, 0.25, NULL, NULL, 1 };

	if (!GetOpts(argc, argv, opts)) {
		PrintOpts(opts);
		return 1;
	}

	MSAclass MSA(opts.a3m);

	MSA.CleanMsa(opts.grow, opts.gcol);

	ProblemFull P(MSA);

	/*
	 srand(time(NULL));
	 std::vector<std::pair<int, int> > e;
	 for (int i = 0; i < 500; i++) {
	 int a = rand() % MSA.GetLen();
	 int b = rand() % MSA.GetLen();
	 if (a > b) {
	 e.push_back(std::make_pair(b, a));
	 }
	 if (a < b) {
	 e.push_back(std::make_pair(a, b));
	 }
	 }
	 P.UnmaskEdges(e);
	 */

	MRFclass MRF = Minimizer::MinimizeLBFGS(P, opts.niter);

	if (opts.mrf != NULL) {
		MRF.Save(opts.mrf);
	}

	if (opts.mtx != NULL) {
		MRFprocessor::MTX result;
		switch (opts.rmode) {
		case 1:
			MRFprocessor::APC(MRF, result);
			break;
		case 2:
			MRFprocessor::FN(MRF, result);
			break;
		default:
			printf("!!! ACHTUNG !!! (this should never happen)\n");
			return 1;
		}
		MRFprocessor::SaveMTX(result, opts.mtx);
	}

	return 0;

}

bool GetOpts(int argc, char *argv[], OPTS &opts) {

	char tmp;
	while ((tmp = getopt(argc, argv, "hi:o:f:n:r:c:m:u:R:")) != -1) {
		switch (tmp) {
		/*option h show the help infomation*/
		case 'h':
			printf("!!! HELP !!!\n");
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
			if (strcmp(optarg, "APC") == 0) {
				opts.rmode = 1;
			} else if (strcmp(optarg, "FN") == 0) {
				opts.rmode = 2;
			} else {
				printf("Error: wrong matrix correction mode '%s'\n", optarg);
				return false;
			}
			opts.umask = optarg;
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
//	printf("          -u list2.txt - residue pairs to be unmasked\n");
	printf("          -R contact matrix correction {APC,FN,...} (APC)\n");

}
