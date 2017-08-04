#include <cstdio>
#include <unistd.h>
#include <ctime>
#include <cstring>

#include <string>
#include <map>
#include <vector>
#include <utility>

#include "MSAclass.h"
#include "ProblemFull.h"
#include "Minimizer.h"
#include "MRFprocessor.h"

/* TODO: separate programs
 * 1) gremlin3 - produces an MRF
 * 2) score_patchdock
 * 3) score_gramm
 * 4) score_decoy */

/* TODO:
 *   - GREMLIN2 symmetric minimizer
 *   - don't store masked edges (simplified problem with less variables)
 *   - GREMLIN1 routine with L1 penalty (ADMM solver needed)
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

typedef std::map<std::string, std::vector<std::pair<size_t, size_t> > > EDGE_LIST;
EDGE_LIST ReadEdges(const char *name);

bool GetOpts(int argc, char *argv[], OPTS &opts);
void PrintOpts(const OPTS &opts);

int main(int argc, char *argv[]) {

	/*
	 * (0) process input parameters
	 */
	OPTS opts = { NULL, NULL, NULL, 25, 0.25, 0.25, NULL, NULL, 1 };
	if (!GetOpts(argc, argv, opts)) {
		PrintOpts(opts);
		return 1;
	}

	/*
	 * (1) read & clean MSA
	 */
	MSAclass MSA(opts.a3m);
	MSA.CleanMsa(opts.grow, opts.gcol);

	/*
	 * (2) read edge constraints
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
	 * (3) prepare the problem
	 */
	ProblemFull P(MSA);
	if (opts.mask != NULL) {
		P.UnmaskAllEdges();
		for (EDGE_LIST::iterator it = L.begin(); it != L.end(); it++) {
			P.MaskEdges(it->second);
		}
	} else if (opts.umask != NULL) {
		for (EDGE_LIST::iterator it = L.begin(); it != L.end(); it++) {
			P.UnmaskEdges(it->second);
		}
	}

	/*
	 * (4) solve P
	 */
	MRFclass MRF = Minimizer::MinimizeLBFGS(P, opts.niter);

	/*
	 * (5) save results
	 */
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
	printf("          -m list1.txt - residue pairs to be masked\n");
	printf("          -u list2.txt - residue pairs to be unmasked\n");
	printf("          -R contact matrix correction {APC,FN,...} (APC)\n");

}

EDGE_LIST ReadEdges(const char *name) {

	EDGE_LIST L;

	FILE *F = fopen(name, "r");
	if (F == NULL) {
		printf("Error: cannot open file for reading '%s'\n", name);
		exit(1);
	}

	size_t a, b;
	const size_t SIZE = 256;
	char buf[SIZE];
	while (fscanf(F, "%lu %lu %s\n", &a, &b, buf) == 3) {
		L[buf].push_back(std::make_pair(a, b));
	}
	fclose(F);

	return L;

}
