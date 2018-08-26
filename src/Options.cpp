/*
 * Options.cpp
 *
 *  Created on: Aug 21, 2018
 *      Author: ivan
 */

#include <string>
#include <cstring>

#include "Options.h"

bool GetOpts(int argc, char *argv[], OPTS &opts) {

	char tmp;
	while ((tmp = getopt(argc, argv, "hi:o:f:n:r:c:R:p:s:t:b:")) != -1) {
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
		case 'R': /* regularization mode */
			if (strcmp(optarg, "FN") == 0) {
				opts.rmode = 1;
			} else if (strcmp(optarg, "APC") == 0) {
				opts.rmode = 2;
			} else if (strcmp(optarg, "PROB5") == 0) {
				opts.rmode = 3;
			} else if (strcmp(optarg, "PROB8") == 0) {
				opts.rmode = 4;
			} else if (strcmp(optarg, "ZILCH") == 0) {
				opts.rmode = 0;
			} else {
				printf("Error: wrong matrix correction mode '%s'\n", optarg);
				return false;
			}
			break;
		case 't': /* number of threads to use */
			opts.nthreads = atoi(optarg);
			break;
		case 'b': /* APC-corrected contact map (for bbconacts) */
			opts.apc = optarg;
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

	printf("\nUsage:   ./gremlin3 [-option] [argument]\n\n");
	printf("Options:  -i alignment.a3m               - input, required\n");
	printf("          -o matrix.txt                  - output, optional\n");
	printf("          -b apcmatrix.txt               - output, optional\n");
//	printf("          -f mrf.txt                     - output, optional\n");
	printf("          -n number of iterations          %ld\n", opts.niter);
	printf("          -r max gaps per row [0;1)        %.2lf\n", opts.grow);
	printf("          -c max gaps per column [0;1)     %.2lf\n", opts.gcol);
	printf("          -R contact matrix correction\n");
	printf("             {FN,APC,PROB5,PROB8}          PROB8\n");
	printf("          -t number of threads             %d\n", opts.nthreads);

}

void PrintCap(const OPTS &opts) {

	time_t timer;
	time(&timer);
	struct tm* tm_info = localtime(&timer);
	char buf[100];
	strftime(buf, 26, "%Y:%m:%d / %H:%M:%S", tm_info);

	printf("# %s\n", std::string(70, '-').c_str());
	printf("# gremlin - a program to predict protein contact maps %18s\n",
	VERSION);
	printf("#           from multiple sequence alignments\n");
	printf("# %s\n", std::string(70, '-').c_str());

	printf("# %20s : %s\n", "start date/time", buf);
	printf("# %20s : %s\n", "MSA file", opts.a3m);
	if (opts.mtx != NULL) {
		printf("# %20s : %s\n", "contact matrix", opts.mtx);
	}
	if (opts.mrf != NULL) {
		printf("# %20s : %s\n", "MRF", opts.mrf);
	}
	printf("# %20s : %d\n", "threads", opts.nthreads);

	printf("# %s\n", std::string(70, '-').c_str());

//	printf("# %10s %15s %10s %10s %10s %10s %10s %10s %5s %5s %5s\n", "TMPLT",
//			"best_params", "cont_sco", "gap_sco", "max_scoA", "max_scoB",
//			"tot_scoA", "tot_scoB", "Nali", "lenA", "lenB");
//	printf("#\n");

}

