/*
 * Options.h
 *
 *  Created on: Aug 21, 2018
 *      Author: ivan
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <unistd.h>

#define VERSION "V20180829"

struct OPTS {
	char *a3m; /* A3M file */
	char *mtx; /* file with the computed contact matrix */
	char *mrf; /* file to save MRF */
	size_t niter; /* number of iterations */
	double grow; /* gaps per row */
	double gcol; /* gaps per col */
	int rmode; /* regularization mode */
	int nthreads; /* number of threads to use */
	char *apc; /* APC-corrected contact map */
};

bool GetOpts(int argc, char *argv[], OPTS &opts);
void PrintOpts(const OPTS &opts);
void PrintCap(const OPTS &opts);

#endif /* OPTIONS_H_ */
