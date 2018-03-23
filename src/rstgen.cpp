/*
 * rstgen.cpp
 *
 *  Created on: Mar 23, 2018
 *      Author: aivan
 */

/*
 * Input:
 * 		-i alignment.a3m (take sequence 1 only)
 * 		-m file.mtx
 * 		-f fraction (f*Len)
 * 		-m mode {SIG,BND}
 * 		-p probability cutoff (for BND)
 * 		-k sequence separation
 *
 * Output:
 * 		-o rosetta.rst
 */

#include <string>

#include "MSAclass.h"

struct OPTS {

	std::string a3m; /* A3M file */
	std::string mtx; /* file with the computed contact matrix */
	std::string rst; /* Rosetta restraints file */

	double f; /* fraction of top pairs */
	double p; /* probability cutoff */

	int kmin; /* sequence separation */

	int type; /* SIG, BND */

};

bool GetOpts(int argc, char *argv[], OPTS &opts);
void PrintOpts(const OPTS &opts);
void PrintCap(const OPTS &opts);

int main(int argc, char *argv[]) {

	return 0;

}
