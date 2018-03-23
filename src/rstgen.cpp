/*
 * rstgen.cpp
 *
 *  Created on: Mar 23, 2018
 *      Author: aivan
 */

#include <cstring>
#include <cassert>
#include <unistd.h>

#include <string>
#include <map>
#include <tuple>
#include <algorithm>

#include "MSAclass.h"

struct OPTS {

	std::string a3m; /* A3M file */
	std::string mtx; /* file with the computed contact matrix */
	std::string rst; /* Rosetta restraints file */

	std::string type; /* SIG, BND */

	double f; /* fraction of top pairs */
	double p; /* probability cutoff */

	int kmin; /* sequence separation */

	std::string seq;

};

bool GetOpts(int argc, char *argv[], OPTS &opts);
void PrintOpts(const OPTS &opts);

using namespace std;

typedef map<string, tuple<string, string, double, double> > MAP;

MAP ReadTable();

struct Contact {
	int i, j;
	char a, b;
	double s;
};

vector<Contact> ReadContacts(const OPTS&);

int main(int argc, char *argv[]) {

	/*
	 * (0) process input parameters
	 */
	OPTS opts = { "", "", "", "SIG", 1.5, 0.95, 3 };
	if (!GetOpts(argc, argv, opts)) {
		PrintOpts(opts);
		return 1;
	}

	/*
	 * (1) read MSA
	 */
	MSAclass MSA(opts.a3m);
	opts.seq = MSA.GetSequence(0);

	/*
	 * (2) read parameters table
	 */
	MAP TBL = ReadTable();

	/*
	 * (3) read contacts & sort
	 */
	vector<Contact> cont = ReadContacts(opts);
	sort(cont.begin(), cont.end(),
			[](const Contact &c1, const Contact &c2)->bool {
				return c1.s > c2.s;
			});

	/*
	 * (4) output
	 */
	FILE *F = fopen(opts.rst.c_str(), "w");
	if (F == NULL) {
		printf("Error: cannot open restraints file '%s'\n", opts.rst.c_str());
		exit(1);
	}

	unsigned N = opts.seq.length() * opts.f;
	unsigned counter = 0;
	for (auto &c : cont) {
		string key = { c.a, c.b };
		MAP::iterator it = TBL.find(key);
		if (it == TBL.end()) {
			printf("# WARNING: key '%s' not found\n", key.c_str());
			continue;
		}
		if (opts.type == "SIG") {
			fprintf(F, "AtomPair %s %d %s %d SCALARWEIGHTEDFUNC %.6f "
					"SUMFUNC 2 SIGMOID %.3f %.3f CONSTANTFUNC -0.5\n",
					get<0>(it->second).c_str(), c.i, get<1>(it->second).c_str(),
					c.j, c.s, get<2>(it->second), get<3>(it->second));
		} else if (opts.type == "BND" && c.s > opts.p) {
			fprintf(F, "AtomPair %s %d %s %d SCALARWEIGHTEDFUNC %.6f "
					"BOUNDED 0 %.3f 1 0.5 #\n", get<0>(it->second).c_str(), c.i,
					get<1>(it->second).c_str(), c.j, c.s,
					get<2>(it->second) + 1.0 / get<3>(it->second));
		}
		counter++;
		if (counter > N) {
			break;
		}
	}
	fclose(F);

	return 0;

}

bool GetOpts(int argc, char *argv[], OPTS &opts) {

	char tmp;
	while ((tmp = getopt(argc, argv, "hi:m:f:t:p:k:o:")) != -1) {
		switch (tmp) {
		case 'h': /* help */
			printf("!!! HELP !!!\n");
			return false;
			break;
		case 'i': /* A3M file (in) */
			opts.a3m = std::string(optarg);
			break;
		case 'm': /* MTX file (out) */
			opts.mtx = std::string(optarg);
			break;
		case 'o': /* Rosetta restraints file */
			opts.rst = std::string(optarg);
			break;
		case 'f': /* fraction of top residues */
			opts.f = atof(optarg);
			break;
		case 'p': /* probability cutoff */
			opts.p = atof(optarg);
			break;
		case 'k': /* sequence separation */
			opts.kmin = atof(optarg);
			break;
		case 't': /* restraints type */
			opts.type = std::string(optarg);
			break;
		default:
			return false;
			break;
		}
	}

	if (opts.a3m == "") {
		printf("Error: A3M file not specified '-i'\n");
		return false;
	}

	if (opts.mtx == "") {
		printf("Error: MTX file not specified '-m'\n");
		return false;
	}

	if (opts.rst == "") {
		printf("Error: restraints file not specified '-o'\n");
		return false;
	}

	if (opts.type != "SIG" && opts.type != "BND") {
		printf("Error: either sigmoid (SIG) or bounded (BND) restraints "
				"are supported '-t'\n");
		return false;
	}

	return true;

}

void PrintOpts(const OPTS &opts) {

	printf("\nUsage:   ./rstgen [-option] [argument]\n\n");
	printf("Options:  -i alignment.a3m               - input, required\n");
	printf("          -m matrix.txt                  - input, required\n");
	printf("          -o restraints.txt              - output, required\n");
	printf("          -t restraint type {SIG, BND}     %s\n",
			opts.type.c_str());
	printf("          -f fraction of top contacts      %.2f * Len\n", opts.f);
	printf("          -p probability cutoff            %.2f\n", opts.p);
	printf("          -k sequence separation           %d\n", opts.kmin);
	printf("\n");

}

MAP ReadTable() {

	char *DATADIR = getenv("GREMLINDAT");
	if (DATADIR == NULL) {
		printf("Error: environment variable 'GREMLINDAT' not set\n");
		exit(1);
	}

	char name[1024];
	sprintf(name, "%s/cb-cb.txt", DATADIR);

	MAP MAP_;

	FILE *F = fopen(name, "r");
	char s1[5], s2[5], s3[5];
	double p1, p2;
	while (fscanf(F, "%s %s %s %lf %lf\n", s1, s2, s3, &p1, &p2) == 5) {
		MAP_[string(s1)] = make_tuple(s2, s3, p1, p2);
	}
	fclose(F);

	return MAP_;

}

vector<Contact> ReadContacts(const OPTS &opts) {

	vector<Contact> vec;
	FILE *F = fopen(opts.mtx.c_str(), "r");
	if (F == NULL) {
		printf("Error: cannot open matrix file '%s'\n", opts.mtx.c_str());
		exit(1);
	}

	const unsigned SIZE = 65536;
	char buf[SIZE];

	int dim = 0;

	Contact c;
	while (fgets(buf, SIZE, F) != NULL) {
		dim++;
	}
	rewind(F);

	assert(dim == (int )opts.seq.length()); /* length mismatch */

	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			fscanf(F, "%lf", &(c.s));
			if (j - i >= opts.kmin && c.s > 1.0e-6) {
				c.i = i + 1;
				c.j = j + 1;
				c.a = opts.seq[i];
				c.b = opts.seq[j];
				vec.push_back(c);
			}
		}
	}

	fclose(F);

	return vec;

}
