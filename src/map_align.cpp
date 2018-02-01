/*
 * map_align.cpp
 *
 *  Created on: Jan 26, 2018
 *      Author: ivan
 */

#include <unistd.h>
#include <string>

#include "CMap.h"
#include "MSAclass.h"
#include "Chain.h"
#include "MapAlign.h"

#define DMAX 5.0
#define KMIN 3

struct OPTS {
	std::string seq; /* sequence file */
	std::string con; /* contacts file */
	std::string pdb; /* PDB file */
	std::string out; /* output file */
};

bool GetOpts(int argc, char *argv[], OPTS &opts);
void PrintOpts(const OPTS &opts);
CMap MapFromPDB(const std::string& name);

int main(int argc, char *argv[]) {

	/*
	 * (0) process input parameters
	 */
	OPTS opts = { "", "", "", "" };
	if (!GetOpts(argc, argv, opts)) {
		PrintOpts(opts);
		return 1;
	}

	/*
	 * (1) create contact map object
	 *     from contacts file and sequence
	 */
	std::string seqA;
	{
		MSAclass msa(opts.seq.c_str());
		seqA = msa.GetSequence(0);
	}
	CMap mapA(opts.con, seqA);
//	mapA.Print();

	/*
	 * (2) create contact map object from PDB
	 */
	CMap mapB = MapFromPDB(opts.pdb);
//	mapB.Print();

	/*
	 * (3) align two maps
	 */
	std::vector<int> a2b(mapA.Size(), -1);
	MapAlign::PARAMS params = { -1.0, -0.01, 3, 20 };
	MapAlign::Align(mapA, mapB, params, a2b);

	return 0;

}

void PrintOpts(const OPTS &opts) {

	printf("Usage:   ./map_align [-option] [argument]\n");
	printf("Options:  -s sequence.fas (input, required)\n");
	printf("          -c contacts.txt (input, required)\n");
	printf("          -p template.pdb (input, required)\n");
	printf("          -o match.pdb (output)\n");

}

bool GetOpts(int argc, char *argv[], OPTS &opts) {

	char tmp;
	while ((tmp = getopt(argc, argv, "hs:c:p:o:")) != -1) {
		switch (tmp) {
		case 'h': /* help */
			printf("!!! HELP !!!\n");
			return false;
			break;
		case 's': /* sequence file */
			opts.seq = std::string(optarg);
			break;
		case 'c': /* contacts file */
			opts.con = std::string(optarg);
			break;
		case 'p': /* PDB file */
			opts.pdb = std::string(optarg);
			break;
		case 'o': /* match file */
			opts.out = std::string(optarg);
			break;
		default:
			return false;
			break;
		}
	}

	if (opts.seq == "") {
		printf("Error: sequence file not specified ('-s')\n");
		return false;
	}

	if (opts.con == "") {
		printf("Error: contacts file not specified ('-c')\n");
		return false;
	}

	if (opts.pdb == "") {
		printf("Error: PDB file not specified ('-p')\n");
		return false;
	}

	return true;

}

CMap MapFromPDB(const std::string& name) {

	Chain C(name.c_str());

	std::string seq(C.nRes, 'X');
	for (int i = 0; i < C.nRes; i++) {
		seq[i] = MSAclass::itoaa(C.residue[i].type);
	}

	/* temp. adjacency matrix */
	double **mtx = (double**) malloc(C.nRes * sizeof(double*));
	for (int i = 0; i < C.nRes; i++) {
		mtx[i] = (double*) calloc(C.nRes, sizeof(double));
	}

	/* find neighbors */
	kdres *res;
	double pos[3];
	for (int i = 0; i < C.nAtoms; i++) {

		Atom *A = C.atom[i];
		int a = A->residue - C.residue;

		if (A->type == 'H') { /* exclude hydrogens */
			continue;
		}

		res = kd_nearest_range3f(C.kd, A->x, A->y, A->z, DMAX);
		while (!kd_res_end(res)) {
			Atom *B = *((Atom**) kd_res_item(res, pos));
			int b = B->residue - C.residue;
			mtx[a][b] = 1.0;
			mtx[b][a] = 1.0;
			kd_res_next(res);
		}
		kd_res_free(res);

	}

	/* create CMap object */
	AListT adj(C.nRes);
	for (int i = 0; i < C.nRes; i++) {
		for (int j = 0; j < C.nRes; j++) {
			int sep = abs(i - j);
			if (sep < KMIN) {
				continue;
			}
			if (mtx[i][j] > 1.0e-10) {
				adj[i].push_back( { j, 1.0, sep });
			}
		}
	}
	CMap map(adj, seq);

	/* free */
	for (int i = 0; i < C.nRes; i++) {
		free(mtx[i]);
	}
	free(mtx);

	return map;

}
