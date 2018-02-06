/*
 * map_align.cpp
 *
 *  Created on: Jan 26, 2018
 *      Author: ivan
 */

#include <unistd.h>
#include <string>
//#include <thread>
#include <future>

#include <omp.h>

#include "MSAclass.h"
#include "Chain.h"
#include "Info.h"
#include "CMap.h"
#include "MapAlign.h"

#define DMAX 5.0
#define KMIN 3

struct OPTS {
	std::string seq; /* sequence file */
	std::string con; /* contacts file */
	std::string pdb; /* PDB file */
	std::string out; /* output file */
	std::string dir; /* folder where templates are stored */
	std::string list; /* list of template IDs (file) */
	std::string prefix; /* prefix to output matches */
	int num; /* number of models to save */
	int nthreads; /* number of threads to use */
};

bool GetOpts(int argc, char *argv[], OPTS &opts);
void PrintOpts(const OPTS &opts);

CMap MapFromPDB(const Chain &C);
void SaveMatch(std::string, const Chain&, const std::vector<int>&,
		const std::string&);
void SaveAtom(FILE *F, Atom *A, int atomNum, int resNum, char type);

std::pair<std::string, MP_RESULT> Align(const CMap&, const OPTS&,
		const MapAlign::PARAMS& params, const std::string&);

int main(int argc, char *argv[]) {

	/*
	 * (0) process input parameters
	 */
	OPTS opts = { "", "", "", "", "", "", "", 0, 1 };
	if (!GetOpts(argc, argv, opts)) {
		PrintOpts(opts);
		return 1;
	}
	omp_set_num_threads(opts.nthreads);
	printf("# %20s : %d\n", "number of threads", opts.nthreads);
	printf("# %20s : %s\n", "sequence file", opts.seq.c_str());
	printf("# %20s : %s\n", "contacts file", opts.con.c_str());
	MapAlign::PARAMS params = { -1.0, -0.01, 3, 10 };

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
	printf("# %20s : %.3f\n", "max query score", MapAlign::MaxScore(mapA));

	/*
	 * (2) process single PDB input file (if any)
	 */
	if (opts.pdb != "") {
		Chain C(opts.pdb);
		CMap mapB(MapFromPDB(C));
		std::vector<int> a2b;
		MP_RESULT result = MapAlign::Align(mapA, mapB, params);
		if (opts.out != "") {
			SaveMatch(opts.out, C, a2b, seqA);
		}

		return 0;

	}

	/*
	 * (3) process template library
	 */

	/* read IDs */
	std::vector<std::string> listB;
	{
		FILE *F = fopen(opts.list.c_str(), "r");
		if (F == NULL) {
			printf("Error: cannot open list file '%s'\n", opts.list.c_str());
			return 1;
		}
		char id[1024];
		while (fscanf(F, "%s\n", id) == 1) {
			listB.push_back(id);
		}
		fclose(F);
	}
	printf("# %20s : %lu\n", "IDs read", listB.size());

	/* TODO: make separate loops for single CPU
	 *       and multithreaded executions */

	/* TODO: make output of Align(...) function
	 *       more meaningful */

	/* read PDBs one by one and calculate alignments */
#pragma omp parallel for
	for (unsigned i = 0; i < listB.size(); i++) {
		std::pair<std::string, MP_RESULT> result = Align(mapA, opts, params,
				listB[i]);
#pragma omp critical
		{
			printf("--> %s %s", result.first.c_str(),
					result.second.label.c_str());
			for (auto &s : result.second.sco) {
				printf(" %.3f", s);
			}
			for (auto &s : result.second.len) {
				printf(" %d", s);
			}
			printf("\n");
		}

	}

	/*
	 * (4) save top hits
	 */

	return 0;

}

void PrintOpts(const OPTS &opts) {

	printf("\nUsage:   ./map_align [-option] [argument]\n\n");
	printf("Options:  -s sequence.fas (input, required)\n");
	printf("          -c contacts.txt (input, required)\n\n");
	printf("          ************* single template ************\n");
	printf("          -p template.pdb (input)\n");
	printf("          -o match.pdb (output)\n\n");
	printf("          ********** library of templates **********\n");
	printf("          -D PATH_TO_TEMPLATES (input)\n");
	printf("          -L list.txt (input)\n");
	printf("          -O PREFIX for saving top hits (input)\n");
	printf("          -N number of top hits to save (input)\n");
	printf("          -M max number of residues in the template\n\n");
	printf("          ****************** misc ******************\n");
	printf("          -t number of threads\n\n");

}

bool GetOpts(int argc, char *argv[], OPTS &opts) {

	char tmp;
	while ((tmp = getopt(argc, argv, "hs:c:p:o:D:L:O:N:v:t:")) != -1) {
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
		case 'D': /* path to templates */
			opts.dir = std::string(optarg);
			break;
		case 'L': /* list file of template IDs */
			opts.list = std::string(optarg);
			break;
		case 'O': /* prefix to save top hits */
			opts.prefix = std::string(optarg);
			break;
		case 'N': /* list file of template IDs */
			opts.num = atoi(optarg);
			break;
		case 't': /* number of threads */
			opts.nthreads = atoi(optarg);
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

	bool fl_both = (opts.pdb != "" && (opts.dir != "" || opts.list != ""));
	bool fl_none = (opts.pdb == "" && (opts.dir == "" || opts.list == ""));

	if (fl_both || fl_none) {
		printf("Error: set either a PDB file ('-p') or "
				"a library of templates ('-D' and '-L')\n");
		return false;
	}

	return true;

}

CMap MapFromPDB(const Chain &C) {

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
				adj[i].push_back(std::make_tuple(j, 1.0, sep));
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

void SaveAtom(FILE *F, Atom *A, int atomNum, int resNum, char type) {

	const char *resName = AAA3[MSAclass::aatoi(type)];

	fprintf(F,
			"ATOM  %5d  %-3s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
			atomNum, A->name, A->altLoc, resName, 'A', resNum, ' ', A->x, A->y,
			A->z, A->occup, A->temp, A->element, A->charge);

}

void SaveMatch(std::string name, const Chain& C, const std::vector<int>& a2b,
		const std::string& seq) {

	FILE *F = fopen(name.c_str(), "w");
	if (F == NULL) {
		printf("Error: cannot open PDB file for saving '%s'\n", name.c_str());
		return;
	}

	for (unsigned i = 0; i < a2b.size(); i++) {
		int idx = a2b[i];
		if (idx > -1) {
			Residue &R = C.residue[idx];
			SaveAtom(F, R.N, i * 4 + 1, i + 1, seq[i]);
			SaveAtom(F, R.CA, i * 4 + 2, i + 1, seq[i]);
			SaveAtom(F, R.C, i * 4 + 3, i + 1, seq[i]);
			SaveAtom(F, R.O, i * 4 + 4, i + 1, seq[i]);
		}
	}

	fclose(F);

}

std::pair<std::string, MP_RESULT> Align(const CMap& mapA, const OPTS& opts,
		const MapAlign::PARAMS& params, const std::string& id) {

	std::string name = opts.dir + "/" + id + ".pdb";
	Chain B(name);
	MP_RESULT result;
	if (B.nRes > 20 && B.nRes < 1000) {
		CMap mapB = MapFromPDB(B);
		result = MapAlign::Align(mapA, mapB, params);
		return std::make_pair(id, result);
	} else {
		return std::make_pair(id, result);
	}

}
