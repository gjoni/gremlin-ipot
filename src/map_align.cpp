/*
 * map_align.cpp
 *
 *  Created on: Jan 26, 2018
 *      Author: ivan
 */

#include <unistd.h>
#include <string>
#include <algorithm>

#include <omp.h>

#include "MSAclass.h"
#include "Chain.h"
#include "Info.h"
#include "CMap.h"
#include "MapAlign.h"

#define DMAX 5.0
#define KMIN 3
#define VERSION "V20180206"

struct OPTS {
	std::string seq; /* sequence file */
	std::string con; /* contacts file */
	std::string pdb; /* PDB file */
	std::string out; /* output file */
	std::string dir; /* folder where templates are stored */
	std::string list; /* list of template IDs (file) */
	std::string prefix; /* prefix to output matches */
	unsigned num; /* number of models to save */
	int nthreads; /* number of threads to use */
};

bool GetOpts(int argc, char *argv[], OPTS &opts);
void PrintOpts(const OPTS &opts);

void PrintCap(const OPTS &opts);

CMap MapFromPDB(const Chain &C);
void SaveMatch(std::string, const Chain&, const std::vector<int>&,
		const std::string&);
void SaveAtom(FILE *F, Atom *A, int atomNum, int resNum, char type);

MP_RESULT Align(const CMap&, const OPTS&, const MapAlign::PARAMS& params,
		const std::string&);

bool compare(const std::pair<std::string, MP_RESULT> &a,
		const std::pair<std::string, MP_RESULT> &b) {
	return (a.second.sco[0] > b.second.sco[0]);
}

int main(int argc, char *argv[]) {

	/*
	 * (0) process input parameters
	 */
	OPTS opts = { "", "", "", "", "", "", "", 0, 1 };
	if (!GetOpts(argc, argv, opts)) {
		PrintOpts(opts);
		return 1;
	}

#if defined(_OPENMP)
	omp_set_num_threads(opts.nthreads);
#endif

	MapAlign::PARAMS params = { -1.0, -0.01, 3, 10 };
	PrintCap(opts);

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

	/* read PDBs one by one and calculate alignments */
	std::vector<std::pair<std::string, MP_RESULT> > hits;

#if defined(_OPENMP)
#pragma omp parallel for
#endif
	for (unsigned i = 0; i < listB.size(); i++) {
		MP_RESULT result = Align(mapA, opts, params, listB[i]);

#if defined(_OPENMP)
#pragma omp critical
#endif
		{
			if (result.sco.size()) {
				printf("# %10s %15s", listB[i].c_str(), result.label.c_str());
				for (auto &s : result.sco) {
					printf(" %10.3f", s);
				}
				for (auto &s : result.len) {
					printf(" %5d", s);
				}
				printf("\n");
				hits.push_back(std::make_pair(listB[i], result));
			} else {
				printf("# %10s %15s\n", listB[i].c_str(), "...skipped...");
			}
		}

	}

	/*
	 * (4) save top hits
	 */
	if (opts.num) {

		printf("# %s\n", std::string(70, '-').c_str());

		/* make sure that topN results exist */
		opts.num = hits.size() < opts.num ? hits.size() : opts.num;

		/* sort in decreasing order of
		 * contact_score + gap_score */
		std::sort(hits.begin(), hits.end(), compare);

		/* print info about topN matches */
		for (unsigned i = 0; i < opts.num; i++) {
			std::string &id = hits[i].first;
			MP_RESULT &result = hits[i].second;
			printf("T %10s %15s", id.c_str(), result.label.c_str());
			for (auto &s : result.sco) {
				printf(" %10.3f", s);
			}
			for (auto &s : result.len) {
				printf(" %5d", s);
			}
			printf("\n");

			/* save partial matches (if requested) */
			if (opts.prefix != "") {
				std::string name = opts.dir + "/" + id + ".pdb";
				Chain B(name);
				name = opts.prefix + id + ".pdb";
				SaveMatch(name, B, result.a2b, seqA);
			}
		}
	}

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

void PrintCap(const OPTS &opts) {

	printf("# %s\n", std::string(70, '-').c_str());
	printf("# map_align - a program to align protein contact maps %18s\n",
	VERSION);
	printf("# %s\n", std::string(70, '-').c_str());

	printf("# %20s : %s\n", "sequence file", opts.seq.c_str());
	printf("# %20s : %s\n", "contacts file", opts.con.c_str());
	printf("# %20s : %s\n", "list file", opts.list.c_str());
	printf("# %20s : %s\n", "path to templates", opts.dir.c_str());
	printf("# %20s : %d\n", "threads", opts.nthreads);

	printf("# %s\n", std::string(70, '-').c_str());

	printf("#\n");
	printf("# %10s %15s %10s %10s %10s %10s %10s %10s %5s %5s %5s\n", "TMPLT",
			"best_params", "cont_sco", "gap_sco", "max_scoA", "max_scoB",
			"tot_scoA", "tot_scoB", "Nali", "lenA", "lenB");
	printf("#\n");

}

bool GetOpts(int argc, char *argv[], OPTS &opts) {

	char tmp;
	while ((tmp = getopt(argc, argv, "hs:c:p:o:D:L:O:N:v:t:M:")) != -1) {
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

MP_RESULT Align(const CMap& mapA, const OPTS& opts,
		const MapAlign::PARAMS& params, const std::string& id) {

	std::string name = opts.dir + "/" + id + ".pdb";
	Chain B(name);
	MP_RESULT result;
	if (B.nRes > 20 && B.nRes < 1000) {
		CMap mapB = MapFromPDB(B);
		result = MapAlign::Align(mapA, mapB, params);
	}

	return result;

}
