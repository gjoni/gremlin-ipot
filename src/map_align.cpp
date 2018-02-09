/*
 * map_align.cpp
 *
 *  Created on: Jan 26, 2018
 *      Author: ivan
 */

#include <unistd.h>
#include <string>
#include <algorithm>
#include <ctime>

#include <omp.h>

#include "MSAclass.h"
#include "Chain.h"
#include "Info.h"
#include "CMap.h"
#include "MapAlign.h"
#include "TMalign.h"

#define DMAX 5.0
#define KMIN 3
#define VERSION "V20180209"

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
	double tmmax; /* TM-score cut-off for cleaning of top matches */
	int maxres; /* max template size */
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

double TMscore(const Chain&, const Chain&, const std::vector<int>&,
		const std::vector<int>&);

bool compare(const std::pair<std::string, MP_RESULT> &a,
		const std::pair<std::string, MP_RESULT> &b) {
	return (a.second.sco[0] + a.second.sco[1]
			> b.second.sco[0] + b.second.sco[1]);
}

int main(int argc, char *argv[]) {

	/*
	 * (0) process input parameters
	 */
	OPTS opts = { "", "", "", "", "", "", "", 10, 1, 0.8, 1000 };
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
	int nskipped = 0;

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic) num_threads(opts.nthreads)
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
				nskipped++;
			}
			fflush(stdout);
		}

	}
	printf("# %s\n", std::string(70, '-').c_str());
	printf("# %20s : %d\n", "skipped", nskipped);

	/*
	 * (4) process top hits
	 */

	/* make sure that topN results exist */
	opts.num = hits.size() < opts.num ? hits.size() : opts.num;

	/* sort in decreasing order of
	 * contact_score + gap_score */
	std::sort(hits.begin(), hits.end(), compare);

	/* store info about top hits in these vectors */
	std::vector<Chain> chains;
	std::vector<std::pair<std::string, MP_RESULT> > top_hits;

	/*
	 * (4a) clean based on TM-score
	 */

	/* clean top hits by TM-score only if requested by user
	 * (and cut-off is reasonable) */
	if (opts.tmmax < 1.0 && opts.tmmax > 0.1) {

		/* loop over all hits */
		int counter = 0;
		for (auto &result : hits) {

			/* load template */
			Chain B(opts.dir + "/" + result.first.c_str() + ".pdb");
			std::vector<int> &b2ref = result.second.a2b;

			/* check whether current template is similar
			 * to any of the already processed templates */
			bool fl = 1;
			for (unsigned i = 0; i < chains.size(); i++) {
				Chain &C = chains[i];
				std::vector<int> &c2ref = top_hits[i].second.a2b;
				double tm = TMscore(B, C, b2ref, c2ref);
				if (tm > opts.tmmax) {
					fl = 0;
					break;
				}
			}

			/* save this hit if no similarity to
			 * previous hits detected */
			if (fl) {
				chains.push_back(B);
				top_hits.push_back(result);
			} else {
				counter++;
			}

			/* stop if enough hits collected */
			if (chains.size() >= opts.num) {
				break;
			}
		}

		printf("# %20s : %d\n", "similarities", counter);
		printf("# %s\n", std::string(70, '-').c_str());

	} else {

		/* process topN hits without clustering */
		top_hits.assign(hits.begin(), hits.begin() + opts.num);
		for (auto &result : top_hits) {
			std::string &id = result.first;
			chains.push_back(Chain(opts.dir + "/" + id.c_str() + ".pdb"));
		}

	}

	/*
	 * (4b) print info about topN matches
	 */

	/* make sure that topN results exist after clustering */
	opts.num = top_hits.size() < opts.num ? top_hits.size() : opts.num;

	for (unsigned i = 0; i < opts.num; i++) {
		std::string &id = top_hits[i].first;
		MP_RESULT &result = top_hits[i].second;
		printf("T %10s %15s", id.c_str(), result.label.c_str());
		for (auto &s : result.sco) {
			printf(" %10.3f", s);
		}
		for (auto &l : result.len) {
			printf(" %5d", l);
		}
		printf("\n");

		/* save partial matches (if requested by user) */
		if (opts.prefix != "") {
			SaveMatch(opts.prefix + id + ".pdb", chains[i], result.a2b, seqA);
		}
	}

	/*
	 * (5) finish date/time
	 */
	time_t timer;
	time(&timer);
	struct tm* tm_info = localtime(&timer);
	char buf[100];
	strftime(buf, 26, "%Y:%m:%d / %H:%M:%S", tm_info);
	printf("# %s\n", std::string(70, '-').c_str());
	printf("# %20s : %s\n", "end date/time", buf);
	printf("# %s\n", std::string(70, '-').c_str());

	return 0;

}

void PrintOpts(const OPTS &opts) {

	printf("\nUsage:   ./map_align [-option] [argument]\n\n");
	printf("Options:  -s sequence.fas                - input, required\n");
	printf("          -c contacts.txt                - input, required\n\n");
	printf("          ***************** single template ****************\n");
	printf("          -p template.pdb                - input, required\n");
	printf("          -o match.pdb                   - output, optional\n\n");
	printf("                                  OR                        \n");
	printf("          ************** library of templates **************\n");
	printf("          -D path to templates           - input, required\n");
	printf("          -L list.txt with template IDs  - input, required\n");
	printf("          -O prefix for saving top hits  - output, optional\n");
	printf("          -N number of top hits to save    (%u)\n", opts.num);
	printf("          -T TM-score cleaning cut-off     (%.2f)\n", opts.tmmax);
	printf("          -M max template size             (%d)\n\n", opts.maxres);
	printf("          ********************** misc **********************\n");
	printf("          -t number of threads             (%d)\n", opts.nthreads);
	printf("\n");

}

void PrintCap(const OPTS &opts) {

	time_t timer;
	time(&timer);
	struct tm* tm_info = localtime(&timer);
	char buf[100];
	strftime(buf, 26, "%Y:%m:%d / %H:%M:%S", tm_info);

	printf("# %s\n", std::string(70, '-').c_str());
	printf("# map_align - a program to align protein contact maps %18s\n",
	VERSION);
	printf("# %s\n", std::string(70, '-').c_str());

	printf("# %20s : %s\n", "start date/time", buf);
	printf("# %20s : %s\n", "sequence file", opts.seq.c_str());
	printf("# %20s : %s\n", "contacts file", opts.con.c_str());
	printf("# %20s : %s\n", "list file", opts.list.c_str());
	printf("# %20s : %s\n", "path to templates", opts.dir.c_str());
	printf("# %20s : %d\n", "threads", opts.nthreads);
	printf("# %20s : %.3f\n", "TM-score cut-off", opts.tmmax);

	printf("# %s\n", std::string(70, '-').c_str());

	printf("#\n");
	printf("# %10s %15s %10s %10s %10s %10s %10s %10s %5s %5s %5s\n", "TMPLT",
			"best_params", "cont_sco", "gap_sco", "max_scoA", "max_scoB",
			"tot_scoA", "tot_scoB", "Nali", "lenA", "lenB");
	printf("#\n");

}

bool GetOpts(int argc, char *argv[], OPTS &opts) {

	char tmp;
	while ((tmp = getopt(argc, argv, "hs:c:p:o:D:L:O:N:v:t:M:T:")) != -1) {
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
		case 'T': /* TM-score clustering cut-off */
			opts.tmmax = atof(optarg);
			break;
		case 'M': /* max template size */
			opts.tmmax = atof(optarg);
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
			//int sep = abs(i - j);
			int sep = abs(C.residue[i].seqNum - C.residue[j].seqNum);
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
			if (R.N == NULL || R.C == NULL || R.O == NULL) {
				continue;
			}
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
	if (B.nRes > 20 && B.nRes < opts.maxres) {
		CMap mapB = MapFromPDB(B);
		result = MapAlign::Align(mapA, mapB, params);
	}

	return result;

}

double TMscore(const Chain& A, const Chain& B, const std::vector<int>& a2ref,
		const std::vector<int>& b2ref) {

	/* allocate memory */
	unsigned dim = a2ref.size();
	double **x = (double**) malloc(dim * sizeof(double*));
	double **y = (double**) malloc(dim * sizeof(double*));
	for (unsigned i = 0; i < dim; i++) {
		x[i] = (double*) malloc(3 * sizeof(double));
		y[i] = (double*) malloc(3 * sizeof(double));
	}

	/* get aligned residues */
	dim = 0;
	unsigned dimA = 0, dimB = 0;
	for (unsigned i = 0; i < a2ref.size(); i++) {

		int idxa = a2ref[i];
		int idxb = b2ref[i];

		if (idxa > -1 && idxb > -1) {

			x[dim][0] = A.residue[idxa].CA->x;
			x[dim][1] = A.residue[idxa].CA->y;
			x[dim][2] = A.residue[idxa].CA->z;

			y[dim][0] = B.residue[idxb].CA->x;
			y[dim][1] = B.residue[idxb].CA->y;
			y[dim][2] = B.residue[idxb].CA->z;

			dim++;

		}

		dimA += (idxa > -1);
		dimB += (idxb > -1);

	}

	double tm = 0.0;
	if (dim >= 5) {
		TMalign TM;
		tm = TM.GetTMscore(x, y, dim);
		tm = tm * dim / std::min(dimA, dimB);
	}

	/* free */
	for (unsigned i = 0; i < a2ref.size(); i++) {
		free(x[i]);
		free(y[i]);
	}
	free(x);
	free(y);

	return tm;

}
