/*
 * MapAlign.cpp
 *
 *  Created on: Jan 26, 2018
 *      Author: ivan
 */

#include <cstring>
#include <cmath>
#include <unordered_map>

#include "MapAlign.h"

MapAlign::MapAlign() {
	// TODO Auto-generated constructor stub

}

MapAlign::~MapAlign() {
	// TODO Auto-generated destructor stub
}

using namespace std;

void MapAlign::Alloc(SWDATA *swdata) {

	unsigned M = swdata->M;
	unsigned N = swdata->N;

	swdata->mtx = (double**) malloc(M * sizeof(double*));
	swdata->sco = (double**) malloc((M + 1) * sizeof(double*));
	swdata->label = (char**) malloc((M + 1) * sizeof(char*));
	for (unsigned i = 0; i < M; i++) {
		swdata->mtx[i] = (double*) malloc(N * sizeof(double));
		swdata->sco[i] = (double*) malloc((N + 1) * sizeof(double));
		swdata->label[i] = (char*) malloc((N + 1) * sizeof(char));
	}
	swdata->sco[M] = (double*) malloc((N + 1) * sizeof(double));
	swdata->label[M] = (char*) malloc((N + 1) * sizeof(char));

}

void MapAlign::Free(SWDATA* swdata) {

	for (unsigned i = 0; i < swdata->M; i++) {
		free(swdata->mtx[i]);
		free(swdata->sco[i]);
		free(swdata->label[i]);
	}
	free(swdata->sco[swdata->M]);
	free(swdata->label[swdata->M]);

	free(swdata->mtx);
	free(swdata->sco);
	free(swdata->label);

}

void MapAlign::InitMTX(SWDATA& swdata, const CMap& A, const CMap& B,
		double sep_x, double sep_y) {

	/* init mtx[][] */
	for (unsigned i = 0; i < swdata.M; i++) {
		memset(swdata.mtx[i], 0, swdata.N * sizeof(double));
	}

	/* left of diagonal */
	for (auto &idxa : A.GetLeftMap()) {
		const NListT& listA = A.GetLeftList(idxa);
		for (auto &idxb : B.GetLeftMap()) {
			const NListT& listB = B.GetLeftList(idxb);
			swdata.mtx[idxa][idxb] += SW1(listA, listB, sep_x, sep_y);
		}
	}

	/* right of diagonal */
	for (auto &idxa : A.GetRightMap()) {
		const NListT& listA = A.GetRightList(idxa);
		for (auto &idxb : B.GetRightMap()) {
			const NListT& listB = B.GetRightList(idxb);
			swdata.mtx[idxa][idxb] += SW1(listA, listB, sep_x, sep_y);
		}
	}

}

double MapAlign::Intersect(const NListT& listA, const NListT& listB,
		const std::vector<int>& a2b) {

	double score = 0.0;

	for (auto &a : listA) {
		int mapb = a2b[get<0>(a)];
		if (mapb < 0) {
			continue;
		}
		double sco_a = get<1>(a);
		unsigned sep_a = get<2>(a);
		for (auto &b : listB) {
			if ((unsigned) mapb == get<0>(b)) {
				score += sco_a * get<1>(b) * sepw(min(sep_a, get<2>(b)));
			}
		}
	}

	return score;

}

void MapAlign::UpdateMTX(SWDATA& swdata, const CMap& A, const CMap& B,
		double gap_e, int iter, std::vector<int>& a2b) {

	/* temp mtx[][] matrix */
	double **mtx = (double**) malloc(swdata.M * sizeof(double*));
	for (unsigned i = 0; i < swdata.M; i++) {
		mtx[i] = (double*) malloc(swdata.N * sizeof(double));
		memcpy(mtx[i], swdata.mtx[i], swdata.M);
	}

	/* iterate iter times */
	for (int it = 0; it < iter; it++) {

		/* align */
		SW2(swdata, gap_e, a2b);

		/*
		 * update similarity matrix
		 */
		double IT = 1.0 + it;
		double s1 = IT / (IT + 1.0);
		double s2 = 1.0 - s1;

		/* left of diagonal */
		for (auto &idxa : A.GetLeftMap()) {
			const NListT& listA = A.GetLeftList(idxa);
			for (auto &idxb : B.GetLeftMap()) {
				const NListT& listB = B.GetLeftList(idxb);
				double s = Intersect(listA, listB, a2b);
				swdata.mtx[idxa][idxb] += s * s2 / s1;
			}
		}

		/* right of diagonal */
		for (auto &idxa : A.GetRightMap()) {
			const NListT& listA = A.GetRightList(idxa);
			for (auto &idxb : B.GetRightMap()) {
				const NListT& listB = B.GetRightList(idxb);
				double s = Intersect(listA, listB, a2b);
				swdata.mtx[idxa][idxb] += s * s2 / s1;
			}
		}

		/* rescale scoring matrix */
		for (unsigned i = 0; i < swdata.M; i++) {
			for (unsigned j = 0; j < swdata.N; j++) {
				if (swdata.mtx[i][j] > 1.0e-10) {
					swdata.mtx[i][j] *= s1;
				}
			}
		}

	}

	/* copy temp sco[][] back & free */
	for (unsigned i = 0; i < swdata.M; i++) {
		memcpy(swdata.mtx[i], mtx[i], swdata.M);
		free(mtx[i]);
	}
	free(mtx);

}
vector<double> MapAlign::Assess(const CMap& A, const CMap& B,
		const std::vector<int>& a2b, double gap_e_w) {

	vector<double> scores;

	/* score matched contacts */
	double sco_con = 0.0;
	for (auto &c : A.edges) {
		int i = a2b[c.first.first];
		int j = a2b[c.first.second];
		if (i < 0 || j < 0) {
			continue;
		}
		EListT::const_iterator it = B.edges.find( { i, j });
		if (it != B.edges.end()) {
			sco_con += c.second.first * it->second.first
					* sepw(min(c.second.second, it->second.second));
		}
	}
	scores.push_back(sco_con);

	/* gap penalty score */

	return scores;

}

double MapAlign::Align(const CMap& A, const CMap& B, PARAMS& par,
		vector<int>& a2b) {

	double score = 0.0;

	/*
	 * (1) init alignment workspace
	 */
	SWDATA swdata = { A.Size(), B.Size(), NULL, NULL, NULL, vector<double>(
			A.Size(), par.gap_open), vector<double>(B.Size(), par.gap_open),
			vector<int>(A.Size()), vector<int>(B.Size()) };
	Alloc(&swdata);

	/*
	 * (2) alignment routine
	 */

	/* try different sep (sequence separation difference) penalties */
	for (auto &sep_x : vector<double> { 0, 1, 2 }) {

		/* try different scaling factors for sep penalties */
		for (auto &sep_y : vector<double> { 1, 2, 4, 8, 16, 32 }) {

			/* get initial score matrix */
			InitMTX(swdata, A, B, sep_x, sep_y);

			/* try different gap_ext penalties */
			for (auto &gap_e : vector<double> { 0.2, 0.1, 0.01, 0.001 }) {

				vector<int> a2b(swdata.M, -1);
				UpdateMTX(swdata, A, B, gap_e, par.iter, a2b);
				vector<double> scores = Assess(A, B, a2b, 0.0);

				printf("PAR: %.1e %.1e %.1e %12.3f\n", sep_x, sep_y, gap_e,
						scores[0]);

			}

		}

	}

	/*
	 * (10) free
	 */
	Free(&swdata);

	return score;

}

double MapAlign::gaussian(double mean, double std, double x) {

	double f = (x - mean) / std;
	return exp(-0.5 * f * f);

}

double MapAlign::sepw(double sep) {

	if (sep <= 4) {
		return 0.50;
	} else if (sep == 5) {
		return 0.75;
	} else {
		return 1.00;
	}

}

double MapAlign::SW1(const NListT& A, const NListT& B, double sep_x,
		double sep_y) {

	double score = 0.0;

	unsigned rows = A.size();
	unsigned cols = B.size();

	/* init DP matrix */
	double **sco = (double**) malloc((rows + 1) * sizeof(double*));
	for (unsigned i = 0; i < (rows + 1); i++) {
		sco[i] = (double*) calloc(cols + 1, sizeof(double));
	}

	/* run DP forward step
	 * keeping track of the best score */
	for (unsigned i = 1; i <= rows; i++) {

		/* probability and separation for contact A[i] */
		double sco_a = get<1>(A[i - 1]);
		unsigned sep_a = get<2>(A[i - 1]);

		for (unsigned j = 1; j <= cols; j++) {

			/* probability and separation for contact B[j] */
			double sco_b = get<1>(B[j - 1]);
			unsigned sep_b = get<2>(B[j - 1]);

			/* score of the match between
			 * contacts A[i] and B[j] */
			double sep_D = fabs((double) sep_a - (double) sep_b);
			double sep_M = min((double) sep_a, (double) sep_b);
			double sep_std = sep_y * (1.0 + pow(sep_M - 2.0, sep_x));

			double s = 0.0;
			if (sep_D / sep_std < 6) {
				s = sco_a * sco_b * sepw(sep_M) * gaussian(0, sep_std, sep_D);
			}

			/* update DP matrix */
			double A = sco[i - 1][j - 1] + s;
			double D = sco[i - 1][j];
			double R = sco[i][j - 1];

			if (A >= R) {
				if (A >= D) {
					sco[i][j] = A;
				} else {
					sco[i][j] = D;
				}
			} else {
				if (R >= D) {
					sco[i][j] = R;
				} else {
					sco[i][j] = D;
				}
			}

			/* update best score */
			score = sco[i][j] > score ? sco[i][j] : score;

		}

	}

	/* free */
	for (unsigned i = 0; i < (rows + 1); i++) {
		free(sco[i]);
	}
	free(sco);

	return score;

}

double MapAlign::SW2(SWDATA& swdata, double gap_e, std::vector<int>& a2b) {

	double score = 0.0;

	unsigned rows = swdata.M;
	unsigned cols = swdata.N;

	a2b.assign(rows, -1);

	/* init DP space */
	memset(swdata.sco[0], 0, (cols + 1) * sizeof(double));
	memset(swdata.label[0], 0, (cols + 1) * sizeof(char));
	for (unsigned i = 0; i < rows + 1; i++) {
		swdata.sco[i][0] = 0.0;
		swdata.label[i][0] = 0;
	}

	/* DP forward step */
	unsigned max_i = 0, max_j = 0;
	for (unsigned i = 1; i <= rows; i++) {

		for (unsigned j = 1; j <= cols; j++) {
			double A = swdata.sco[i - 1][j - 1] + swdata.mtx[i - 1][j - 1];
			double D = swdata.sco[i - 1][j];
			double R = swdata.sco[i][j - 1];

			if (swdata.label[i - 1][j] == 1) {
				D += swdata.gap_b[j - 1];
			} else {
				D += swdata.gap_b[j - 1] * gap_e;
			}
			if (swdata.label[i][j - 1] == 1) {
				R += swdata.gap_a[i - 1];
			} else {
				R += swdata.gap_a[i - 1] * gap_e;
			}

			if (A <= 0 && D <= 0 && R <= 0) {
				swdata.label[i][j] = 0;
				swdata.sco[i][j] = 0;
			} else {
				if (A >= R) {
					if (A >= D) {
						swdata.label[i][j] = 1;
						swdata.sco[i][j] = A;
					} else {
						swdata.label[i][j] = 2;
						swdata.sco[i][j] = D;
					}
				} else {
					if (R >= D) {
						swdata.label[i][j] = 3;
						swdata.sco[i][j] = R;
					} else {
						swdata.label[i][j] = 2;
						swdata.sco[i][j] = D;
					}
				}
				if (swdata.sco[i][j] > score) {
					max_i = i;
					max_j = j;
					score = swdata.sco[i][j];
				}
			}
		}
	}

	/* DP backward step to get the mapping */
	unsigned i = max_i, j = max_j;
	while (1) {
		if (swdata.label[i][j] == 0) {
			break;
		} else if (swdata.label[i][j] == 1) {
			a2b[i - 1] = j - 1;
			i--;
			j--;
		} else if (swdata.label[i][j] == 2) {
			i--;
		} else if (swdata.label[i][j] == 3) {
			j--;
		}
	}

	return score;

}

