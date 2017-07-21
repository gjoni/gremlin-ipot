/*
 * ProblemBase.cpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ivan
 */

#include <cstring>
#include <cmath>

#include "ProblemBase.h"

ProblemBase::ProblemBase() :
		MSA(NULL), dim(0), w(NULL), we(NULL) {

	/* */

}

ProblemBase::ProblemBase(MSAclass &MSA_) :
		MSA(&MSA_), dim(0), w(NULL), we(NULL) {

	AllocateBase();

	Reweight();

	UnmaskAllEdges();

}

ProblemBase::~ProblemBase() {

	FreeBase();

}

void ProblemBase::Reweight(double t) {

	size_t ncol = MSA->ncol;
	size_t nrow = MSA->nrow;

	unsigned char const *msa = MSA->msa;

	size_t idthres = (size_t) ceil(t * ncol);

	memset(w, 0, sizeof(double) * nrow);

	size_t nij = nrow * (nrow + 1) / 2;

	for (size_t ij = 0; ij < nij; ij++) {

		// compute i and j from ij
		// http://stackoverflow.com/a/244550/1181102
		size_t i, j;
		{
			size_t ii = nrow * (nrow + 1) / 2 - 1 - ij;
			size_t K = floor((sqrt(8 * ii + 1) - 1) / 2);
			i = nrow - 1 - K;
			j = ij - nrow * i + i * (i + 1) / 2;
		}

		size_t ids = 0;
		for (size_t k = 0; k < ncol; k++) {
			if (msa[i * ncol + k] == msa[j * ncol + k]) {
				ids++;
			}
		}

		if (ids > idthres) {
			w[i]++;
			w[j]++;
		}
	}

	for (size_t i = 0; i < nrow; i++) {
		w[i] = 1. / (w[i] - 1);
	}

	double wsum = 0;
	double wmin = w[0], wmax = w[0];
	for (size_t i = 0; i < nrow; i++) {
		double wt = w[i];
		wsum += wt;
		if (wt > wmax) {
			wmax = wt;
		}
		if (wt < wmin) {
			wmin = wt;
		}
	}

	printf("# threshold= %.1f Beff= %g mean= %g min= %g max= %g\n", t, wsum,
			wsum / nrow, wmin, wmax);

}

void ProblemBase::MaskEdges(const std::vector<std::pair<int, int> > &e) {

	/* set all we[][] to 1.0 */
	UnmaskAllEdges();

	/* mask the specified edges */
	for (auto const& edge : e) {
		size_t i = MSA->GetMsaIdx(edge.first);
		size_t j = MSA->GetMsaIdx(edge.second);
		if (i < SIZE_MAX && j < SIZE_MAX) {
			we[i][j] = we[j][i] = 0.0;
		}

	}

}

void ProblemBase::UnmaskEdges(const std::vector<std::pair<int, int> > &e) {

	/* set all we[][] to 0.0 */
	for (size_t i = 0; i < MSA->ncol; i++) {
		for (size_t j = 0; j < MSA->ncol; j++) {
			we[i][j] = 0.0;
		}
	}

	/* mask the specified edges */
	for (auto const& edge : e) {
		size_t i = MSA->GetMsaIdx(edge.first);
		size_t j = MSA->GetMsaIdx(edge.second);
		if (i < SIZE_MAX && j < SIZE_MAX) {
			we[i][j] = we[j][i] = 1.0;
		}

	}

}

void ProblemBase::UnmaskAllEdges() {

	/* set all we[][] to 1.0 */
	for (size_t i = 0; i < MSA->ncol; i++) {
		for (size_t j = 0; j < MSA->ncol; j++) {
			we[i][j] = 1.0;
		}
	}

}

void ProblemBase::AllocateBase() {

	w = (double*) malloc(MSA->nrow * sizeof(double));
	we = (double**) malloc(MSA->ncol * sizeof(double*));
	for (size_t i = 0; i < MSA->ncol; i++) {
		we[i] = (double*) malloc(MSA->ncol * sizeof(double));
	}

}

void ProblemBase::FreeBase() {

	if (we != NULL) {
		for (size_t i = 0; i < MSA->ncol; i++) {
			free(we[i]);
		}
		free(we);
	}
	free(w);

	w = NULL;
	we = NULL;

}

size_t ProblemBase::GetDim() {

	return dim;

}
