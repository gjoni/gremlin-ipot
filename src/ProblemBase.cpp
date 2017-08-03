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
		MSA(NULL), msa(NULL), dim(0), w(NULL), we(NULL) {

	/* */

}

ProblemBase::ProblemBase(const MSAclass &MSA_) :
		MSA(&MSA_), msa(NULL), dim(MSA_.GetNcol() * MSA_.GetNrow()), w(NULL), we(
		NULL) {

	AllocateBase();

	msa = MSA->GetMsa();

	Reweight();

	MSAclass::aatoi(msa, dim);

	UnmaskAllEdges();

}

ProblemBase::ProblemBase(const ProblemBase &source) :
		MSA(source.MSA), msa(NULL), dim(source.MSA->ncol * source.MSA->nrow), w(
		NULL), we(
		NULL) {

	AllocateBase();

	size_t msa_dim = MSA->nrow * MSA->ncol;

	msa = (unsigned char *) malloc(msa_dim * sizeof(unsigned char));

	memcpy(msa, source.msa, msa_dim * sizeof(unsigned char));
	memcpy(w, source.w, MSA->nrow * sizeof(double));
	memcpy(we, source.we, MSA->ncol * MSA->ncol * sizeof(bool));

}

ProblemBase::~ProblemBase() {

	FreeBase();

}

void ProblemBase::Reweight(double t) {

	size_t ncol = MSA->ncol;
	size_t nrow = MSA->nrow;

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

void ProblemBase::MaskEdges(const std::vector<std::pair<size_t, size_t> > &e) {

	/* set all we[][] to 1.0 */
	UnmaskAllEdges();

	/* mask the specified edges */
	for (auto const& edge : e) {
		size_t i = MSA->GetMsaIdx(edge.first);
		size_t j = MSA->GetMsaIdx(edge.second);
		if (i < SIZE_MAX && j < SIZE_MAX) {
			we[i * MSA->ncol + j] = false;
			we[j * MSA->ncol + i] = false;
		}

	}

}

void ProblemBase::UnmaskEdges(
		const std::vector<std::pair<size_t, size_t> > &e) {

	/* set all we[][] to false */
	memset(we, 0, MSA->ncol * MSA->ncol * sizeof(bool));

	/* mask the specified edges */
	for (auto const& edge : e) {
		size_t i = MSA->GetMsaIdx(edge.first);
		size_t j = MSA->GetMsaIdx(edge.second);
		if (i < SIZE_MAX && j < SIZE_MAX) {
			we[i * MSA->ncol + j] = true;
			we[j * MSA->ncol + i] = true;
		}
	}

}

void ProblemBase::UnmaskAllEdges() {

	/* set all we[][] to true */
	for (size_t i = 0; i < MSA->ncol * MSA->ncol; i++) {
		we[i] = true;
	}

}

void ProblemBase::AllocateBase() {

	w = (double*) malloc(MSA->nrow * sizeof(double));
	we = (bool*) malloc(MSA->ncol * MSA->ncol * sizeof(bool));

}

void ProblemBase::FreeBase() {

	free(w);
	free(we);
	free(msa);

}

size_t ProblemBase::GetDim() {

	return dim;

}
