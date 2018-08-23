/*
 * Problem1Body.cpp
 *
 *  Created on: Aug 23, 2018
 *      Author: aivan
 */

#include <cstring>
#include <cmath>
#include <cassert>

#include <omp.h>

#include "Problem1Body.h"

Problem1Body::Problem1Body() :
		ProblemBase(), lsingle(0.0), dim1body(0), ea(NULL), pa(NULL), lpa(NULL) {

	/* nothing to be done */

}

Problem1Body::Problem1Body(const MSAclass &MSA_) :
		ProblemBase(MSA_), ea(NULL), pa(NULL), lpa(NULL) {

	lsingle = 0.01;
	dim1body = MSA->ncol * MSAclass::NAA;

	dim = dim1body;

	Allocate();

	size_t ntemp = 2 * MSA->nrow * MSA->NAA * MSA->ncol + MSA->nrow * MSA->ncol;

	printf("#  vars to minimize: %lu (%.1fMB)\n", dim, 8.0 * dim / 1024 / 1024);
	printf("#         temp vars: %lu (%.1fMB)\n", ntemp,
			8.0 * ntemp / 1024 / 1024);

}

Problem1Body::Problem1Body(const Problem1Body &source) :
		ProblemBase(source), lsingle(source.lsingle), dim1body(source.dim1body), ea(
		NULL), pa(NULL), lpa(NULL) {

	Allocate();

	memcpy(ea, source.ea, MSA->nrow * MSA->NAA * MSA->ncol * sizeof(double));
	memcpy(pa, source.pa, MSA->nrow * MSA->NAA * MSA->ncol * sizeof(double));
	memcpy(lpa, source.lpa, MSA->nrow * MSA->ncol * sizeof(double));

}

Problem1Body::~Problem1Body() {

	Free();

}

void Problem1Body::Allocate() {

	ea = (double*) malloc(MSA->nrow * MSA->NAA * MSA->ncol * sizeof(double));
	pa = (double*) malloc(MSA->nrow * MSA->NAA * MSA->ncol * sizeof(double));
	lpa = (double*) malloc(MSA->nrow * MSA->ncol * sizeof(double));

}

void Problem1Body::Free() {

	free(ea);
	free(pa);
	free(lpa);

}

Problem1Body& Problem1Body::operator=(const Problem1Body &source) {

	assert(this != &source);

	FreeBase();

	dim = source.dim;
	dim1body = source.dim1body;

	MSA = source.MSA;

	MSA = source.MSA;

	AllocateBase();
	Allocate();

	memcpy(we, source.we, MSA->ncol * MSA->ncol * sizeof(double));

	memcpy(ea, source.ea, MSA->nrow * MSA->NAA * MSA->ncol * sizeof(double));
	memcpy(pa, source.pa, MSA->nrow * MSA->NAA * MSA->ncol * sizeof(double));
	memcpy(lpa, source.lpa, MSA->nrow * MSA->ncol * sizeof(double));

	return *this;

}

double Problem1Body::f(const double *x) {

	double f = 0.0;

	size_t ncol = MSA->ncol;
	size_t nrow = MSA->nrow;
	size_t NAA = MSAclass::NAA;

	/* loop over all sequences in the MSA */
#if defined(_OPENMP)
#pragma omp parallel for reduction (+:f)
#endif
	for (size_t i = 0; i < nrow; i++) {

		/* sequence weight */
		double weight = MSA->weight[i];

		/* current sequence */
		unsigned char *seq = msa + i * ncol;

		/* precompute energies of every letter
		 * at every position in the sequence */
		double *e = (double*) malloc(NAA * ncol * sizeof(double));

		/* logarithm of local partition functions
		 * (aka one-site pseudo-log-likelihoods
		 * or local free energies) */
		double *lp = (double*) malloc(ncol * sizeof(double));

		/* initialize energies with local fields */
		memcpy(e, x, ncol * (NAA - 1) * sizeof(double));

		/* fix the local fields for gaps at zero */
		memset(e + (NAA - 1) * ncol, 0, ncol * sizeof(double));

		/* compute local partition functions */
		memset(lp, 0, sizeof(double) * ncol);
		for (size_t a = 0; a < NAA; a++) {
			for (size_t s = 0; s < ncol; s++) {
				lp[s] += exp(e[a * ncol + s]);
			}
		}

		for (size_t s = 0; s < ncol; s++) {
			lp[s] = log(lp[s]);
		}

		/* update the objective function */
		for (size_t k = 0; k < ncol; k++) {
			f += weight * (-e[seq[k] * ncol + k] + lp[k]);
		}

		free(e);
		free(lp);

	}

	/* regularization */
	double reg = 0.0;
#if defined(_OPENMP)
#pragma omp parallel for reduction (+:reg)
#endif
	for (size_t v = 0; v < dim1body; v++) {
		reg += lsingle * x[v] * x[v];
	}

	f += reg;

	return f;

}

void Problem1Body::fdf(const double *x, double *f, double *g) {

	size_t ncol = MSA->ncol;
	size_t nrow = MSA->nrow;
	size_t NAA = MSAclass::NAA;

	memset(lpa, 0, nrow * ncol * sizeof(double));

	/* loop over all sequences in the MSA */

#if defined(_OPENMP)
#pragma omp parallel for
#endif
	for (size_t i = 0; i < nrow; i++) {

		/* precompute energies of every letter
		 * at every position in the sequence */
		double *e = ea + i * NAA * ncol;

		/* logarithm of local partition functions
		 * (aka one-site pseudo-log-likelihoods
		 * or local free energies) */
		double *lp = lpa + i * ncol;

		/* local probabilities of a every letter
		 * at every position in the sequence*/
		double *p = pa + i * NAA * ncol;

		/* initialize energies with local fields */
		memcpy(e, x, ncol * (NAA - 1) * sizeof(double));

		/* fix the local fields for gaps at zero */
		memset(e + (NAA - 1) * ncol, 0, ncol * sizeof(double));

		/* compute local partition functions */
		for (size_t a = 0; a < NAA; a++) {
			for (size_t s = 0; s < ncol; s++) {
				lp[s] += exp(e[a * ncol + s]);
			}
		}

		for (size_t s = 0; s < ncol; s++) {
			lp[s] = log(lp[s]);
		}

		/* compute local probabilities */
		for (size_t a = 0; a < NAA; a++) {
			for (size_t s = 0; s < ncol; s++) {
				p[a * ncol + s] = exp(e[a * ncol + s] - lp[s]);
			}
		}

	}

	/* compute objective function */
	{
		double obj = 0;

#if defined(_OPENMP)
#pragma omp parallel for reduction (+:obj)
#endif
		for (size_t i = 0; i < nrow; i++) {

			double *lp = lpa + i * ncol;
			unsigned char *s = msa + i * ncol;
			double *e = ea + i * NAA * ncol;
			double f = 0.0;
			for (size_t k = 0; k < ncol; k++) {
				f += MSA->weight[i] * (*lp++ - e[*s++ * ncol + k]);
			}
			obj += f;
		}
		*f = obj;
	}

	/* set gradient to 0 */
	memset(g, 0, sizeof(double) * dim);

	/* compute f and derivatives of h[] */
#if defined(_OPENMP)
#pragma omp parallel for
#endif
	for (size_t k = 0; k < ncol; k++) {

		for (size_t i = 0; i < nrow; i++) {

			double weight = MSA->weight[i];
			double *p = pa + i * NAA * ncol;
			unsigned char xik = (msa + i * ncol)[k];

			if (xik < NAA - 1) {
				g[xik * ncol + k] -= weight;
			}

			for (size_t a = 0; a < NAA - 1; a++) {
				g[a * ncol + k] += weight * p[a * ncol + k];
			}
		}
	}

	double reg = 0.0;

	/* regularize h */

#if defined(_OPENMP)
#pragma omp parallel for reduction (+:reg)
#endif
	for (size_t v = 0; v < dim1body; v++) {
		reg += lsingle * x[v] * x[v];
		g[v] += 2.0 * lsingle * x[v];
	}

	*f += reg;

}

void Problem1Body::df(const double *x, double *g) {

	double f;
	fdf(x, &f, g);

}

void Problem1Body::GetMRFvector(const double *x, double *mrfx) {

	/* */

}
