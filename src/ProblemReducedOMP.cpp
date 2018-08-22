/*
 * ProblemReducedOMPOMP.cpp
 *
 *  Created on: Mar 22, 2018
 *      Author: aivan
 */

#include <cstring>
#include <cmath>
#include <cassert>

#include "ProblemReducedOMP.h"

#include <omp.h>

ProblemReducedOMP::ProblemReducedOMP() :
		ProblemBase(), lsingle(0.0), lpair(0.0), dim1body(0), dim2body(0), dim2reduced(
				0), gaux(
		NULL), ea(NULL), pa(NULL), lpa(NULL) {

	/* nothing to be done */

}

ProblemReducedOMP::ProblemReducedOMP(const MSAclass &MSA_) :
		ProblemBase(MSA_), gaux(NULL), ea(NULL), pa(NULL), lpa(NULL) {

	lsingle = 0.01;
	lpair = 0.2 * (MSA->ncol - 1);

	dim1body = MSA->ncol * MSAclass::NAA;
	dim2body = dim1body * dim1body;

	dim2reduced = dim1body * (MSA->ncol - 1) * MSAclass::NAA / 2;

	dim = dim1body + dim2reduced;

	Allocate();

	size_t ntemp = dim2body + 2 * MSA->nrow * MSA->NAA * MSA->ncol
			+ MSA->nrow * MSA->ncol;

	printf("#  vars to minimize: %lu (%.1fMB)\n", dim, 8.0 * dim / 1024 / 1024);
	printf("#         temp vars: %lu (%.1fMB)\n", ntemp,
			8.0 * ntemp / 1024 / 1024);

}

ProblemReducedOMP::ProblemReducedOMP(const ProblemReducedOMP &source) :
		ProblemBase(source), lsingle(source.lsingle), lpair(source.lpair), dim1body(
				source.dim1body), dim2body(source.dim2body), dim2reduced(
				source.dim2reduced), gaux(NULL), ea(
		NULL), pa(NULL), lpa(NULL) {

	Allocate();

	memcpy(gaux, source.gaux, dim2body * sizeof(double));
	memcpy(ea, source.ea, MSA->nrow * MSA->NAA * MSA->ncol * sizeof(double));
	memcpy(pa, source.pa, MSA->nrow * MSA->NAA * MSA->ncol * sizeof(double));
	memcpy(lpa, source.lpa, MSA->nrow * MSA->ncol * sizeof(double));

}

ProblemReducedOMP::~ProblemReducedOMP() {

	Free();

}

void ProblemReducedOMP::Allocate() {

	gaux = (double*) malloc(dim2body * sizeof(double));
	ea = (double*) malloc(MSA->nrow * MSA->NAA * MSA->ncol * sizeof(double));
	pa = (double*) malloc(MSA->nrow * MSA->NAA * MSA->ncol * sizeof(double));
	lpa = (double*) malloc(MSA->nrow * MSA->ncol * sizeof(double));

}

void ProblemReducedOMP::Free() {

	free(gaux);
	free(ea);
	free(pa);
	free(lpa);

}

ProblemReducedOMP& ProblemReducedOMP::operator=(
		const ProblemReducedOMP &source) {

	assert(this != &source);

	FreeBase();
	Free();

	dim = source.dim;
	dim1body = source.dim1body;
	dim2body = source.dim2body;
	dim2reduced = source.dim2reduced;

	MSA = source.MSA;

	AllocateBase();
	Allocate();

	memcpy(we, source.we, MSA->ncol * MSA->ncol * sizeof(double));

	memcpy(gaux, source.gaux, dim2body * sizeof(double));
	memcpy(ea, source.ea, MSA->nrow * MSA->NAA * MSA->ncol * sizeof(double));
	memcpy(pa, source.pa, MSA->nrow * MSA->NAA * MSA->ncol * sizeof(double));
	memcpy(lpa, source.lpa, MSA->nrow * MSA->ncol * sizeof(double));

	return *this;

}

void ProblemReducedOMP::df(const double *x, double *g) {

	double f;
	fdf(x, &f, g);

}

double ProblemReducedOMP::f(const double *x) {

	double f = 0.0;

	size_t ncol = MSA->ncol;
	size_t nrow = MSA->nrow;
	size_t NAA = MSAclass::NAA;

	const double *x1 = x; /* local fields Vi */
	const double *x2 = x + dim1body; /* couplings Wij */

	/* unroll 2-body terms into symmetric 2d matrix (redundant) */
	memset(gaux, 0, dim2body * sizeof(double));
#if defined(_OPENMP)
#pragma omp parallel for
#endif
	for (size_t k = 0; k < ncol * (ncol - 1) / 2; k++) {
		size_t i, j;
		To2D(k, i, j);
		const double *x2_k = x2 + k * NAA * NAA;
		for (size_t a = 0; a < NAA; a++) {
			for (size_t b = 0; b < NAA; b++) {
				double v = *x2_k++;
				gaux[((a * ncol + j) * NAA + b) * ncol + i] = v;
				gaux[((b * ncol + i) * NAA + a) * ncol + j] = v;
			}
		}
	}

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
		memcpy(e, x1, ncol * (NAA - 1) * sizeof(double));

		/* fix the local fields for gaps at zero */
		memset(e + (NAA - 1) * ncol, 0, ncol * sizeof(double));

		/* add interactions with all other positions */
		for (size_t k = 0; k < ncol; k++) {

			/* wp[] - array of interaction energies of res k
			 * of identity seq[k] with all other positions
			 * of varying identities*/
			const double *wp = gaux + (seq[k] * ncol + k) * NAA * ncol;
			double *ep = e;
			for (size_t j = 0; j < NAA * ncol; j++) {
				*ep++ += *wp++;
			}
		}

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

#if defined(_OPENMP)
#pragma omp parallel for reduction (+:reg)
#endif
	for (size_t v = 0; v < dim2body; v++) {
		reg += 0.5 * lpair * gaux[v] * gaux[v];
	}

	f += reg;

	return f;

}

void ProblemReducedOMP::fdf(const double *x, double *f, double *g) {

	size_t ncol = MSA->ncol;
	size_t nrow = MSA->nrow;
	size_t NAA = MSAclass::NAA;

	const double *x1 = x; /* local fields Vi */
	const double *x2 = x + dim1body; /* couplings Wij */

	double *g1 = g;
	double *g2 = g + dim1body;

	memset(lpa, 0, nrow * ncol * sizeof(double));

	/* unroll 2-body terms into symmetric 2d matrix (redundant) */
	memset(gaux, 0, dim2body * sizeof(double));
#if defined(_OPENMP)
#pragma omp parallel for
#endif
	for (size_t k = 0; k < ncol * (ncol - 1) / 2; k++) {
		size_t i, j;
		To2D(k, i, j);
		const double *x2_k = x2 + k * NAA * NAA;
		for (size_t a = 0; a < NAA; a++) {
			for (size_t b = 0; b < NAA; b++) {
				double v = *x2_k++;
				gaux[((a * ncol + i) * NAA + b) * ncol + j] = v;
				gaux[((b * ncol + j) * NAA + a) * ncol + i] = v;
			}
		}
	}

	/* loop over all sequences in the MSA */

#if defined(_OPENMP)
#pragma omp parallel for
#endif
	for (size_t i = 0; i < nrow; i++) {

		unsigned char *seq = msa + i * ncol;

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
		memcpy(e, x1, ncol * (NAA - 1) * sizeof(double));

		/* fix the local fields for gaps at zero */
		memset(e + (NAA - 1) * ncol, 0, ncol * sizeof(double));

		/* add interactions with all other positions */
		for (size_t k = 0; k < ncol; k++) {
			const double *wp = gaux + (seq[k] * ncol + k) * NAA * ncol;
			double *ep = e;
			for (size_t j = 0; j < NAA * ncol; j++) {
				*ep++ += *wp++;
			}
		}

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
	memset(gaux, 0, dim2body * sizeof(double));

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
				g1[xik * ncol + k] -= weight;
			}

			for (size_t a = 0; a < NAA - 1; a++) {
				g1[a * ncol + k] += weight * p[a * ncol + k];
			}
		}
	}

	/* derivatives of J[][] */
	for (size_t i = 0; i < nrow; i++) {

		double weight = MSA->weight[i];
		unsigned char *seq = msa + i * ncol;

#if defined(_OPENMP)
#pragma omp parallel for
#endif
		for (size_t k = 0; k < ncol; k++) {

			double *gaux_p = gaux + (seq[k] * ncol + k) * NAA * ncol;

			for (size_t j = 0; j < ncol; j++) {
				gaux_p[seq[j] * ncol + j] -= weight;
			}

			double *pp = pa + i * NAA * ncol;
			for (size_t j = 0; j < NAA * ncol; j++) {
				*gaux_p++ += weight * *pp++;
			}

		}

	}

	/* add transposed J[][] with untransposed and save
	 * the upper triangle in g[]
	 * (roll back to the reduced representation) */

#if defined(_OPENMP)
#pragma omp parallel for
#endif
	for (size_t k = 0; k < ncol * (ncol - 1) / 2; k++) {
		size_t i, j;
		To2D(k, i, j);
		double *g2_k = g2 + k * NAA * NAA;
		for (size_t a = 0; a < NAA; a++) {
			for (size_t b = 0; b < NAA; b++) {
				*g2_k++ = gaux[((a * ncol + i) * NAA + b) * ncol + j]
						+ gaux[((b * ncol + j) * NAA + a) * ncol + i];
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

	/* regularize J */

#if defined(_OPENMP)
#pragma omp parallel for reduction (+:reg)
#endif
	for (size_t v = 0; v < dim2reduced; v++) {
		reg += lpair * x2[v] * x2[v];
		g2[v] += 2.0 * lpair * x2[v];
	}

	*f += reg;

}

void ProblemReducedOMP::GetMRFvector(const double *x, double *mrfx) {

	memset(mrfx, 0, (dim1body + dim2body) * sizeof(double));

	size_t NAA = MSAclass::NAA;
	size_t ncol = MSA->ncol;

	for (size_t i = 0; i < ncol; i++) {
		for (size_t a = 0; a < NAA; a++) {
			mrfx[i * NAA + a] = x[a * ncol + i];
		}
	}

	double *J = mrfx + dim1body;

	for (size_t k = 0; k < ncol * (ncol - 1) / 2; k++) {
		size_t i, j;
		To2D(k, i, j);
		const double *x2_k = x + dim1body + k * NAA * NAA;
		for (size_t a = 0; a < NAA; a++) {
			for (size_t b = 0; b < NAA; b++) {
				double v = *x2_k++;
				J[(i * ncol + j) * NAA * NAA + a * NAA + b] = v;
				J[(j * ncol + i) * NAA * NAA + b * NAA + a] = v;
			}
		}
	}

}

size_t ProblemReducedOMP::To1D(size_t i, size_t j) {

	size_t n = MSA->ncol;

	assert(i != j);

	/* make sure i < j */
	if (i > j) {
		size_t tmp = j;
		j = i;
		i = tmp;
	}

// https://stackoverflow.com/questions/27086195/linear-index-upper-triangular-matrix
	return (n * (n - 1) / 2) - (n - i) * ((n - i) - 1) / 2 + j - i - 1;

}

void ProblemReducedOMP::To2D(size_t k, size_t &i, size_t &j) {

	size_t n = MSA->ncol;

// https://stackoverflow.com/questions/27086195/linear-index-upper-triangular-matrix
	i = n - 2 - floor(sqrt(-8 * k + 4 * n * (n - 1) - 7) / 2.0 - 0.5);
	j = k + i + 1 - n * (n - 1) / 2 + (n - i) * ((n - i) - 1) / 2;

}

