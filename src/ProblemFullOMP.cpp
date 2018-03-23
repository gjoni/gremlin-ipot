/*
 * ProblemFullOMPOMP.cpp
 *
 *  Created on: Mar 22, 2018
 *      Author: aivan
 */

#include <cstring>
#include <cmath>
#include <cassert>

#include "ProblemFullOMP.h"

#include <omp.h>

ProblemFullOMP::ProblemFullOMP() :
		ProblemBase(), lsingle(0.0), lpair(0.0), dim1body(0), dim2body(0), gaux(
		NULL), ea(NULL), pa(NULL), lpa(NULL) {

	/* nothing to be done */

}

ProblemFullOMP::ProblemFullOMP(const MSAclass &MSA_) :
		ProblemBase(MSA_), gaux(NULL), ea(NULL), pa(NULL), lpa(NULL) {

	lsingle = 0.01;
	lpair = 0.2 * (MSA->ncol - 1);

	dim1body = MSA->ncol * MSAclass::NAA;
	dim2body = MSA->ncol * MSAclass::NAA * MSA->ncol * MSAclass::NAA;
	dim = dim1body + dim2body;

	Allocate();

}

ProblemFullOMP::ProblemFullOMP(const ProblemFullOMP &source) :
		ProblemBase(source), lsingle(source.lsingle), lpair(source.lpair), dim1body(
				source.dim1body), dim2body(source.dim2body), gaux(NULL), ea(
		NULL), pa(NULL), lpa(NULL) {

	Allocate();

	memcpy(gaux, source.gaux, dim2body * sizeof(double));
	memcpy(ea, source.ea, MSA->nrow * MSA->NAA * MSA->ncol * sizeof(double));
	memcpy(pa, source.pa, MSA->nrow * MSA->NAA * MSA->ncol * sizeof(double));
	memcpy(lpa, source.lpa, MSA->nrow * MSA->ncol * sizeof(double));

}

ProblemFullOMP::~ProblemFullOMP() {

	Free();

}

void ProblemFullOMP::Allocate() {

	gaux = (double*) malloc(dim2body * sizeof(double));
	ea = (double*) malloc(MSA->nrow * MSA->NAA * MSA->ncol * sizeof(double));
	pa = (double*) malloc(MSA->nrow * MSA->NAA * MSA->ncol * sizeof(double));
	lpa = (double*) malloc(MSA->nrow * MSA->ncol * sizeof(double));

}

void ProblemFullOMP::Free() {

	free(gaux);
	free(ea);
	free(pa);
	free(lpa);

}

ProblemFullOMP& ProblemFullOMP::operator=(const ProblemFullOMP &source) {

	assert(this != &source); /* an attempt to assign Residue to itself */

	FreeBase();
	Free();

	dim = source.dim;
	dim1body = source.dim1body;
	dim2body = source.dim2body;

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

void ProblemFullOMP::df(const double *x, double *g) {

	double f;
	fdf(x, &f, g);

}

double ProblemFullOMP::f(const double *x) {

	double f = 0.0;

	size_t ncol = MSA->ncol;
	size_t nrow = MSA->nrow;
	size_t NAA = MSAclass::NAA;

	const double *x1 = x; /* local fields Vi */
	const double *x2 = x + dim1body; /* couplings Wij */

	/* loop over all sequences in the MSA */
	for (size_t i = 0; i < nrow; i++) {

		/* sequence weight */
		double weight = MSA->weight[i];

		/* current sequence */
		unsigned char *seq = msa + i * ncol;

		/* precomputed energies of every letter
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
			const double *wp = x2 + (seq[k] * ncol + k) * NAA * ncol;
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
	for (size_t v = 0; v < dim1body; v++) {
		reg += lsingle * x[v] * x[v];
	}

	for (size_t v = dim1body; v < dim; v++) {
		reg += 0.5 * lpair * x[v] * x[v];
	}

	f += reg;

	return f;

}

void ProblemFullOMP::fdf(const double *x, double *f, double *g) {

	size_t ncol = MSA->ncol;
	size_t nrow = MSA->nrow;
	size_t NAA = MSAclass::NAA;

	const double *x1 = x; /* local fields Vi */
	const double *x2 = x + dim1body; /* couplings Wij */

	double *g1 = g;
	double *g2 = g + dim1body;

	/* set gradient to 0 */
	memset(g, 0, sizeof(double) * dim);
	memset(gaux, 0, dim2body * sizeof(double));

	memset(lpa, 0, nrow * ncol * sizeof(double));

	/* loop over all sequences in the MSA */
#if defined(_OPENMP)
#pragma omp parallel for
#endif
	for (size_t i = 0; i < nrow; i++) {

		unsigned char *seq = msa + i * ncol;

		/* precomputed energies of every letter
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
			const double *wp = x2 + (seq[k] * ncol + k) * NAA * ncol;
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
		*f = 0.0;
		double *lp = lpa;
		unsigned char *s = msa;
		for (size_t i = 0; i < nrow; i++) {
			double *e = ea + i * NAA * ncol;
			for (size_t k = 0; k < ncol; k++) {
				*f += MSA->weight[i] * (*lp++ - e[*s++ * ncol + k]);
			}
		}
	}

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

	/* make derivatives of J[][] symmetric -
	 * add transposed onto untransposed */
	double *gaux_p = gaux;
	double *g2_p = g2;
	for (size_t b = 0; b < NAA; b++) {
		for (size_t k = 0; k < ncol; k++) {
			for (size_t a = 0; a < NAA; a++) {
				for (size_t j = 0; j < ncol; j++) {
					*g2_p++ = *gaux_p++
							+ gaux[((a * ncol + j) * NAA + b) * ncol + k];
				}
			}
		}
	}

//#if defined(_OPENMP)
//#pragma omp parallel for
//#endif
	for (size_t b = 0; b < NAA; b++) {
		for (size_t k = 0; k < ncol; k++) {
			for (size_t a = 0; a < NAA; a++) {

				/* set gradients to zero for self-edges */
				g2[((b * ncol + k) * NAA + a) * ncol + k] = 0;

				/* set gradient for masked edges to zero */
//				for (size_t j = 0; j < ncol; j++) {
//					if (we[k * ncol + j] == false) {
//						g2[((b * ncol + k) * NAA + a) * ncol + j] = 0.0;
//					}
//				}
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
	for (size_t v = dim1body; v < dim; v++) {
		reg += 0.5 * lpair * x[v] * x[v];
		g[v] += 2.0 * lpair * x[v];
	}

	*f += reg;

}

void ProblemFullOMP::GetMRFvector(const double *x, double *mrfx) {

	memset(mrfx, 0, dim * sizeof(double));

	size_t NAA = MSAclass::NAA;
	size_t ncol = MSA->ncol;

	for (size_t i = 0; i < ncol; i++) {
		for (size_t a = 0; a < NAA; a++) {
			mrfx[i * NAA + a] = x[a * ncol + i];
		}
	}

	double *J = mrfx + dim1body;
	const double *x2 = x + dim1body;

	for (size_t i = 0; i < ncol; i++) {
		for (size_t a = 0; a < NAA; a++) {
			for (size_t j = 0; j < ncol; j++) {
				for (size_t b = 0; b < NAA; b++) {
					J[(i * ncol + j) * NAA * NAA + a * NAA + b] = x2[((a * ncol
							+ i) * NAA + b) * ncol + j];
				}
			}
		}
	}

}
