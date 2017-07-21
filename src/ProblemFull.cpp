/*
 * ProblemFull.cpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ivan
 */

#include <cstring>
#include <cmath>

#include "ProblemFull.h"

ProblemFull::ProblemFull() :
		ProblemBase(), lsingle(0.0), lpair(0.0), gaux(NULL) {

	/* nothing to be done */

}

ProblemFull::ProblemFull(MSAclass &MSA_) :
		ProblemBase(MSA_) {

	lsingle = 0.01;
	lpair = 0.2 * (MSA->ncol - 1);

	dim = MSA->ncol * MSAclass::NAA * (1 + MSA->ncol * MSAclass::NAA)
	/* - MSA->ncol */;

	Allocate();

}

ProblemFull::~ProblemFull() {

	Free();

}

void ProblemFull::Allocate() {

	gaux = (double*) malloc(dim * sizeof(double));

}

void ProblemFull::Free() {

	free(gaux);

}

/* alignment of variables
 *
 * local fields V={vi(AA)}, 1D array of (NAA x ncol) size
 * i - position in sequence, AA - amino acid identity
 *
 *    v0('A'), v1('A'), ..., vL-1('A'),
 *    v0('R'), v1('R'), ..., vL-1('R'),
 *    ...,
 *    v0('-'), v1('-'), ..., vL-1('-')
 *
 * couplings W={wi,j(AAi,AAj)}, 1D array of (NAA * NAA * ncol * ncol) size
 * i,j - positions in sequence
 * AAi,AAj - amino acid identities at positions i,j
 *
 *    w0,0('A','A'), w0,1('A','A'), ..., w0,L-1('A','A'),
 *    w0,0('A','R'), w0,1('A','R'), ..., w0,L-1('A','R'),
 *    ...,
 *    w0,0('A','-'), w0,1('A','-'), ..., w0,L-1('A','-'),
 *
 *    w0,0('R','A'), w0,1('R','A'), ..., w0,L-1('R','A'),
 *    w0,0('R','R'), w0,1('R','R'), ..., w0,L-1('R','R'),
 *    ...,
 *    w0,0('R','-'), w0,1('R','-'), ..., w0,L-1('R','-'),
 *
 *    w0,0('-','A'), w0,1('-','A'), ..., w0,L-1('-','A'),
 *    w0,0('-','R'), w0,1('-','R'), ..., w0,L-1('-','R'),
 *    ...,
 *    w0,0('-','-'), w0,1('-','-'), ..., w0,L-1('-','-'),
 *
 *    w1,0('A','A'), w1,1('A','A'), ..., w1,L-1('A','A'),
 *    w1,0('A','R'), w1,1('A','R'), ..., w1,L-1('A','R'),
 *    ...,
 *    w1,0('A','-'), w1,1('A','-'), ..., w1,L-1('A','-'),
 *    ...
 */

void ProblemFull::df(const gsl_vector *x, gsl_vector *g) {

	double f;
	fdf(x, &f, g);

}

double ProblemFull::f(const gsl_vector *x) {

	double f = 0.0;

	size_t ncol = MSA->ncol;
	size_t nrow = MSA->nrow;
	size_t NAA = MSAclass::NAA;

	size_t nsingle = ncol * NAA;

	const double *x1 = x->data; /* local fields Vi */
	const double *x2 = x->data + nsingle; /* couplings Wij */

	/* loop over all sequences in the MSA */
	for (size_t i = 0; i < nrow; i++) {

		/* sequence weight */
		double weight = w[i];

		/* current sequence */
		unsigned char *seq = MSA->msa + i * ncol;

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
	for (size_t v = 0; v < nsingle; v++) {
		reg += lsingle * x->data[v] * x->data[v];
	}

	for (size_t v = nsingle; v < dim; v++) {
		reg += 0.5 * lpair * x->data[v] * x->data[v];
	}

	f += reg;

	return f;

}

void ProblemFull::fdf(const gsl_vector *x, double *f, gsl_vector *g) {

	size_t ncol = MSA->ncol;
	size_t nrow = MSA->nrow;
	size_t NAA = MSAclass::NAA;

	size_t nsingle = ncol * NAA;
	size_t nvar = nsingle + ncol * ncol * NAA * NAA;

	const double *x1 = x->data; /* local fields Vi */
	const double *x2 = x->data + nsingle; /* couplings Wij */

	double *g1 = g->data;
	double *g2 = g->data + nsingle;

	/* set fx and gradient to 0 initially */
	*f = 0.0;
	memset(g->data, 0, sizeof(double) * nvar);
	memset(gaux, 0, sizeof(double) * (NAA * NAA * ncol * ncol));

	/* loop over all sequences in the MSA */
	for (size_t i = 0; i < nrow; i++) {

		double weight = w[i];
		unsigned char *seq = MSA->msa + i * ncol;

		/* precomputed energies of every letter
		 * at every position in the sequence */
		double *e = (double*) malloc(NAA * ncol * sizeof(double));

		/* logarithm of local partition functions
		 * (aka one-site pseudo-log-likelihoods
		 * or local free energies) */
		double *lp = (double*) malloc(ncol * sizeof(double));

		/* local probabilities of specific AA
		 * at every position in the sequence*/
		double *p = (double*) malloc(NAA * ncol * sizeof(double));

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
		memset(lp, 0, sizeof(double) * ncol);
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

		/* compute f and derivatives of h[] */
		for (size_t k = 0; k < ncol; k++) {

			unsigned char xik = seq[k];

			*f += weight * (-e[xik * ncol + k] + lp[k]);

			if (xik < NAA - 1) {
				g1[xik * ncol + k] -= weight;
			}

			for (size_t a = 0; a < NAA - 1; a++) {
				g1[a * ncol + k] += weight * p[a * ncol + k];
			}

		}

		/* derivatives of J[][] */
		for (size_t k = 0; k < ncol; k++) {

			double *gaux_p = gaux + (seq[k] * ncol + k) * NAA * ncol;

			for (size_t j = 0; j < ncol; j++) {
				gaux_p[seq[j] * ncol + j] -= weight;
			}

			double *pp = p;
			for (size_t j = 0; j < NAA * ncol; j++) {
				*gaux_p++ += weight * *pp++;
			}

		}

		free(e);
		free(lp);
		free(p);

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

	/* set gradients to zero for self-edges */
	for (size_t b = 0; b < NAA; b++) {
		for (size_t k = 0; k < ncol; k++) {
			for (size_t a = 0; a < NAA; a++) {
				g2[((b * ncol + k) * NAA + a) * ncol + k] = 0;
			}
		}
	}

	double reg = 0.0;

	/* regularize h */
	for (size_t v = 0; v < nsingle; v++) {
		reg += lsingle * x->data[v] * x->data[v];
		g->data[v] += 2.0 * lsingle * x->data[v];
	}

	/* regulariz J */
	for (size_t v = nsingle; v < nvar; v++) {
		reg += 0.5 * lpair * x->data[v] * x->data[v];
		g->data[v] += 2.0 * lpair * x->data[v];
	}

	*f += reg;

}
