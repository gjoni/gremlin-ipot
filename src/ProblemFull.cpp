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
		ProblemBase(), lsingle(0.0), lpair(0.0), gsymm(NULL) {

	/* nothing to be done */

}

ProblemFull::ProblemFull(MSAclass &MSA_) :
		ProblemBase(MSA_) {

	lsingle = 0.01;
	lpair = 0.2 * (MSA->ncol - 1);

	dim = MSA->ncol * MSAclass::NAA * (1 + MSA->ncol * MSAclass::NAA);

	Allocate();

}

ProblemFull::~ProblemFull() {

	Free();

}

void ProblemFull::Allocate() {

	gsymm = (double*) malloc(dim * sizeof(double));

}

void ProblemFull::Free() {

	free(gsymm);

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

double ProblemFull::f(const gsl_vector *x) {

	double f = 0.0;

	size_t ncol = MSA->ncol;
	size_t nrow = MSA->nrow;

	size_t NAA = MSAclass::NAA;

	size_t nsingle = ncol * (NAA - 1);
//	size_t nvar = nsingle + ncol * ncol * NAA * NAA;

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

		/* local probabilities of specific AA
		 * at every position in the sequence*/
//		double *p = (double*) malloc(NAA * ncol * sizeof(double));
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
//		for (size_t a = 0; a < NAA; a++) {
//			for (size_t s = 0; s < ncol; s++) {
//				p[a * ncol + s] = exp(e[a * ncol + s] - lp[s]);
//			}
//		}
		/* update the objective function */
		for (size_t k = 0; k < ncol; k++) {
			f += weight * (-e[seq[k] * ncol + k] + lp[k]);
		}

		free(e);
		free(lp);
//		free(p);

	}

	return f;

}

void ProblemFull::df(const gsl_vector *x, gsl_vector *g) {

	double f;
	fdf(x, &f, g);

}

void ProblemFull::fdf(const gsl_vector *x, double *f, gsl_vector *g) {

	*f = 0.0;

	size_t ncol = MSA->ncol;
	size_t nrow = MSA->nrow;

	size_t NAA = MSAclass::NAA;

	size_t nsingle = ncol * (MSAclass::NAA - 1);
	size_t nvar = nsingle + ncol * ncol * MSAclass::NAA * MSAclass::NAA;

	const double *x1 = x->data;
	const double *x2 = x->data + nsingle;

	double *g1 = g->data;
	double *g2 = g->data + nsingle;

// set fx and gradient to 0 initially
	memset(g->data, 0, sizeof(double) * nvar);
	memset(gsymm, 0, dim * sizeof(double));

	for (size_t i = 0; i < nrow; i++) {

		double weight = w[i];
		unsigned char *seq = MSA->msa + i * ncol;

		double *e = (double*) malloc(NAA * ncol * sizeof(double));
		double *lp = (double*) malloc(ncol * sizeof(double));
		double *p = (double*) malloc(NAA * ncol * sizeof(double));

		memcpy(e, x1, ncol * (NAA - 1) * sizeof(double));
		memset(e + (NAA - 1) * ncol, 0, ncol * sizeof(double));

		for (size_t k = 0; k < ncol; k++) {
			unsigned char xik = seq[k];
			const double *w = x2 + (xik * ncol + k) * NAA * ncol;
			double *ep = e;

			for (size_t j = 0; j < NAA * ncol; j++) {
				*ep++ += *w++;
			}
		}

		// compute precomp_sum(s) = log( sum(a=1..21) exp(PC(a,s)) )
		memset(lp, 0, sizeof(double) * ncol);
		for (size_t a = 0; a < NAA; a++) {
			for (size_t s = 0; s < ncol; s++) {
				lp[s] += exp(e[a * ncol + s]);
			}
		}

		for (size_t s = 0; s < ncol; s++) {
			lp[s] = log(lp[s]);
		}

		for (size_t a = 0; a < NAA; a++) {
			for (size_t s = 0; s < ncol; s++) {
				p[a * ncol + s] = exp(e[a * ncol + s] - lp[s]);
			}
		}

		// actually compute fx and gradient
		for (size_t k = 0; k < ncol; k++) {

			unsigned char xik = seq[k];

			*f += weight * (-e[seq[k] * ncol + k] + lp[k]);

			if (xik < NAA - 1) {
				g1[xik * ncol + k] -= weight;
			}

			for (size_t a = 0; a < NAA - 1; a++) {
				g1[a * ncol + k] += weight * p[a * ncol + k];
			}

		}

		for (size_t k = 0; k < ncol; k++) {

			double *g2p = g2 + (seq[k] * ncol + k) * NAA * ncol;

			for (size_t j = 0; j < ncol; j++) {
				g2p[seq[j] * ncol + j] -= weight;
			}

//			double *g2p = g2 + (seq[k] * ncol + k) * NAA * ncol;
			double *pp = p;
			for (size_t j = 0; j < NAA * ncol; j++) {
				*g2p++ += weight * *pp++;
			}
		}

		free(e);
		free(lp);
		free(p);

	} // i

	// add transposed onto un-transposed
	double *gsymm_p = gsymm;
	double *g2_p = g2;
	for (size_t b = 0; b < NAA; b++) {
		for (size_t k = 0; k < ncol; k++) {
			for (size_t a = 0; a < NAA; a++) {
				for (size_t j = 0; j < ncol; j++) {
//					G2L(b, k, a, j) = G2(b, k, a, j) + G2(a, j, b, k);
					*gsymm_p++ = *g2_p++
							+ g2[((a * ncol + j) * NAA + b) * ncol + k];
				}
			}
		}
	}

	// set gradients to zero for self-edges
	for (size_t b = 0; b < NAA; b++) {
		for (size_t k = 0; k < ncol; k++) {
			for (size_t a = 0; a < NAA; a++) {
//				G2L(b, k, a, k) = 0;
				gsymm[((b * ncol + k) * NAA + a) * ncol + k] = 0.0;
			}
		}
	}

	// regularization
	double reg = 0.0;
	for (size_t v = 0; v < nsingle; v++) {
		reg += lsingle * x->data[v] * x->data[v];
		g->data[v] += 2.0 * lsingle * x->data[v]; // F2 is 2.0
	}

	for (size_t v = nsingle; v < nvar; v++) {
		reg += 0.5 * lpair * x->data[v] * x->data[v]; // F05 is 0.5
		g->data[v] += 2.0 * lpair * x->data[v]; // F2 is 2.0
	}

	*f += reg;
}
