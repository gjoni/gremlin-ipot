/*
 * ProblemRRCE.cpp
 *
 *  Created on: Jul 26, 2017
 *      Author: ivan
 */

#include <cstdlib>
#include <cstring>
#include <cmath>

#include "ProblemRRCE.h"
#include "EigenRRCE.h"

ProblemRRCE::ProblemRRCE() :
		ProblemBase(), pot(), nmodes(0), em(NULL), ev(NULL), lsingle(0), lpair(
				0), dim1body(0), dim2body(0) {

	/* */

}

ProblemRRCE::~ProblemRRCE() {

	Free();

}

void ProblemRRCE::Allocate() {

	em = (double*) malloc(nmodes * 20 * 20 * sizeof(double));
	ev = (double*) malloc(nmodes * sizeof(double));

}

void ProblemRRCE::Free() {

	free(em);
	free(ev);

}

ProblemRRCE::ProblemRRCE(const MSAclass &MSA_, size_t n) :
		ProblemBase(MSA_), pot("data/RRCE20SCC/table.7.0A_k5"), nmodes(n), em(
		NULL), ev(NULL), lsingle(0), lpair(0), dim1body(0), dim2body(0) {

	Allocate();

	/*
	 * do eigendecomposition of J[][]
	 */
	double **J = (double**) malloc(20 * sizeof(double*));
	for (int i = 0; i < 20; i++) {
		J[i] = (double*) malloc(20 * sizeof(double));
	}

	pot.GetCouplings(J);
	EigenRRCE Eigen(J);

	printf("# RRCE reconstruction error / correl:  %.5f / %.5f\n",
			Eigen.GetReconstructionError(nmodes),
			Eigen.GetReconstructionCorrel(nmodes));

	printf("# Eigenvalues (squared):");
	for (size_t n = 0; n < nmodes; n++) {
		Eigen.GetEigenmatrix(n, em + n * 400);
		ev[n] = Eigen.GetEigenvalue(n);
		printf(" %.5f", ev[n] * ev[n]);
	}
	printf("\n");

	/*
	 * set regularization params
	 */
	lsingle = 0.01;
	lpair = 0.01; //0.2 * (MSA->ncol - 1);

	/*
	 * set dimensions
	 */
	dim1body = MSA->ncol * MSAclass::NAA;
	dim2body = MSA->ncol * MSA->ncol * nmodes;
	dim = dim1body + dim2body;

	/*
	 * free
	 */
	for (int i = 0; i < 20; i++) {
		free(J[i]);
	}
	free(J);

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
 * couplings W={wi,j(nmode)}, 1D array of (NAA * NAA * nmodes) size
 * i,j - positions in sequence
 * AAi,AAj - amino acid identities at positions i,j
 *
 *    w0,0(mode_1), w0,1(mode_1), ..., w0,L-1(mode_1),
 *    w0,0(mode_2), w0,1(mode_2), ..., w0,L-1(mode_2),
 *    ...,
 *    w0,0(mode_n), w0,1(mode_n), ..., w0,L-1(mode_n),
 *
 *    w1,0(mode_1), w1,1(mode_1), ..., w1,L-1(mode_1),
 *    w1,0(mode_2), w1,1(mode_2), ..., w1,L-1(mode_2),
 *    ...,
 *    w1,0(mode_n), w1,1(mode_n), ..., w1,L-1(mode_n),
 *
 *    ...
 */

double ProblemRRCE::f(const double *x) {

	double f = 0.0;

	size_t ncol = MSA->ncol;
	size_t nrow = MSA->nrow;
	size_t NAA = MSAclass::NAA;

	const double *x1 = x; /* local fields Vi */
	const double *x2 = x + dim1body; /* couplings Wij */

	/* loop over all sequences in the MSA */
	for (size_t i = 0; i < nrow; i++) {

		/* sequence weight */
		double weight = w[i];

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

			for (size_t j = 0; j < ncol; j++) {

				/* interaction coefficients for edge k -> j */
				const double *ckj = x2 + (k * ncol + j) * nmodes;

				/* scan through all AA identities at position j */
				for (size_t a = 0; a < NAA - 1; a++) {

					/* add up contributions from all eigenmodes */
					for (size_t n = 0; n < nmodes; n++) {
						e[a * ncol + k] += ckj[n]
								* em[n * 400 + a * 20 + seq[j]];
					}
				}
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

	/* regularization - 1-body terms */
	double reg = 0.0;
	for (size_t v = 0; v < dim1body; v++) {
		reg += lsingle * x[v] * x[v];
	}

	/* regularization - 2-body terms */
	for (size_t ij = 0; ij < ncol * ncol; ij++) {
		const double *c = x2 + ij * nmodes;
		for (size_t k = 0; k < nmodes * 400; k++) {
			double r = em[k] * c[k / 400];
			reg += 0.5 * lpair * r * r;
		}
	}

	f += reg;

	return f;

}

void ProblemRRCE::df(const double *x, double *g) {

	double f = 0.0;
	fdf(x, &f, g);

}

void ProblemRRCE::fdf(const double *x, double *f, double *g) {

	size_t ncol = MSA->ncol;
	size_t nrow = MSA->nrow;
	size_t NAA = MSAclass::NAA;

	const double *x1 = x; /* local fields Vi */
	const double *x2 = x + dim1body; /* couplings Wij */

	double *g1 = g;
	double *g2 = g + dim1body;

	/* set fx and gradient to 0 initially */
	*f = 0.0;
	memset(g, 0, sizeof(double) * dim);

	/* aux array to store asymmetric 2-body gradient */
	double *gaux = (double*) calloc(dim2body, sizeof(double));

	/* loop over all sequences in the MSA */
	for (size_t i = 0; i < nrow; i++) {

		double weight = w[i];
		unsigned char *seq = msa + i * ncol;

		/* precomputed energies of every letter
		 * at every position in the sequence */
		double *e = (double*) malloc(NAA * ncol * sizeof(double));
		if (e == NULL) {
			printf("Error: not enough memory\n");
		}

		/* logarithm of local partition functions
		 * (aka one-site pseudo-log-likelihoods
		 * or local free energies) */
		double *lp = (double*) malloc(ncol * sizeof(double));
		if (lp == NULL) {
			printf("Error: not enough memory\n");
		}

		/* local probabilities of a every letter
		 * at every position in the sequence*/
		double *p = (double*) malloc(NAA * ncol * sizeof(double));
		if (p == NULL) {
			printf("Error: not enough memory\n");
		}

		/* initialize energies with local fields */
		memcpy(e, x1, ncol * (NAA - 1) * sizeof(double));

		/* fix the local fields for gaps at zero */
		memset(e + (NAA - 1) * ncol, 0, ncol * sizeof(double));

		/* add interactions with all other positions */
		for (size_t k = 0; k < ncol; k++) {

			for (size_t j = 0; j < ncol; j++) {

				/* interaction coefficients for edge k -> j */
				const double *ckj = x2 + (k * ncol + j) * nmodes;

				/* scan through all AA identities at position j */
				for (size_t a = 0; a < NAA - 1; a++) {

					/* add up contributions from all eigenmodes */
					for (size_t n = 0; n < nmodes; n++) {
						if (seq[j] < 20) {
							e[a * ncol + k] += ckj[n]
									* em[n * 400 + a * 20 + seq[j]];
						}
					}
				}
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

			if (seq[k] >= 20) {
				continue;
			}

			for (size_t j = 0; j < ncol; j++) {

				if (seq[j] >= 20) {
					continue;
				}

				for (size_t n = 0; n < nmodes; n++) {

					double *gg = gaux + (k * ncol + j) * nmodes + n;
					double *JJ = em + n * 400;
					*gg -= weight * JJ[seq[k] * 20 + seq[j]];

					for (size_t c = 0; c < 20; c++) {
						*gg += weight * JJ[c * 20 + seq[j]] * p[c * ncol + k];
					}

				}
			}
		}

		free(e);
		free(lp);
		free(p);

	}

	/* make derivatives of J[][] symmetric -
	 * add transposed onto untransposed */
	for (size_t i = 0; i < ncol; i++) {
		for (size_t j = 0; j < ncol; j++) {
			for (size_t k = 0; k < nmodes; k++) {
				g2[(i * ncol + j) * nmodes + k] = gaux[(i * ncol + j) * nmodes
						+ k] + gaux[(j * ncol + i) * nmodes + k];
			}
		}
		for (size_t k = 0; k < nmodes; k++) {
			g2[(i * ncol + i) * nmodes + k] = 0.0;
		}
	}

	free(gaux);

	double reg = 0.0;

	/* regularize h */
	for (size_t v = 0; v < dim1body; v++) {
		reg += lsingle * x[v] * x[v];
		g[v] += 2.0 * lsingle * x[v];
	}

	/* regularize J */
	for (size_t v = dim1body; v < dim; v++) {
		reg += 0.5 * lpair * x[v] * x[v];
		g[v] += 2.0 * lpair * x[v];
	}

//	for (size_t ij = 0; ij < ncol * ncol; ij++) {
//
//		const double *c = x2 + ij * nmodes;
//
//		double *rr = (double*) calloc(400, sizeof(double));
//		for (size_t pq = 0; pq < 400; pq++) {
//			for (size_t n = 0; n < nmodes; n++) {
//				rr[pq] += em[n * 400 + pq] * c[n];
//			}
//		}
//
//		for (size_t pq = 0; pq < 400; pq++) {
//			reg += 0.5 * lpair * rr[pq] * rr[pq];
//		}
//
//		for (size_t n = 0; n < nmodes; n++) {
//			for (size_t ab = 0; ab < 400; ab++) {
//				g2[ij * nmodes + n] += 2.0 * lpair * em[400 * nmodes + ab]
//						* rr[ab];
//			}
//		}
//
//		free(rr);
//	}

	*f += reg;

}

void ProblemRRCE::GetMRFvector(const double *x, double *mrfx) {

	size_t NAA = MSAclass::NAA;
	size_t ncol = MSA->ncol;

	size_t dim_mrf = ncol * NAA + ncol * ncol * NAA * NAA;

	memset(mrfx, 0, dim_mrf * sizeof(double));

	for (size_t i = 0; i < ncol; i++) {
		for (size_t a = 0; a < NAA; a++) {
			mrfx[i * NAA + a] = x[a * ncol + i];
		}
	}

	for (size_t i = 0; i < ncol; i++) {
		for (size_t j = 0; j < ncol; j++) {
			double *Jij = mrfx + ncol * NAA + (i * ncol + j) * NAA * NAA;
			const double *cij = x + dim1body + (i * ncol + j) * nmodes;
			for (size_t a = 0; a < NAA - 1; a++) {
				for (size_t b = 0; b < NAA - 1; b++) {
					for (size_t n = 0; n < nmodes; n++) {
						Jij[a * NAA + b] += cij[n] * em[n * 400 + (a * 20 + b)];
					}
				}
			}
		}
	}

}
