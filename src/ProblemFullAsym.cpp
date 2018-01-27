/*
 * ProblemFullAsym.cpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ivan
 */

#include <cstring>
#include <cmath>
#include <cassert>

#include "ProblemFullAsym.h"
#include <omp.h>

ProblemFullAsym::ProblemFullAsym() :
		ProblemBase(), lsingle(0.0), lpair(0.0), lskew(0.0), dim1body(0), dim2body(
				0), J(
		NULL), dx(NULL), dy(
		NULL) {

	Allocate();

}

ProblemFullAsym::ProblemFullAsym(const MSAclass &MSA_) :
		ProblemBase(MSA_), J(NULL), dx(NULL), dy(
		NULL) {

	lsingle = 0.01;
	lpair = 0.2 * (MSA->ncol - 1);
	lskew = 0.2 * (MSA->ncol - 1);

	dim1body = MSA->ncol * MSAclass::NAA;
	dim2body = MSA->ncol * MSAclass::NAA * MSA->ncol * MSAclass::NAA;
	dim = dim1body + dim2body;

	Allocate();

	SetPriors();

}

ProblemFullAsym::ProblemFullAsym(const ProblemFullAsym &source) :
		ProblemBase(source), lsingle(source.lsingle), lpair(source.lpair), lskew(
				source.lskew), dim1body(source.dim1body), dim2body(
				source.dim2body), J(NULL), dx(NULL), dy(
		NULL) {

	Allocate();

	for (int i = 0; i < 21; i++) {
		memcpy(J[i], source.J[i], 21 * sizeof(double));
		memcpy(dx[i], source.dx[i], 21 * sizeof(double));
		memcpy(dy[i], source.dy[i], 21 * sizeof(double));
	}

}

ProblemFullAsym::~ProblemFullAsym() {

	Free();

}

ProblemFullAsym& ProblemFullAsym::operator=(const ProblemFullAsym &source) {

	/* TODO: does not look like this constructor does the job */

	assert(this != &source); /* an attempt to assign Residue to itself */

	FreeBase();

	dim = source.dim;
	dim1body = source.dim1body;
	dim2body = source.dim2body;

	MSA = source.MSA;

	AllocateBase();

	memcpy(we, source.we, MSA->ncol * MSA->ncol * sizeof(double));

	for (int i = 0; i < 21; i++) {
		memcpy(J[i], source.J[i], 21 * sizeof(double));
		memcpy(dx[i], source.dx[i], 21 * sizeof(double));
		memcpy(dy[i], source.dy[i], 21 * sizeof(double));
	}

	return *this;

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
 * TODO: ??? isn't it W={wi,j(AAj,AAi)} ???
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

void ProblemFullAsym::df(const double *x, double *g) {

	double f;
	fdf(x, &f, g);

}

double ProblemFullAsym::f(const double *x) {

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

//	for (size_t v = dim1body; v < dim; v++) {
//		reg += 0.5 * lpair * x[v] * x[v];
//	}

	size_t v = dim1body;
	for (size_t i = 0; i < ncol; i++) {
		for (size_t a = 0; a < NAA; a++) {
			for (size_t j = 0; j < ncol; j++) {
				for (size_t b = 0; b < NAA; b++) {
					reg += 0.5 * Prior(x[v], a, b);
					v++;
				}
			}
		}
	}

	f += reg;

	return f;

}

void ProblemFullAsym::fdf(const double *x, double *f, double *g) {

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
//#pragma omp parallel for
	for (size_t i = 0; i < nrow; i++) {

		double weight = MSA->weight[i];
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
//#pragma omp critical
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
//#pragma omp critical
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

//#pragma omp parallel for ordered
	for (size_t b = 0; b < NAA; b++) {
		for (size_t k = 0; k < ncol; k++) {
			for (size_t a = 0; a < NAA; a++) {

				/* set gradients to zero for self-edges */
				g2[((b * ncol + k) * NAA + a) * ncol + k] = 0;

				/* set gradient for masked edges to zero */
				for (size_t j = 0; j < ncol; j++) {
					if (we[k * ncol + j] == false) {
						g2[((b * ncol + k) * NAA + a) * ncol + j] = 0.0;
					}
				}

			}
		}
	}

	free(gaux);

	/* symmetrize 21x21 submatrices */
//	for (size_t b = 0; b < NAA; b++) {
//		for (size_t k = 0; k < ncol; k++) {
//			for (size_t a = b + 1; a < NAA; a++) {
//				for (size_t j = 0; j < ncol; j++) {
//					double s =
//							0.5
//									* (g2[((b * ncol + k) * NAA + a) * ncol + j]
//											+ g2[((a * ncol + k) * NAA + b)
//													* ncol + j]);
//					g2[((b * ncol + k) * NAA + a) * ncol + j] = s;
//					g2[((a * ncol + k) * NAA + b) * ncol + j] = s;
//				}
//			}
//
//		}
//	}
	double reg = 0.0;

	/* regularize h */
//#pragma omp parallel for ordered reduction (+:reg)
	for (size_t v = 0; v < dim1body; v++) {
		reg += lsingle * x[v] * x[v];
		g[v] += 2.0 * lsingle * x[v];
	}

	/* regularize J */
//#pragma omp parallel for ordered reduction (+:reg)
//	for (size_t v = dim1body; v < dim; v++) {
//		reg += 0.5 * lpair * x[v] * x[v];
//		g[v] += 2.0 * lpair * x[v];
//	}
	size_t v = dim1body;
	for (size_t i = 0; i < ncol; i++) {
		for (size_t a = 0; a < NAA; a++) {
			for (size_t j = 0; j < ncol; j++) {
				for (size_t b = 0; b < NAA; b++) {
					reg += 0.5 * Prior(x[v], a, b);
					g[v] += dPrior(x[v], a, b);
//					g[v] += 2.0 * lpair * x[v];
					v++;
				}
			}
		}
	}

	*f += reg;

}

void ProblemFullAsym::GetMRFvector(const double *x, double *mrfx) {

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

void ProblemFullAsym::Allocate() {

	J = (double**) malloc(21 * sizeof(double*));
	dx = (double**) malloc(21 * sizeof(double*));
	dy = (double**) malloc(21 * sizeof(double*));

	for (int i = 0; i < 21; i++) {
		J[i] = (double*) malloc(21 * sizeof(double));
		dx[i] = (double*) malloc(21 * sizeof(double));
		dy[i] = (double*) malloc(21 * sizeof(double));
	}

}

void ProblemFullAsym::Free() {

	if (J != NULL) {
		for (int i = 0; i < 21; i++) {
			free(J[i]);
		}
		free(J);
	}

	if (dx != NULL) {
		for (int i = 0; i < 21; i++) {
			free(dx[i]);
		}
		free(dx);
	}

	if (dy != NULL) {
		for (int i = 0; i < 21; i++) {
			free(dy[i]);
		}
		free(dy);
	}

}

void ProblemFullAsym::SetPriors() {

	/* set to zero */
	for (int i = 0; i < 21; i++) {
		for (int j = 0; j < 21; j++) {
			J[i][j] = 0.0;
			dx[i][j] = 0.0;
			dy[i][j] = 0.0;
		}
	}

	/* load statistical potential */
	RRCE RRCE_(RRCE::RRCE20RC, 7.8, 5);

	/* copy J[][] */
	RRCE_.GetCouplings(J);

	/* set dx[][], dy[][] */
	for (size_t a = 0; a < 20; a++) {
		for (size_t b = 0; b < 20; b++) {

			/* reverse coupling because we learn preferences,
			 * not energies */
			J[a][b] = -J[a][b];

			double x0 = 0.0;

			/* Newton iterations */
			for (int i = 0; i < 100; i++) {
				double x = x0 - dfPrior(x0, J[a][b]) / ddfPrior(x0, J[a][b]);
//				printf("%d %f %f\n", i, x, x0);
				x0 = x;
			}

			dx[a][b] = x0;
			dy[a][b] = fPrior(x0, J[a][b]);

//			printf("%6.3f%6.3f%7.3f ", J[a][b], dx[a][b],
//					fPrior(dx[a][b], J[a][b]) - dy[a][b]);

		}
//		printf("\n");
	}

}

double ProblemFullAsym::fPrior(double x, double mu) {

	return (lpair * x * x + lskew * exp(mu * x));

}

double ProblemFullAsym::dfPrior(double x, double mu) {

	return (2.0 * lpair * x + lskew * mu * exp(mu * x));

}

double ProblemFullAsym::ddfPrior(double x, double mu) {

	return (2.0 * lpair + lskew * mu * mu * exp(mu * x));

}

double ProblemFullAsym::Prior(double x, size_t a, size_t b) {

//	return (lpair * x * x);
	if (a >= 20 || b >= 20) {
		return (lpair * x * x);
	} else {
		return (fPrior(x + dx[a][b], J[a][b]) - dy[a][b]);
	}

}

double ProblemFullAsym::dPrior(double x, size_t a, size_t b) {

//	return (2.0 * lpair * x);
	if (a >= 20 || b >= 20) {
		return (2.0 * lpair * x);
	} else {
		return (dfPrior(x + dx[a][b], J[a][b]));
	}

}

void ProblemFullAsym::SetLpair(double l) {

	lpair = l * (MSA->ncol - 1);
	SetPriors();

}

void ProblemFullAsym::SetLskew(double l) {

	lskew = l * (MSA->ncol - 1);
	SetPriors();

}

void ProblemFullAsym::SetL2(double lpair_, double lskew_) {

	lpair = lpair_ * (MSA->ncol - 1);
	lskew = lskew_ * (MSA->ncol - 1);
	printf("# lpair= %f (%f)\n", lpair_, lpair);
	printf("# lskew= %f (%f)\n", lskew_, lskew);
	SetPriors();

}
