/*
 * ProblemSmoothL1.cpp
 *
 *  Created on: Aug 30, 2018
 *      Author: aivan
 */

#include <cassert>
#include <cmath>
#include <cstring>

#include "ProblemSmoothL1.h"

ProblemSmoothL1::ProblemSmoothL1() :
		ProblemL2(), eps(100.0) {

	/* */

}

ProblemSmoothL1::ProblemSmoothL1(const MSAclass &MSA) :
		ProblemL2(MSA), eps(100.0) {

	/* */

}

ProblemSmoothL1::ProblemSmoothL1(const ProblemL2 &source) :
		ProblemL2(source), eps(100.0) {

	/* */

}

ProblemSmoothL1::~ProblemSmoothL1() {

	/* */

}

ProblemSmoothL1& ProblemSmoothL1::operator=(const ProblemL2 &source) {

	assert(this != &source);

	/*
	 * TODO: too much memory overhead (3 copies of same data)
	 */

	ProblemSmoothL1 tmp(source);

	FreeBase();
	Free();

	dim = tmp.dim;
	dim1body = tmp.dim1body;
	dim2body = tmp.dim2body;
	dim2reduced = tmp.dim2reduced;
	eps = 10;

	MSA = tmp.MSA;

	size_t ncol = MSA->GetNcol();
	size_t nrow = MSA->GetNrow();
	size_t NAA = MSAclass::NAA;

	SetLsingle(tmp.lsingle);
	SetLpair(tmp.lpair);

	AllocateBase();
	Allocate();

	memcpy(we, tmp.we, ncol * ncol * sizeof(double));

	memcpy(gaux, tmp.gaux, dim2body * sizeof(double));
	memcpy(ea, tmp.ea, nrow * NAA * ncol * sizeof(double));
	memcpy(pa, tmp.pa, nrow * NAA * ncol * sizeof(double));
	memcpy(lpa, tmp.lpa, nrow * ncol * sizeof(double));

	return *this;

}

double ProblemSmoothL1::Reg_f(const double *x) {

	printf("Reg_f() call ...\n");

	double reg = 0.0;

#if defined(_OPENMP)
#pragma omp parallel for reduction (+:reg)
#endif
	for (size_t v = 0; v < dim1body; v++) {

//		double r = log(1.0 + exp(-eps * x[v])) + log(1.0 + exp(eps * x[v]));
//		reg += lsingle / eps * r;

		reg += lsingle * sqrt(x[v] * x[v] + eps);
	}

#if defined(_OPENMP)
#pragma omp parallel for reduction (+:reg)
#endif
	for (size_t v = 0; v < dim2body; v++) {

//		double r = log(1.0 + exp(-eps * gaux[v]))
//				+ log(1.0 + exp(eps * gaux[v]));
//		reg += 0.5 * lpair / eps * r;

		reg += 0.5 * lpair * sqrt(gaux[v] * gaux[v] + eps);

	}

	return reg;

	return 0;

}

double ProblemSmoothL1::Reg_fdf(const double *x, double *g) {

//	printf("Reg_fdf() call: eps= %f ...\n", eps);

	double reg = 0.0;

	const double *x2 = x + dim1body;
	double *g2 = g + dim1body;

	/* regularize h */
#if defined(_OPENMP)
#pragma omp parallel for reduction (+:reg)
#endif
	for (size_t v = 0; v < dim1body; v++) {
		double r1 = 1.0 + exp(-eps * x[v]);
		double r2 = 1.0 + exp(eps * x[v]);
		reg += lsingle / eps * (log(r1) + log(r2));
		g[v] += lsingle * (1.0 / r1 - 1.0 / r2);

//		reg += lsingle * sqrt(x[v] * x[v] + eps);
//		g[v] += lsingle * x[v] / sqrt(x[v] * x[v] + eps);
	}

	/* regularize J */
#if defined(_OPENMP)
#pragma omp parallel for reduction (+:reg)
#endif
	for (size_t v = 0; v < dim2reduced; v++) {
		double r1 = 1.0 + exp(-eps * x2[v]);
		double r2 = 1.0 + exp(eps * x2[v]);
		reg += lpair / eps * (log(r1) + log(r2));
		g2[v] += lpair * (1.0 / r1 - 1.0 / r2);

//		reg += lpair * sqrt(x2[v] * x2[v] + eps);
//		g2[v] += lpair * x2[v] / sqrt(x2[v] * x2[v] + eps);

	}

	return reg;

}

void ProblemSmoothL1::SetSmoothingFactor(double eps_) {

	eps = eps_;

}
