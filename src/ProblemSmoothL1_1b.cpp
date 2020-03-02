/*
 * ProblemSmoothL1_1b.cpp
 *
 *  Created on: Aug 30, 2018
 *      Author: aivan
 */

#include <cassert>
#include <cmath>
#include <cstring>

#include "ProblemSmoothL1_1b.h"

ProblemSmoothL1_1b::ProblemSmoothL1_1b() :
		ProblemL2_1b(), eps(100.0) {

	/* */

}

ProblemSmoothL1_1b::ProblemSmoothL1_1b(const MSAclass &MSA) :
		ProblemL2_1b(MSA), eps(100.0) {

	/* */

}

ProblemSmoothL1_1b::ProblemSmoothL1_1b(const ProblemL2_1b &source) :
		ProblemL2_1b(source), eps(100.0) {

	/* */

}

ProblemSmoothL1_1b::~ProblemSmoothL1_1b() {

	/* */

}

ProblemSmoothL1_1b& ProblemSmoothL1_1b::operator=(const ProblemL2_1b &source) {

	assert(this != &source);

	//printf("hello\n");

	ProblemSmoothL1_1b tmp(source);

	FreeBase();
	Free();

	dim = tmp.dim;
	dim1body = tmp.dim1body;

	SetLsingle(tmp.lsingle);

	MSA = tmp.MSA;

	AllocateBase();
	Allocate();

	memcpy(we, tmp.we, MSA->ncol * MSA->ncol * sizeof(double));

	memcpy(ea, tmp.ea, MSA->nrow * MSA->NAA * MSA->ncol * sizeof(double));
	memcpy(pa, tmp.pa, MSA->nrow * MSA->NAA * MSA->ncol * sizeof(double));
	memcpy(lpa, tmp.lpa, MSA->nrow * MSA->ncol * sizeof(double));

	return *this;

}

double ProblemSmoothL1_1b::Reg_f(const double *x) {

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

	return 0;

}

double ProblemSmoothL1_1b::Reg_fdf(const double *x, double *g) {

	//printf("Reg_fdf() call: eps= %f ...\n", eps);

	double reg = 0.0;

	/* regularize h */
#if defined(_OPENMP)
#pragma omp parallel for reduction (+:reg)
#endif
	for (size_t v = 0; v < dim1body; v++) {
//		double r1 = 1.0 + exp(-eps * x[v]);
//		double r2 = 1.0 + exp(eps * x[v]);
//		reg += lsingle / eps * (log(r1) + log(r2));
//		g[v] += lsingle * (1.0 / r1 - 1.0 / r2);

		reg += lsingle * sqrt(x[v] * x[v] + eps);
		g[v] += lsingle * x[v] / sqrt(x[v] * x[v] + eps);
	}

	return reg;

}

void ProblemSmoothL1_1b::SetSmoothingFactor(double eps_) {

	eps = eps_;

}
