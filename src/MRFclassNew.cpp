/*
 * MRFclassNew.cpp
 *
 *  Created on: Aug 24, 2018
 *      Author: aivan
 */

#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cstdio>
#include <cmath>

#include "MRFclassNew.h"

MRFclassNew::MRFclassNew() :
		dim(0), nvar1b(0), nvar2b(0), x(NULL), MSA(NULL) {

	/* */

}

MRFclassNew::MRFclassNew(const MSAclass &MSA_) :
		dim(MSA_.GetNcol()), nvar1b(0), nvar2b(0), x(NULL), MSA(&MSA_) {

	nvar1b = dim * MSA_.NAA;
	nvar2b = dim * (dim - 1) / 2 * MSA_.NAA * MSA_.NAA;

	Allocate();

}

MRFclassNew::MRFclassNew(const MRFclassNew &source) :
		dim(source.dim), nvar1b(source.nvar1b), nvar2b(source.nvar2b), x(NULL), MSA(
				source.MSA) {

	Allocate();
	memcpy(x, source.x, (nvar1b + nvar2b) * sizeof(double));

}

MRFclassNew::MRFclassNew(const std::string &name) :
		dim(0), nvar1b(0), nvar2b(0), x(NULL), MSA(NULL) {

	FILE *F = fopen(name.c_str(), "r");
	if (F == NULL) {
		printf("Error: cannot read MRF from file '%s'\n", name.c_str());
		exit(1);
	}

	size_t NAA = MSAclass::NAA;

	/* read dimensions */
	fscanf(F, "%lu %lu %lu\n", &dim, &nvar1b, &nvar2b);

	Allocate();

	/* read fields */
	for (size_t i = 0; i < dim; i++) {
		fscanf(F, "%*s\n");
		for (size_t a = 0; a < NAA; a++) {
//			fscanf(F, "%lf ", x + (i * NAA + a));
			fscanf(F, "%lf ", x + (a * dim + i));
		}
	}

	/* couplings */
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = i + 1; j < dim; j++) {
			double *Jp = x + nvar1b + (i * dim + j) * NAA * NAA;
			fscanf(F, "%*s\n");
			for (size_t aa = 0; aa < NAA * NAA; aa++) {
				fscanf(F, "%lf ", Jp++);
			}
		}
	}

	fclose(F);

}

MRFclassNew::~MRFclassNew() {

	Free();

}

MRFclassNew& MRFclassNew::operator=(const MRFclassNew &source) {

	assert(this != &source);

	Free();

	dim = source.dim;
	nvar1b = source.nvar1b;
	nvar2b = source.nvar2b;

	MSA = source.MSA;

	Allocate();

	memcpy(x, source.x, (nvar1b + nvar2b) * sizeof(double));

	return *this;

}

void MRFclassNew::Allocate() {

	if (MSA != NULL) {
		x = (double*) calloc(nvar1b + nvar2b, sizeof(double));
	}

}

void MRFclassNew::Free() {

	free(x);

}

size_t MRFclassNew::GetDim() const {

	return dim;

}

double* MRFclassNew::GetX() const {

	return x;

}

double MRFclassNew::FNorm(const double *mat, size_t dim) {

	double norm = 0.0, mean = 0.0;

	for (size_t aa = 0; aa < dim; aa++) {
		mean += mat[aa];
	}
	mean /= dim * dim;

	/* last symbol (gap) is omitted */
	for (size_t a = 0; a < dim - 1; a++) {
		for (size_t b = 0; b < dim - 1; b++) {
			double m = mat[a * dim + b] - mean;
			norm += m * m;
		}
	}

	return sqrt(norm);

}

void MRFclassNew::FN(double **mtx) const {

	size_t NAA = MSAclass::NAA;

	/* initialize mtx[][] with zeroes */
	for (size_t i = 0; i < dim; i++) {
		memset(mtx[i], 0, dim * sizeof(double));
	}

	/* Frobenius norms of (NAA-1) x (NAA-1) submatrices */
	for (size_t k = 0; k < dim * (dim - 1) / 2; k++) {

		size_t i, j;
		To2D(k, i, j);

		double *w = x + nvar1b + k * NAA * NAA;
		mtx[i][j] = FNorm(w, NAA);
		mtx[j][i] = mtx[i][j];

	}

}

void MRFclassNew::APC(double **mtx) const {

	FN(mtx);
	APC(dim, mtx);

}

void MRFclassNew::APC(size_t dim, double **mtx) {

	/* row/col averages */
	double *means = (double*) calloc(dim, sizeof(double));
	double meansum = 0.0;
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			double w = mtx[i][j];
			means[j] += w / dim;
			meansum += w;
		}
	}
	meansum /= dim * dim;

	/* apply APC to mtx[][] */
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			mtx[i][j] -= (means[i] * means[j]) / meansum;
		}
		mtx[i][i] = 0.0;
	}

	/* free */
	free(means);

}

void MRFclassNew::Zscore(size_t dim, double **mtx) {

	/* calculate average */
	double av = 0.0;
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = i + 1; j < dim; j++) {
			av += mtx[i][j];
		}
	}
	av = 2.0 * av / dim / (dim - 1);

	/* calculate standard deviation */
	double std = 0.0;
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = i + 1; j < dim; j++) {
			double d = mtx[i][j] - av;
			std += d * d;
		}
	}
	std = sqrt(2.0 * std / dim / (dim - 1));

	/* update matrix */
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = i + 1; j < dim; j++) {
			double x = (mtx[i][j] - av) / std;
			mtx[i][j] = mtx[j][i] = x;
		}
	}

}

void MRFclassNew::SetMSA(const MSAclass &MSA_) {

	assert(MSA_.GetNcol() == dim);

	MSA = &MSA_;

}

void MRFclassNew::Save(const std::string &name) const {

	FILE *F = fopen(name.c_str(), "w");
	if (F == NULL) {
		printf("Error: cannot open '%s' file to save MRF\n", name.c_str());
		exit(1);
	}

	size_t NAA = MSAclass::NAA;

	/* save dimensions */
	fprintf(F, "%ld %ld %ld\n", dim, nvar1b, nvar2b);

	/* local fields */
	for (size_t i = 0; i < dim; i++) {
		fprintf(F, "V[%ld]\n", i);
		for (size_t a = 0; a < NAA; a++) {
//			fprintf(F, "%.5e ", x[i * NAA + a]);
			fprintf(F, "%.5e ", x[a * dim + i]);
		}
		fprintf(F, "\n");
	}

	/* couplings */
	for (size_t k = 0; k < dim * (dim - 1) / 2; k++) {

		size_t i, j;
		To2D(k, i, j);

		double *Jp = x + nvar1b + k * NAA * NAA;

		fprintf(F, "W[%ld][%ld]\n", i, j);
		for (size_t aa = 0; aa < NAA * NAA; aa++) {
			fprintf(F, "%.5e ", *Jp++);
		}
		fprintf(F, "\n");

	}

	fclose(F);

}

void MRFclassNew::To2D(size_t k, size_t &i, size_t &j) const {

	size_t n = MSA->GetNcol();

// https://stackoverflow.com/questions/27086195/linear-index-upper-triangular-matrix
	i = n - 2 - floor(sqrt(-8 * k + 4 * n * (n - 1) - 7) / 2.0 - 0.5);
	j = k + i + 1 - n * (n - 1) / 2 + (n - i) * ((n - i) - 1) / 2;

}
