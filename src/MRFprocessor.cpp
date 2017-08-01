/*
 * MRFprocessor.cpp
 *
 *  Created on: Jul 21, 2017
 *      Author: ivan
 */

#include <cstring>
#include <cmath>

#include "MRFprocessor.h"
#include "MSAclass.h"

MRFprocessor::MRFprocessor() {
	/* */
}

MRFprocessor::~MRFprocessor() {
	/* */
}

double MRFprocessor::FNorm(const double *mat, size_t dim) {

	double norm = 0.0, mean = 0.0;

	for (size_t aa = 0; aa < dim; aa++) {
		mean += mat[aa];
	}
	mean /= dim * dim;

	/* TODO: subtracting of means gives almost no effect */

	for (size_t a = 0; a < dim - 1; a++) {
		for (size_t b = 0; b < dim - 1; b++) {
			double m = mat[a * dim + b] - mean;
			norm += m * m;
		}
	}

	return sqrt(norm);

}

void MRFprocessor::FN(const MRFclass &MRF, double **mtx) {

	size_t dim = MRF.GetDim();
	size_t NAA = MSAclass::NAA;

	/* initialize mtx with zeroes */
	for (size_t i = 0; i < dim; i++) {
		memset(mtx[i], 0, dim * sizeof(double));
	}

	/* Frobenius norms of 21x21 submatrices */
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			mtx[i][j] = FNorm(MRF.J + (i * dim + j) * NAA * NAA, NAA);
		}
	}

}

void MRFprocessor::APC(const MRFclass &MRF, double **mtx) {

	size_t dim = MRF.GetDim();

	FN(MRF, mtx);

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
	double min = 1.0;
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			mtx[i][j] -= (means[i] * means[j]) / meansum;
			if (mtx[i][j] < min) {
				min = mtx[i][j];
			}
		}
	}

	/* shift all elements by min */
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			mtx[i][j] -= min;
		}
		mtx[i][i] = 0.0;
	}

	/* free */
	free(means);

}

void MRFprocessor::APC(const MRFclass &MRF, MTX &result) {

	size_t dim = MRF.GetDim();
	result.dim = dim;

	result.mtx1d.resize(dim * dim);

	double **mtx = (double**) malloc(dim * sizeof(double*));
	for (size_t i = 0; i < dim; i++) {
		mtx[i] = (double*) malloc(dim * sizeof(double));
	}

	APC(MRF, mtx);

	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			result.mtx1d[i * dim + j] = mtx[i][j];
		}
	}

	for (size_t i = 0; i < dim; i++) {
		free(mtx[i]);
	}
	free(mtx);

}

void MRFprocessor::FN(const MRFclass &MRF, MTX &result) {

	size_t dim = MRF.GetDim();
	result.dim = dim;

	result.mtx1d.resize(dim * dim);

	double **mtx = (double**) malloc(dim * sizeof(double*));
	for (size_t i = 0; i < dim; i++) {
		mtx[i] = (double*) malloc(dim * sizeof(double));
	}

	FN(MRF, mtx);

	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			result.mtx1d[i * dim + j] = mtx[i][j];
		}
	}

	for (size_t i = 0; i < dim; i++) {
		free(mtx[i]);
	}
	free(mtx);

}

void MRFprocessor::SaveMTX(const MTX &result, const std::string &name) {

	FILE *F = fopen(name.c_str(), "w");
	if (F == NULL) {
		printf("Error: cannot open '%s' file to save MTX\n", name.c_str());
		exit(1);
	}

	for (size_t i = 0; i < result.dim; i++) {
		for (size_t j = 0; j < result.dim; j++) {
			fprintf(F, "%.5e ", result.mtx1d[i * result.dim + j]);
		}
		fprintf(F, "\n");
	}

	fclose(F);

}

double MRFprocessor::GetScore(const MTX &result,
		const std::vector<std::pair<size_t, size_t> > &contacts) {

	double E = 0.0;

	size_t dim = result.dim;

	for (const auto& c : contacts) {
		size_t i = c.first;
		size_t j = c.second;
		if (i < dim && j < dim && i >= 0 && j >= 0) {
			E += result.mtx1d[i * dim + j];
		}
	}

	return E;

}
