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
	// TODO Auto-generated constructor stub

}

MRFprocessor::~MRFprocessor() {
	// TODO Auto-generated destructor stub
}

double MRFprocessor::FNorm(const double *mat, size_t dim) {

	double norm = 0.0, mean = 0.0;

	for (size_t aa = 0; aa < dim; aa++) {
		mean += mat[aa];
	}
	mean /= dim * dim;

	for (size_t a = 0; a < dim; a++) {
		for (size_t b = 0; b < dim; b++) {
			double m = mat[a * dim + b] - mean;
			norm += m * m;
		}
	}

	return sqrt(norm);

}

void MRFprocessor::APC(const MRFclass &MRF, double **mtx) {

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
