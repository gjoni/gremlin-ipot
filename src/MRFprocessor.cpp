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

	/* Frobenius norms of 20x20 submatrices */
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
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			mtx[i][j] -= (means[i] * means[j]) / meansum;
		}
	}

	/* shift all elements by min */
	/*
	 double min = 1.0;
	 for (size_t i = 0; i < dim; i++) {
	 for (size_t j = i + 1; j < dim; j++) {
	 if (mtx[i][j] < min) {
	 min = mtx[i][j];
	 }
	 }
	 }

	 for (size_t i = 0; i < dim; i++) {
	 for (size_t j = 0; j < dim; j++) {
	 mtx[i][j] -= min;
	 }
	 mtx[i][i] = 0.0;
	 }
	 */

	/* free */
	free(means);

}

void MRFprocessor::BAPC(const MRFclass &MRF, MTX &result, size_t shift) {

	size_t dim = MRF.GetDim();
	result.dim = dim;

	result.mtx1d.resize(dim * dim);

	double **mtx = (double**) malloc(dim * sizeof(double*));
	for (size_t i = 0; i < dim; i++) {
		mtx[i] = (double*) malloc(dim * sizeof(double));
	}

	FN(MRF, mtx);

	/*
	 * row/col averages
	 */
	double *means_diag = (double*) calloc(dim, sizeof(double));
	double *means_offd = (double*) calloc(dim, sizeof(double));

	/* diagonal block A */
	double meansum_A = 0.0;
	for (size_t i = 0; i < shift; i++) {
		for (size_t j = 0; j < shift; j++) {
			double w = mtx[i][j];
			means_diag[j] += w / shift;
			meansum_A += w;
		}
	}
	meansum_A /= shift * shift;

	/* diagonal block B */
	double meansum_B = 0.0;
	for (size_t i = shift; i < dim; i++) {
		for (size_t j = shift; j < dim; j++) {
			double w = mtx[i][j];
			means_diag[j] += w / (dim - shift);
			meansum_B += w;
		}
	}
	meansum_B /= (dim - shift) * (dim - shift);

	/* off-diagonal block AB */
	double meansum_AB = 0.0;
	for (size_t i = 0; i < shift; i++) {
		for (size_t j = shift; j < dim; j++) {
			double w = mtx[i][j];
			means_offd[i] += w / (dim - shift);
			means_offd[j] += w / shift;
			meansum_AB += w;
		}
	}
	meansum_AB /= shift * (dim - shift);

	/*
	 * apply APC to mtx[][]
	 */

	/* diagonal block A */
	for (size_t i = 0; i < shift; i++) {
		for (size_t j = 0; j < shift; j++) {
			mtx[i][j] -= (means_diag[i] * means_diag[j]) / meansum_A;
		}
	}

	/* diagonal block B */
	for (size_t i = shift; i < dim; i++) {
		for (size_t j = shift; j < dim; j++) {
			mtx[i][j] -= (means_diag[i] * means_diag[j]) / meansum_B;
		}
	}

	/* off-diagonal block AB */
	for (size_t i = 0; i < shift; i++) {
		for (size_t j = shift; j < dim; j++) {
			mtx[i][j] -= (means_offd[i] * means_offd[j]) / meansum_AB;
			mtx[j][i] = mtx[i][j];
		}
	}

	/*
	 * copy mtx[][] into MRF
	 */
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			result.mtx1d[i * dim + j] = mtx[i][j];
		}
	}

	/*
	 * free
	 */
	for (size_t i = 0; i < dim; i++) {
		free(mtx[i]);
	}
	free(mtx);

	free(means_diag);
	free(means_offd);

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
