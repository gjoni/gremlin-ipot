/*
 * MRFprocessor.cpp
 *
 *  Created on: Jul 21, 2017
 *      Author: ivan
 */

#include <cstring>
#include <cmath>
#include <cassert>

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

	/* TODO: subtracting means gives almost no effect */

	for (size_t a = 0; a < dim - 1; a++) {
		for (size_t b = 0; b < dim - 1; b++) {
			double m = mat[a * dim + b] - mean;
			norm += m * m;
		}
	}

	return sqrt(norm);

}

void MRFprocessor::DI(const MRFclass &MRF, const MSAclass &MSA, double **mtx) {

	assert(MRF.GetDim() == MSA.GetNcol()); /* MRF/MSA size mismatch */

	size_t ncol = MRF.GetDim();
	size_t NAA = MSAclass::NAA;

	/* initialize mtx with zeroes */
	for (size_t i = 0; i < ncol; i++) {
		memset(mtx[i], 0, ncol * sizeof(double));
	}

	/*
	 * 1-site probabilities
	 */
	double **Pi = (double**) malloc(ncol * sizeof(double*));
	for (size_t i = 0; i < ncol; i++) {

		Pi[i] = (double*) malloc((NAA - 1) * sizeof(double));

		double *hi = MRF.h + i * NAA;
		double Z = 0.0;
		for (size_t a = 0; a < NAA - 1; a++) {
			Pi[i][a] = exp(hi[a]);
			Z += Pi[i][a];
		}

		for (size_t a = 0; a < NAA - 1; a++) {
			Pi[i][a] /= Z;
		}

	}

	/*
	 * 2-site probabilities
	 */
	for (size_t i = 0; i < ncol; i++) {
		double *hi = MRF.h + i * NAA;
		for (size_t j = i + 1; j < ncol; j++) {
			double Zij = 0.0;
			double *Jij = MRF.J + (i * ncol + j) * NAA * NAA;
			double *hj = MRF.h + j * NAA;

			double **eab = (double**) malloc((NAA - 1) * sizeof(double*));
			for (size_t a = 0; a < NAA - 1; a++) {
				eab[a] = (double*) malloc((NAA - 1) * sizeof(double));
				for (size_t b = 0; b < NAA - 1; b++) {
					eab[a][b] = Jij[a * NAA + b] + hi[a] + hj[b];
					Zij += exp(eab[a][b]);
				}
			}

//			for (size_t a = 0; a < NAA - 1; a++) {
//				for (size_t b = 0; b < NAA - 1; b++) {
//					double pab = exp(eab[a][b]) / Zij;
//					mtx[i][j] += pab * log(pab / Pi[i][a] / Pi[j][b]);
//				}
//				free(eab[a]);
//			}
//			free(eab);

			for (size_t a = 0; a < NAA - 1; a++) {
				double fa = MSA.GetFi(i, a);
				for (size_t b = 0; b < NAA - 1; b++) {
					double pab = exp(eab[a][b]) / Zij;
					double fb = MSA.GetFi(j, b);
					mtx[i][j] += pab * log(pab / fa / fb);
				}
				free(eab[a]);
			}
			free(eab);

			mtx[j][i] = mtx[i][j];

		}
	}

	/*
	 * free
	 */
	for (size_t i = 0; i < ncol; i++) {
		free(Pi[i]);
	}
	free(Pi);

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
		mtx[i][i] = 0.0;
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

void MRFprocessor::Zscore(size_t dim, double **mtx) {

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
