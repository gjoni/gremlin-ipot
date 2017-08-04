/*
 * EigenRRCE.cpp
 *
 *  Created on: Jun 1, 2017
 *      Author: ivan
 */

#include "EigenRRCE.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>

#include <cstdlib>

EigenRRCE::EigenRRCE(double **J_) {

	Allocate();

	/*
	 * (1) initialize J
	 */
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			J[i][j] = J_[i][j];
		}
	}

	/*
	 * (2) decompose
	 */
	Decompose();
	printf("# e[]     :");
	for (int i = 0; i < dim; i++) {
		printf("%7.3f", e[i]);
	}
	printf("\n");

	for (int i = 0; i < dim; i++) {
		printf("# ev[%2d][]:", i);
		for (int j = 0; j < 20; j++) {
			printf("%7.3f", ev[i][j]);
		}
		printf("\n");
	}

}

EigenRRCE::EigenRRCE() {

	Allocate();

}

EigenRRCE::~EigenRRCE() {

	Free();

}

void EigenRRCE::Allocate() {

	J = (double**) malloc(dim * sizeof(double*));
	ev = (double**) malloc(dim * sizeof(double*));
	for (int i = 0; i < dim; i++) {
		J[i] = (double*) malloc(dim * sizeof(double));
		ev[i] = (double*) malloc(dim * sizeof(double));
	}

	e = (double*) malloc(dim * sizeof(double));

}

void EigenRRCE::Free() {

	free(e);
	for (int i = 0; i < dim; i++) {
		free(J[i]);
		free(ev[i]);
	}
	free(J);
	free(ev);

}

void EigenRRCE::Decompose() {

	gsl_matrix *m = gsl_matrix_alloc(dim, dim);
	gsl_matrix *evec = gsl_matrix_alloc(dim, dim);
	gsl_vector *eval = gsl_vector_alloc(dim);

	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			gsl_matrix_set(m, i, j, J[i][j]);
		}
	}

	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(dim);
	gsl_eigen_symmv(m, eval, evec, w);
	gsl_eigen_symmv_free(w);
	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_DESC);

	for (int i = 0; i < dim; i++) {
		e[i] = eval->data[i];
		for (int j = 0; j < dim; j++) {
			ev[i][j] = gsl_matrix_get(evec, j, i);
		}
	}

	gsl_matrix_free(m);
	gsl_matrix_free(evec);
	gsl_vector_free(eval);

}

double EigenRRCE::GetEigenEnergy(const int a, const int b, const int i) {

	double E = 0.0;

	if (i < 0) {
		for (int k = -i; k <= dim; k++) {
			E += e[k - 1] * ev[k - 1][a] * ev[k][b];
		}
	} else {
		E = e[i] * ev[i][a] * ev[i][b];
	}

	return E;

}

double EigenRRCE::GetReconstructionError(const int i) {

	double dr = 0.0, rr = 0.0;

	for (int a = 0; a < dim; a++) {
		for (int b = a; a < dim; a++) {
			double Eab = 0.0;
			for (int k = 0; k < i; k++) {
				Eab += e[k] * ev[k][a] * ev[k][b];
			}
			double d = Eab - J[a][b];
			double r = Eab + J[a][b];
			dr += d * d;
			rr += r * r;
		}
	}

	double err = sqrt(dr) / sqrt(rr);

	return err;
}

double EigenRRCE::GetReconstructionCorrel(const int i) {

	const int N = dim * (dim + 1) / 2;
	double *Eoriginal = (double*) malloc(N * sizeof(double));
	double *Ereconstr = (double*) malloc(N * sizeof(double));

	int idx = 0;

	for (int a = 0; a < dim; a++) {
		for (int b = a; b < dim; b++) {
			Eoriginal[idx] = J[a][b];
			Ereconstr[idx] = 0.0;
			for (int k = 0; k < i; k++) {
				Ereconstr[idx] += e[k] * ev[k][a] * ev[k][b];
			}
			//printf("%f %f\n", Eoriginal[idx], Ereconstr[idx]);
			idx++;
		}
	}

	double correl = gsl_stats_correlation(Eoriginal, 1, Ereconstr, 1, N);

	free(Eoriginal);
	free(Ereconstr);

	return correl;

}

double EigenRRCE::GetEigenvalue(const int i) {

	return e[i - 1];

}

void EigenRRCE::GetEigenmatrix(int i, double **em) {

	for (int a = 0; a < dim; a++) {
		for (int b = 0; b < dim; b++) {
			em[a][b] = ev[i - 1][a] * ev[i - 1][b];
		}
	}

}
