/*
 * MRFclass.cpp
 *
 *  Created on: Jul 20, 2017
 *      Author: ivan
 */

#include <cstring>
#include <cassert>

#include "MRFclass.h"

MRFclass::MRFclass() :
		dim(0), dimh(0), dimJ(0), h(NULL), J(NULL) {

	/* */

}

MRFclass::MRFclass(const MRFclass &s) :
		dim(s.dim), dimh(s.dimh), dimJ(s.dimJ), h(NULL), J(NULL) {

	Allocate();

	memcpy(h, s.h, dimh * sizeof(double));
	memcpy(J, s.J, dimJ * sizeof(double));

}

MRFclass::MRFclass(double *h_, double *J_, size_t dim_) :
		dim(dim_), h(NULL), J(NULL) {

	size_t NAA = MSAclass::NAA;

	dimh = dim * NAA;
	dimJ = dimh * dimh;

	Allocate();

	for (size_t i = 0; i < dim; i++) {
		for (size_t a = 0; a < NAA; a++) {
			h[i * NAA + a] = h_[a * dim + i];
		}
	}

	for (size_t i = 0; i < dim; i++) {
		for (size_t a = 0; a < NAA; a++) {
			for (size_t j = 0; j < dim; j++) {
				for (size_t b = 0; b < NAA; b++) {
					J[(i * dim + j) * NAA * NAA + a * NAA + b] = J_[((a * dim
							+ i) * NAA + b) * dim + j];
				}
			}
		}
	}

}

MRFclass::MRFclass(const std::string &name) :
		dim(0), dimh(0), dimJ(0), h(NULL), J(NULL) {

	FILE *F = fopen(name.c_str(), "r");
	if (F == NULL) {
		printf("Error: cannot open '%s' file to save MRF\n", name.c_str());
		exit(1);
	}

	size_t NAA = MSAclass::NAA;

	/* dimension (seq length) */
	fscanf(F, "%lu\n", &dim);
	dimh = dim * NAA;
	dimJ = dimh * dimh;

	Allocate();

	/* local fields */
	for (size_t i = 0; i < dim; i++) {
		for (size_t a = 0; a < NAA; a++) {
			fscanf(F, "%lf ", h + (i * NAA + a));
		}
	}

	/* couplings */
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			double *Jp = J + i * dim + j;
			for (size_t aa = 0; aa < NAA * NAA; aa++) {
				fscanf(F, "%lf ", Jp++);
			}
		}
	}

	fclose(F);

}

MRFclass& MRFclass::operator=(const MRFclass &source) {

	assert(this != &source); /* an attempt to assign Residue to itself */

	Free();

	dim = source.dim;
	dimh = source.dimh;
	dimJ = source.dimJ;

	Allocate();

	memcpy(h, source.h, dimh * sizeof(double));
	memcpy(J, source.J, dimJ * sizeof(double));

	return *this;

}

MRFclass::~MRFclass() {

	Free();

}

void MRFclass::Allocate() {

	h = (double*) malloc(dimh * sizeof(double));
	J = (double*) malloc(dimJ * sizeof(double));

}

void MRFclass::Free() {

	free(h);
	free(J);

}

size_t MRFclass::GetDim() const {

	return dim;

}

void MRFclass::Save(const std::string &name) const {

	FILE *F = fopen(name.c_str(), "w");
	if (F == NULL) {
		printf("Error: cannot open '%s' file to save MRF\n", name.c_str());
		exit(1);
	}

	size_t NAA = MSAclass::NAA;

	/* dimension (seq length) */
	fprintf(F, "%ld\n", dim);

	/* local fields */
	for (size_t i = 0; i < dim; i++) {
		for (size_t a = 0; a < NAA; a++) {
			fprintf(F, "%.5e ", h[i * NAA + a]);
		}
		fprintf(F, "\n");
	}

	/* couplings */
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			double *Jp = J + i * dim + j;
			for (size_t aa = 0; aa < NAA * NAA; aa++) {
				fprintf(F, "%.5e ", *Jp++);
			}
			fprintf(F, "\n");
		}
	}

	fclose(F);

}
