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

	dimh = dim * MSAclass::NAA;
	dimJ = dimh * dimh;

	Allocate();

	memcpy(h, h_, dimh * sizeof(double));
	memcpy(J, J_, dimJ * sizeof(double));

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
