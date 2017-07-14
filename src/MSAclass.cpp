/*
 * MSAclass.cpp
 *
 *  Created on: Jul 8, 2016
 *      Author: ivan
 */

#include "MSAclass.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>

#define STRLEN 65536

const unsigned char MSAclass::AMINO_INDICES[26] = {
		//  A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X   Y   Z
		0, 20, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2, 20, 14, 5, 1, 15, 16,
		20, 19, 17, 20, 18, 20 };

const unsigned char MSAclass::CHAR_INDICES[21] = {
		//  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
		//  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   -
		65, 82, 78, 68, 67, 81, 69, 71, 72, 73, 76, 75, 77, 70, 80, 83, 84, 87,
		89, 86, 45 };

MSAclass::MSAclass() :
		msa(NULL), nrow(0), ncol(0) {

	/* nothing to be done */

}

MSAclass::MSAclass(const char *name) :
		msa(NULL), nrow(0), ncol(0) {

	FILE *F = fopen(name, "r");
	if (F == NULL) {
		printf("Error: cannot open %s file\n", name);
		exit(1);
	}

	char seq[STRLEN];
	while (fgets(seq, STRLEN, F)) {
		nrow++;
		TrimRight(seq);
		size_t len = strlen(seq);
		ncol = len > ncol ? len : ncol;
	}

	rewind(F);

	Allocate();

	for (size_t i = 0; i < nrow; i++) {
		fgets(seq, STRLEN, F);
		for (size_t j = 0; j < ncol; j++) {
			msa[i][j] = seq[j];
		}
	}

}

MSAclass::MSAclass(const MSAclass &source) :
		msa(NULL), nrow(source.nrow), ncol(source.ncol) {

	Allocate();

	for (size_t i = 0; i < nrow; i++) {
		memcpy(msa[i], source.msa[i], ncol * sizeof(char));
	}

}

MSAclass & MSAclass::operator =(const MSAclass & source) {

	assert(this != &source); /* an attempt to assign Residue to itself */

	Free();

	nrow = source.nrow;
	ncol = source.ncol;

	Allocate();

	for (size_t i = 0; i < nrow; i++) {
		memcpy(msa[i], source.msa[i], ncol * sizeof(char));
	}

	return *this;

}

void MSAclass::Allocate() {

	msa = (char**) malloc(nrow * sizeof(char*));
	for (size_t i = 0; i < nrow; i++) {
		msa[i] = (char*) malloc(ncol * sizeof(char));
	}

}

void MSAclass::Free() {

	if (msa != NULL) {
		for (size_t i = 0; i < nrow; i++) {
			free(msa[i]);
		}
		free(msa);
	}
	msa = NULL;

}

MSAclass::~MSAclass() {

	Free();

}

void MSAclass::TrimRight(char *str) {

	/*
	 * code from CCMpred
	 */

	// find first non-whitespace character from right
	char *end = str + strlen(str) - 1;
	while (end > str && isspace(*end))
		end--;

	// add new null terminator
	*(end + 1) = 0;

}

char MSAclass::GetResidue(size_t i, size_t j) {

	return msa[i][j];

}

size_t MSAclass::GetNrow() {

	return nrow;

}

size_t MSAclass::GetNcol() {

	return ncol;

}

int MSAclass::aatoi(char aa) {

	if (!isalpha(aa)) {
		return 20;
	}

	aa = toupper(aa);
	if (aa < 65 || aa > 90) {
		return 20;
	}

	return MSAclass::AMINO_INDICES[aa - 65];

}
