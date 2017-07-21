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

#include <algorithm>

/* A B C D E F G H I J K L M N O P Q R S T U V W X Y Z */
const unsigned char MSAclass::AMINO_INDICES[26] = { 0, 20, 4, 3, 6, 13, 7, 8, 9,
		20, 11, 10, 12, 2, 20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 20 };

/* 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 *
 * A R N D C Q E G H I  L  K  M  F  P  S  T  W  Y  V  - */
const unsigned char MSAclass::CHAR_INDICES[NAA] = { 65, 82, 78, 68, 67, 81, 69,
		71, 72, 73, 76, 75, 77, 70, 80, 83, 84, 87, 89, 86, 45 };

unsigned char MSAclass::aatoi(unsigned char aa) {

	if (!isalpha(aa)) {
		return 20;
	}

	aa = toupper(aa);
	if (aa < 65 || aa > 90) {
		return 20;
	}

	return MSAclass::AMINO_INDICES[aa - 65];

}

unsigned char MSAclass::itoaa(unsigned char i) {

	if (i < 0 || i > 20) {
		return '-';
	} else {
		return MSAclass::CHAR_INDICES[i];
	}

}

void MSAclass::aatoi(unsigned char *str, size_t len) {

	for (size_t i = 0; i < len; i++) {
		*str = aatoi(*str);
		str++;
	}

}

MSAclass::MSAclass() :
		len_ref(0), nrow(0), ncol(0) {

	/* nothing to be done */

}

MSAclass::MSAclass(const char *name) :
		len_ref(0), nrow(0), ncol(0) {

	FILE *F = fopen(name, "r");
	if (F == NULL) {
		printf("Error: cannot open %s file\n", name);
		exit(1);
	}

	const size_t STRLEN = 65536;

	char buf[STRLEN];

	fgets(buf, STRLEN, F);
	if (buf[0] != '>') {
		printf("Error: MSA '%s' is not in A3M/FASTA format\n", name);
		exit(1);
	}

	TrimRight(buf);
	a3m.push_back(std::make_pair(std::string(buf), std::string("")));

	while (fgets(buf, STRLEN, F)) {

		TrimRight(buf);

		if (buf[0] == '>') /* create a new entry */{
			a3m.push_back(std::make_pair(std::string(buf), std::string("")));
		} else /* update existing */{
			a3m.back().second += std::string(buf);
		}
	}
	fclose(F);

	/* clean lowercase */
	len_ref = a3m[0].second.size(); /* reference sequence length */
	for (size_t i = 1; i < a3m.size(); i++) {
		size_t len = CleanLowercase(a3m[i].second);
		if (len != len_ref) {
			printf("Error: sequence length mismatch at a3m[%ld]\n", i + 1);
			exit(1);
		}
	}

	/* create row & col maps (no gaps will be cleaned here) */
	nrow = a3m.size();
	row_map.resize(nrow);
	std::generate(row_map.begin(), row_map.end(),
			[] {static int i {0}; return i++;});

	ncol = len_ref;
	col_map.resize(ncol);
	std::generate(col_map.begin(), col_map.end(),
			[] {static int i {0}; return i++;});

	a3m_to_msa = col_map;

	/* populate msa */
	Allocate();

}

MSAclass::MSAclass(const MSAclass &source) :
		a3m(source.a3m), len_ref(source.len_ref), row_map(source.row_map), col_map(
				source.col_map), a3m_to_msa(source.a3m_to_msa), nrow(
				source.nrow), ncol(source.ncol) {

	Allocate();

}

MSAclass & MSAclass::operator =(const MSAclass & source) {

	assert(this != &source); /* an attempt to assign Residue to itself */

	Free();

	a3m = source.a3m;
	len_ref = source.len_ref;

	row_map = source.row_map;
	col_map = source.col_map;

	a3m_to_msa = source.a3m_to_msa;

	nrow = source.nrow;
	ncol = source.ncol;

	Allocate();

	return *this;

}

void MSAclass::Allocate() {

	/* */

}

void MSAclass::Free() {

	/* */

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

size_t MSAclass::CleanLowercase(std::string &str) {

	size_t pos = 0;

	for (size_t i = 0; i < str.size(); i++) {
		if (!islower(str[i])) {
			str[pos] = str[i];
			pos++;
		}
	}
	str.resize(pos);

	return pos;

}

char MSAclass::GetA3Mres(size_t i, size_t j) const {

	return a3m[i].second[j];

}

size_t MSAclass::GetNrow() const {

	return row_map.size();

}

size_t MSAclass::GetNcol() const {

	return col_map.size();

}

size_t MSAclass::GetMsaIdx(size_t idx) const {

	if (idx < 0 || idx >= len_ref) {
		printf("Error: A3M index (%ld) out of range\n", idx);
		exit(1);
	}

	return a3m_to_msa[idx];

}

unsigned char * MSAclass::GetMsa() const {

	unsigned char *msa_ = (unsigned char*) malloc(
			nrow * ncol * sizeof(unsigned char));

	size_t idx = 0;

	for (size_t i = 0; i < row_map.size(); i++) {

		size_t irow = row_map[i];
		const std::string &seq = a3m[irow].second;

		for (size_t j = 0; j < col_map.size(); j++) {
			size_t icol = col_map[j];
			msa_[idx] = seq[icol];
			idx++;
		}
	}

	return msa_;

}

size_t MSAclass::CleanRows(double gaps_frac) {

	row_map.clear();

	for (size_t i = 0; i < a3m.size(); i++) {

		std::string &seq = a3m[i].second;
		size_t ngaps = std::count(seq.begin(), seq.end(), '-');
		if ((double) ngaps / seq.size() < gaps_frac) {
			row_map.push_back(i);
		}

	}

	nrow = row_map.size();

	return nrow;

}

size_t MSAclass::CleanCols(double cols_frac) {

	col_map.clear();
	a3m_to_msa.resize(len_ref);

	/*
	 * calculate numbers of gaps for each column
	 * (only for sequences in row_map[])
	 */
	size_t *counts = (size_t*) malloc(len_ref * sizeof(size_t));
	memset(counts, 0, len_ref * sizeof(size_t));

	for (size_t i = 0; i < row_map.size(); i++) {

		std::string &seq = a3m[row_map[i]].second;
		for (size_t j = 0; j < len_ref; j++) {
			if (seq[j] == '-') {
				counts[j]++;
			}
		}

	}

	/*
	 * populate col_map[] and a3m_to_msa[]
	 */
	size_t count = 0;
	for (size_t i = 0; i < len_ref; i++) {

		if ((double) counts[i] / nrow < cols_frac) {
			col_map.push_back(i);
			a3m_to_msa[i] = count;
			count++;
		} else {
			a3m_to_msa[i] = SIZE_MAX;
		}

	}

	free(counts);

	ncol = col_map.size();

	return ncol;

}

void MSAclass::CleanMsa(double rgaps, double cgaps) {

	CleanRows(rgaps);
	CleanCols(cgaps);

}

void MSAclass::SaveMSA(const std::string &name) const {

	FILE *F = fopen(name.c_str(), "w");

	if (F == NULL) {
		printf("Error: cannot save MSA '%s'\n", name.c_str());
		exit(1);
	}

	for (size_t i = 0; i < row_map.size(); i++) {

		size_t irow = row_map[i];
		const std::string &seq = a3m[irow].second;

		for (size_t j = 0; j < col_map.size(); j++) {
			size_t icol = col_map[j];
			fprintf(F, "%c", seq[icol]);
		}

		fprintf(F, "\n");

	}

	fclose(F);

}

void MSAclass::PrintMSA() const {

	for (size_t i = 0; i < row_map.size(); i++) {

		size_t irow = row_map[i];
		const std::string &seq = a3m[irow].second;

		for (size_t j = 0; j < col_map.size(); j++) {
			size_t icol = col_map[j];
			printf("%c", seq[icol]);
		}

		printf("\n");

	}

}
