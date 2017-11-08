/*
 * ContactList.cpp
 *
 *  Created on: Nov 7, 2017
 *      Author: ivan
 */

#include "ContactList.h"

ContactList::ContactList(size_t _dim) :
		dim(_dim) {

	for (size_t i = 0; i < dim; i++) {
		for (size_t j = i + 1; j < dim; j++) {
			contact.push_back( { i, j, { }, 0.0 });
		}
	}

}

ContactList::ContactList(const EDGE_LIST &EL, size_t _dim) :
		dim(_dim) {

	/* loop over all lists */
	for (EDGE_LIST::const_iterator it = EL.begin(); it != EL.end(); it++) {

		/* loop over contacts within the current list */
		for (const auto &c : it->second) {
			contact.push_back( { c.first, c.second, { }, 0.0 });
		}

	}

}

ContactList::~ContactList() {

}

void ContactList::AddTerm(const std::pair<std::string, double> &t) {

	term.push_back(t);

}

void ContactList::AddFeature(const std::string &name, double **mtx) {

	fname.push_back(name);

	for (auto &c : contact) {
		double s = mtx[c.i][c.j];
		c.feature.push_back(s);
		c.score = s;
	}

	/* initialize total score with the first feature */
//	if (fname.size() == 1) {
//		for (auto &c : contact) {
//			c.score = c.feature.front();
//		}
//
//	}
}

void ContactList::Print() const {

	/* print caption */
	printf("#%9s%10s%12s", "Residue1", "Residue2", "Score");
	for (auto &str : fname) {
		printf("%12s", str.c_str());
	}
	for (auto &t : term) {
		printf("%12s", t.first.c_str());
	}
	printf("\n");

	/* printf features and terms */
	char pattern[] = "%12.5f";
	for (auto &c : contact) {
		printf("%10lu%10lu", c.i, c.j);
		printf(pattern, c.score);
		for (auto &f : c.feature) {
			printf(pattern, f);
		}
		for (auto &t : term) {
			printf(pattern, t.second);
		}
		printf("\n");
	}

}

void ContactList::Print(const MSAclass &MSA) const {

	/* print caption */
	printf("#%9s%10s%12s", "Residue1", "Residue2", "Score");
	for (auto &str : fname) {
		printf("%12s", str.c_str());
	}
	for (auto &t : term) {
		printf("%12s", t.first.c_str());
	}
	printf("\n");

	/* printf features and terms */
	char pattern[] = "%12.5f";
	for (auto &c : contact) {
		size_t i = MSA.GetA3MIdx(c.i);
		size_t j = MSA.GetA3MIdx(c.j);
		char a = MSA.GetA3Mres(0, i);
		char b = MSA.GetA3Mres(0, j);
		printf("%8lu%2c%8lu%2c", i, a, j, b);
		printf(pattern, c.score);
		for (auto &f : c.feature) {
			printf(pattern, f);
		}
		for (auto &t : term) {
			printf(pattern, t.second);
		}
		printf("\n");
	}

}

void ContactList::SaveMTX(const char *name) const {

	FILE *F = fopen(name, "w");
	if (F == NULL) {
		printf("Error: cannot open file to save MTX '%s'\n", name);
		exit(1);
	}

	double **mtx = (double**) malloc(dim * sizeof(double*));
	for (size_t i = 0; i < dim; i++) {
		mtx[i] = (double*) calloc(dim, sizeof(double));
	}

	for (auto &c : contact) {
		mtx[c.i][c.j] = c.score;
		mtx[c.j][c.i] = c.score;
	}

	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim - 1; j++) {
			fprintf(F, "%.5e ", mtx[i][j]);
		}
		fprintf(F, "%.5e\n", mtx[i][dim - 1]);
	}
	fclose(F);

	for (size_t i = 0; i < dim; i++) {
		free(mtx[i]);
	}
	free(mtx);

}

void ContactList::SaveMTX(const char *name, const MSAclass &MSA) const {

	FILE *F = fopen(name, "w");
	if (F == NULL) {
		printf("Error: cannot open file to save MTX '%s'\n", name);
		exit(1);
	}

	size_t len = MSA.GetLen();

	double **mtx = (double**) malloc(len * sizeof(double*));
	for (size_t i = 0; i < len; i++) {
		mtx[i] = (double*) calloc(len, sizeof(double));
	}

	for (auto &c : contact) {
		size_t i = MSA.GetA3MIdx(c.i);
		size_t j = MSA.GetA3MIdx(c.j);
		mtx[i][j] = c.score;
		mtx[j][i] = c.score;
	}

	for (size_t i = 0; i < len; i++) {
		for (size_t j = 0; j < len - 1; j++) {
			fprintf(F, "%.5e ", mtx[i][j]);
		}
		fprintf(F, "%.5e\n", mtx[i][len - 1]);
	}
	fclose(F);

	for (size_t i = 0; i < len; i++) {
		free(mtx[i]);
	}
	free(mtx);

}

EDGE_LIST ReadEdges(const char *name) {

	EDGE_LIST L;

	FILE *F = fopen(name, "r");
	if (F == NULL) {
		printf("Error: cannot open file for reading '%s'\n", name);
		exit(1);
	}

	size_t a, b;
	const size_t SIZE = 256;
	char buf[SIZE];
	while (fscanf(F, "%lu %lu %s\n", &a, &b, buf) == 3) {
		L[buf].push_back(std::make_pair(a, b));
	}
	fclose(F);

	return L;

}
