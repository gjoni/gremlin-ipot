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

void ContactList::Print() {

	/* print caption */
	printf("%12s", "Score");
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
