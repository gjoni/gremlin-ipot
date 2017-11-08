/*
 * ContactList.h
 *
 *  Created on: Nov 7, 2017
 *      Author: ivan
 */

#ifndef CONTACTLIST_H_
#define CONTACTLIST_H_

#include <unistd.h>

#include <vector>
#include <list>
#include <string>
#include <map>
#include <utility>

#include "MSAclass.h"

/*
 * basic structure to store information on a contact
 */
struct Contact {

	/* contacting pair */
	size_t i;
	size_t j;

	/* features/scores */
	std::vector<double> feature;

	/* total score */
	double score;

};

/*
 * compare Contacts based on the total score
 */
//bool operator <(const Contact &a, const Contact &b) {
//	return a.score < b.score;
//}
//
//bool operator >(const Contact &a, const Contact &b) {
//	return a.score > b.score;
//}
//
//bool operator ==(const Contact &a, const Contact &b) {
//	return a.score == b.score;
//}
/*
 * named lists of edges
 */
typedef std::map<std::string, std::vector<std::pair<size_t, size_t> > > EDGE_LIST;
EDGE_LIST ReadEdges(const char *name);

/*
 * container class to store Contacts
 */
class ContactList {
private:

	/* number of residues (cleaned) */
	size_t dim;

	/* common features and their names */
	std::vector<std::pair<std::string, double> > term;

	/* names of individual features */
	std::vector<std::string> fname;

	/* contacts are stored here */
	std::list<Contact> contact;

	/* leave undefined */
	ContactList();

public:

	ContactList(size_t dim);
	ContactList(const EDGE_LIST &EL, size_t dim);

	~ContactList();

	/*
	 * 'features' - are contact-specific
	 * 'terms'    - are shared between all contacts (e.g. Neff)
	 */

	/* add a new common term */
	void AddTerm(const std::pair<std::string, double> &t);

	/* adds a new feature for every contact[] */
	void AddFeature(const std::string &name, double **mtx);

	/*
	 * output
	 */

	void Print() const;
	void Print(const MSAclass &MSA) const;

	void SaveMTX(const char *name) const;
	void SaveMTX(const char *name, const MSAclass &MSA) const;

//	void Save();
//	void SaveMTX();

//	void AddContact(const Contact &c);
//	void AddContacts(const std::vector<Contact> &v);

};

#endif /* CONTACTLIST_H_ */
