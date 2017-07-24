/*
 * MRFprocessor.h
 *
 *  Created on: Jul 21, 2017
 *      Author: ivan
 */

#ifndef MRFPROCESSOR_H_
#define MRFPROCESSOR_H_

#include <vector>
#include <tuple>

#include "MRFclass.h"

class MRFprocessor {
private:
	MRFprocessor();
	~MRFprocessor();

	static double FNorm(const double *x, size_t dim);

public:

	struct Contact {
		size_t a;
		size_t b;
		double score;
	};

	/* TODO:
	 * 1) APC correction
	 * 2) selection of top contacts
	 * 3) conversion into Rosetta constraints */

	static void APC(const MRFclass &MRF, double **mtx);
//	static void APC(const MRFclass &MRF, std::vector<Contact> &contacts);

//	static std::vector<Contact> GetTopContacts(const double **mtx);

	/* TODO: 'energy' score */
//	static double GetScore(const MRFclass &MRF, const MSAclass &MSA,
//			const std::vector<std::pair<size_t, size_t> > &contacts);
	/* TODO: 'contact' score */
//	static double GetScore();
};

#endif /* MRFPROCESSOR_H_ */
