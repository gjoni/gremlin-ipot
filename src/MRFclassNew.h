/*
 * MRFclassNew.h
 *
 *  Created on: Aug 24, 2018
 *      Author: aivan
 */

#ifndef MRFCLASSNEW_H_
#define MRFCLASSNEW_H_

#include "MSAclass.h"

class MRFclassNew {

private:

	/* MRF dimensions (cleaned sequence) */
	size_t dim;
	size_t nvar1b;
	size_t nvar2b;

	/* MRF parameters: both h[] and J[][] */
	double *x;

	const MSAclass *MSA;

	void Allocate();
	void Free();

	void To2D(size_t k, size_t &i, size_t &j) const;

	static double FNorm(const double *x, size_t dim);

public:

	MRFclassNew();
	MRFclassNew(const MSAclass &MSA);
	MRFclassNew(const MRFclassNew &source);

	/* TODO: check correctness of this constructor
	 *       (probably incorrect) */
	MRFclassNew(const std::string &name);

	~MRFclassNew();

	MRFclassNew& operator=(const MRFclassNew &source);

	size_t GetDim() const;

	double* GetX() const;

	void SetMSA(const MSAclass &MSA);

	void Save(const std::string &name) const;

	/*
	 * FUNCTIONS FOR MRF CONVERTION INTO CONTACT MAPS
	 * AND FOR MANIPULATION OF THE LATTER
	 */
	void FN(double **mtx) const;
	void APC(double **mtx) const;
//	void DI(double **mtx);

	static void APC(size_t dim, double **mtx);
	static void Zscore(size_t dim, double **mtx);

//	static void SaveMTX(size_t dim, double **mtx);
//	static void SaveMTX(const MSAclass &MSA, double **mtx);

	/*
	 * FUNCTIONS TO SCORE SEQUENCES
	 * BASED ON THE TRAINED MRF
	 */
//	double E(const std::string &seq);
	/* mode = 1 - chainA, 2 - chainB, 12 - interchain */
//	double E(const std::string &seq, size_t shift, int mode = 12);
};

#endif /* MRFCLASSNEW_H_ */
