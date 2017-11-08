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
#include <map>

#include "MRFclass.h"

class MRFprocessor {
private:

//	const MSAclass *MSA;
//	const MRFclass *MRF;

	MRFprocessor();
	~MRFprocessor();

	static double FNorm(const double *x, size_t dim);

public:

	struct MTX {
		std::vector<double> mtx1d; /* 1d array [dim x dim] */
		size_t dim;
	};

	static void FN(const MRFclass &MRF, double **mtx);
	static void APC(const MRFclass &MRF, double **mtx);
	/* TODO: Zscore*/
	static void Zscore(size_t dim, double **mtx);

	/* TODO */
//	MRFprocessor(const MSAclass &MSA, const MRFclass &MRF);
	/* TODO:
	 * 1) selection of top contacts
	 * 2) scoring of contacts by 'contact' score
	 * 3) conversion into Rosetta constraints */

	/* Frobenius norms of the 20x20 submatrices */
	static void FN(const MRFclass &MRF, MTX &result);

	/* average product correction (APC) */
	static void APC(const MRFclass &MRF, MTX &result);

	/* APC for a submatrix - block-APC */
	static void BAPC(const MRFclass &MRF, MTX &result, size_t shift);

	/* save to file */
	static void SaveMTX(const MTX &result, const std::string &name);

	/* TODO: convert matrix to z-scores */
//	static void ToZscores(MTX &result);

	/* 'contact' score */
	static double GetScore(const MTX &result,
			const std::vector<std::pair<size_t, size_t> > &contacts);

};

#endif /* MRFPROCESSOR_H_ */
