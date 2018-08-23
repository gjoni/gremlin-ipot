/*
 * MRFclass.h
 *
 *  Created on: Jul 20, 2017
 *      Author: ivan
 */

#ifndef MRFCLASS_H_
#define MRFCLASS_H_

#include "MSAclass.h"

class MRFclass {

	friend class Minimizer;
	friend class MRFprocessor;

private:

	size_t dim; /* MRF dimension (sequence length) */

	size_t dimh;
	size_t dimJ;

	size_t len_ref;

	/* TODO: save/read these arrays */
	std::vector<size_t> mtx_to_a3m; /* size() = dim */
	std::vector<size_t> a3m_to_mtx; /* size() = len_ref */

	double *h; /* dim * NAA */
	double *J; /* dim * dim * NAA * NAA */

	bool *we; /* for masking */

	void Allocate();
	void Free();

public:

	MRFclass();
	MRFclass(const MRFclass &source);
	MRFclass(double *x, const MSAclass *MSA);
	MRFclass(double *x, bool *we, const MSAclass *MSA);
	MRFclass(const std::string &name); /* read MRF from file */
	~MRFclass();

	MRFclass& operator=(const MRFclass &source);

	size_t GetDim() const;

	void Save(const std::string &name) const;

	/*
	 * score the native MSA (MSA.ncol == MRF.dim)
	 */

	std::vector<double> ScoreMSA(const MSAclass &MSA);

	/* contacts indices match reference sequence */
	double GetPairEnergies(const MSAclass &MSA,
			const std::vector<std::pair<size_t, size_t> > &contacts) const;

	/* msa[] cleaned and converted to 0..20 characters;
	 * contacts indices should match msa[] */
	double GetPairEnergies(const unsigned char *msa, size_t nrow,
			const std::vector<std::pair<size_t, size_t> > &contacts) const;

};

#endif /* MRFCLASS_H_ */
