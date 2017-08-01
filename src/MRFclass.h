/*
 * MRFclass.h
 *
 *  Created on: Jul 20, 2017
 *      Author: ivan
 */

#ifndef MRFCLASS_H_
#define MRFCLASS_H_

#include "MSAclass.h"

/* TODO: specs
 * - score random sequence
 * - score pair of positions
 * - score a patch
 */

/*
 * TODO: ??? How to deal with gaps in the MSA and MRF ???
 *
 * TODO: ??? What to do with masked edges ???
 */

class MRFclass {

	friend class Minimizer;
	friend class MRFprocessor;

private:

	size_t dim; /* MRF dimension (sequence length) */

	size_t dimh;
	size_t dimJ;

	double *h; /* dim * NAA */
	double *J; /* dim * dim * NAA * NAA */

	bool *we; /* for masking */

	void Allocate();
	void Free();

public:

	MRFclass();
	MRFclass(const MRFclass &source);
	MRFclass(double *h, double *J, size_t dim);
	MRFclass(double *h, double *J, bool *we, size_t dim);
	MRFclass(const std::string &name); /* read MRF from file */
	~MRFclass();

	MRFclass& operator=(const MRFclass &source);

	size_t GetDim() const;

	void Save(const std::string &name) const;

	/*
	 * score the native MSA (MSA.ncol == MRF.dim)
	 */

	/* contacts indices match reference sequence */
	double GetPairEnergies(const MSAclass &MSA,
			const std::vector<std::pair<size_t, size_t> > &contacts) const;

	/* msa[] cleaned and converted to 0..20 characters;
	 * contacts indices should match msa[] */
	double GetPairEnergies(const unsigned char *msa, size_t nrow,
			const std::vector<std::pair<size_t, size_t> > &contacts) const;

	/*
	 * scoring of a nonnative MSA
	 */

//	void SetNonnativeMSA(const MSAclass &MSA);
	/* pair energy for two characters */
//	double GetPairEnergy(unsigned char a, unsigned char b, size_t imrf,
//			size_t jmrf) const;
	/* intra, inter and total energies for a MSA patch */
//	double GetPatchIntraEnergy(const MSAclass &MSA, size_t msabeg,
//			size_t mrfbeg, size_t size) const;
//	double GetPatchInterEnergy(const MSAclass &MSA, size_t msabeg,
//			size_t mrfbeg, size_t size) const;
//	double GetPatchTotalEnergy(const MSAclass &MSA, size_t msabeg,
//			size_t mrfbeg, size_t size);
};

#endif /* MRFCLASS_H_ */
