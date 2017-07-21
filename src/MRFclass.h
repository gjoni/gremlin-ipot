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

private:
	size_t dim; /* MRF dimension (sequence length) */

	size_t dimh;
	size_t dimJ;

	double *h; /* NAA * dim */
	double *J; /* NAA * NAA * dim * dim */

	void Allocate();
	void Free();

public:
	MRFclass();
	MRFclass(const MRFclass &source);
	MRFclass(double *h, double *J, size_t dim);
	~MRFclass();

	MRFclass& operator=(const MRFclass &source);

	/* pair energy for two characters */
//	double GetPairEnergy(unsigned char a, unsigned char b, size_t imrf,
//			size_t jmrf) const;

	/* pair energy for two MSA positions */
//	double GetPairEnergy(const MSAclass &MSA, size_t imsa, size_t jmsa,
//			size_t imrf, size_t jmrf) const;

	/* intra, inter and total energies for a MSA patch */
//	double GetPatchIntraEnergy(const MSAclass &MSA, size_t msabeg,
//			size_t mrfbeg, size_t size) const;
//	double GetPatchInterEnergy(const MSAclass &MSA, size_t msabeg,
//			size_t mrfbeg, size_t size) const;
//	double GetPatchTotalEnergy(const MSAclass &MSA, size_t msabeg,
//			size_t mrfbeg, size_t size);

	/* TODO: save to / read from file */

};

#endif /* MRFCLASS_H_ */
