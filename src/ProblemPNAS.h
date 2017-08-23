/*
 * ProblemPNAS.h
 *
 *  Created on: Aug 23, 2017
 *      Author: ivan
 */

#ifndef PROBLEMPNAS_H_
#define PROBLEMPNAS_H_

#include "ProblemBase.h"

class ProblemPNAS: public ProblemBase {

private:
	size_t nmodes; /* number of eigenmodes */

	double *em; /* 20 x 20 x nmodes - stores eigenmatrices */

	double lsingle;
	double lpair;

	size_t dim1body;
	size_t dim2body;

	void Allocate();
	void Free();

	static const double peak1[];
	static const double peak2[];
	static const double peak3[];

public:

	ProblemPNAS();
	ProblemPNAS(const MSAclass &MSA, size_t n);
	~ProblemPNAS();

	double f(const double *x);
	void df(const double *x, double *g);
	void fdf(const double *x, double *f, double *g);

	void GetMRFvector(const double *x, double *mrfx);

	MRFclass GetMRF(const double *x);

};

#endif /* PROBLEMPNAS_H_ */
