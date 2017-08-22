/*
 * ProblemRRCE.h
 *
 *  Created on: Jul 26, 2017
 *      Author: ivan
 */

#ifndef PROBLEMRRCE_H_
#define PROBLEMRRCE_H_

#include "RRCE.h"
//#include "EigenRRCE.h"
#include "ProblemBase.h"

class ProblemRRCE: public ProblemBase {
private:

	RRCE pot;
	size_t nmodes; /* number of eigenmodes */

	double *em; /* 20 x 20 x nmodes - stores eigenmatrices */
	double *ev; /* eigenvalues */

	double lsingle;
	double lpair;

	size_t dim1body;
	size_t dim2body;

	void Allocate();
	void Free();

public:

	ProblemRRCE();
	ProblemRRCE(const MSAclass &MSA, size_t n);
	~ProblemRRCE();

	double f(const double *x);
	void df(const double *x, double *g);
	void fdf(const double *x, double *f, double *g);

	void GetMRFvector(const double *x, double *mrfx);

	MRFclass GetMRF(const double *x);

};

#endif /* PROBLEMRRCE_H_ */
