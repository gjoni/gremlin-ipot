/*
 * ProblemRRCE.h
 *
 *  Created on: Jul 26, 2017
 *      Author: ivan
 */

#ifndef PROBLEMRRCE_H_
#define PROBLEMRRCE_H_

#include "RRCE.h"
#include "EigenRRCE.h"
#include "ProblemBase.h"

class ProblemRRCE: public ProblemBase {
private:

	RRCE pot;
	size_t nmodes; /* number of eigenmodes */

	double *em; /* 20 x 20 x nmodes - stores eigenmatrices */

	double lsingle;
	double lpair;

	void Allocate();
	void Free();

public:

	ProblemRRCE();
	ProblemRRCE(const MSAclass &MSA, size_t n);
	~ProblemRRCE();

//	double f(const gsl_vector *x);
//	void df(const gsl_vector *x, gsl_vector *g);
//	void fdf(const gsl_vector *x, double *f, gsl_vector *g);

};

#endif /* PROBLEMRRCE_H_ */
