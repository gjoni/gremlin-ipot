/*
 * ProblemReducedOMP.h
 *
 *  Created on: Mar 22, 2018
 *      Author: aivan
 */

#ifndef PROBLEMFULLOMP_H_
#define PROBLEMFULLOMP_H_

#include "ProblemBase.h"

class ProblemReducedOMP: public ProblemBase {

private:

	/* regularization parameters */
	double lsingle;
	double lpair;

	/* Vi and Wij dimensions */
	size_t dim1body;
	size_t dim2body;

	/* */
	double *gaux;
	double *ea;
	double *pa;
	double *lpa;

	void Allocate();
	void Free();

public:

	ProblemReducedOMP();
	ProblemReducedOMP(const MSAclass &MSA);
	ProblemReducedOMP(const ProblemReducedOMP &source);
	~ProblemReducedOMP();

	ProblemReducedOMP& operator=(const ProblemReducedOMP &source);

	double f(const double *x);
	void df(const double *x, double *g);
	void fdf(const double *x, double *f, double *g);

	void GetMRFvector(const double *x, double *mrfx);
};

#endif /* PROBLEMFULLOMP_H_ */
