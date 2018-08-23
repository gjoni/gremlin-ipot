/*
 * Problem1Body.h
 *
 *  Created on: Aug 23, 2018
 *      Author: aivan
 */

#ifndef PROBLEM1BODY_H_
#define PROBLEM1BODY_H_

#include "ProblemBase.h"

class Problem1Body: public ProblemBase {

private:

	/* 1-body regularization parameter */
	double lsingle;

	/* Vi dimensions */
	size_t dim1body;

	/* arrays for storing temp vars */
	double *ea;
	double *pa;
	double *lpa;

	void Allocate();
	void Free();

public:

	Problem1Body();
	Problem1Body(const MSAclass &MSA);
	Problem1Body(const Problem1Body &source);
	~Problem1Body();

	Problem1Body& operator=(const Problem1Body &source);

	double f(const double *x);
	void df(const double *x, double *g);
	void fdf(const double *x, double *f, double *g);

	void GetMRFvector(const double *x, double *mrfx);

};

#endif /* PROBLEM1BODY_H_ */
