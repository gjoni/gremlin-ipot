/*
 * ProblemL2.h
 *
 *  Created on: Mar 22, 2018
 *      Author: aivan
 */

#ifndef PROBLEML2_H_
#define PROBLEML2_H_

#include "ProblemBase.h"

class ProblemL2: public ProblemBase {

protected:

	/* Vi and Wij dimensions */
	size_t dim1body;
	size_t dim2body;
	size_t dim2reduced;

	/* arrays for storing temp vars */
	double *gaux;
	double *ea;
	double *pa;
	double *lpa;

	void Allocate();
	void Free();

	size_t To1D(size_t i, size_t j);
	void To2D(size_t k, size_t &i, size_t &j);

	virtual double Reg_f(const double *x);
	virtual double Reg_fdf(const double *x, double *g);

	void Iterate();

public:

	ProblemL2();
	ProblemL2(const MSAclass &MSA);
	ProblemL2(const ProblemL2 &source);
	~ProblemL2();

	ProblemL2& operator=(const ProblemL2 &source);

	double f(const double *x);
	void df(const double *x, double *g);
	void fdf(const double *x, double *f, double *g);

	void GetMRFvector(const double *x, double *mrfx);

};

#endif /* PROBLEML2_H_ */
