/*
 * ProblemReducedOMP.h
 *
 *  Created on: Mar 22, 2018
 *      Author: aivan
 */

#ifndef PROBLEMREDUCEDOMP_H_
#define PROBLEMREDUCEDOMP_H_

#include "ProblemBase.h"

class ProblemReducedOMP: public ProblemBase {

private:

	/* regularization parameters */
	double lsingle;
	double lpair;

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

	/* TODO: regularization functions
	 * (to allow for easy inheritance) */
	double Reg_f(const double *x);
	double Reg_fdf(const double *x, double *g);

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

#endif /* PROBLEMREDUCEDOMP_H_ */
