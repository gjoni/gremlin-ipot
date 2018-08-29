/*
 * ProblemL2_1b.h
 *
 *  Created on: Aug 23, 2018
 *      Author: aivan
 */

#ifndef PROBLEML2_1B_H_
#define PROBLEML2_1B_H_

#include "ProblemBase.h"

class ProblemL2_1b: public ProblemBase {

private:

	/* Vi dimensions */
	size_t dim1body;

	/* arrays for storing temp vars */
	double *ea;
	double *pa;
	double *lpa;

	void Allocate();
	void Free();

	/* TODO: implement */
	double Reg_f(const double *x);
	double Reg_fdf(const double *x, double *g);

	void Iterate();

public:

	ProblemL2_1b();
	ProblemL2_1b(const MSAclass &MSA);
	ProblemL2_1b(const ProblemL2_1b &source);
	~ProblemL2_1b();

	ProblemL2_1b& operator=(const ProblemL2_1b &source);

	double f(const double *x);
	void df(const double *x, double *g);
	void fdf(const double *x, double *f, double *g);

	void GetMRFvector(const double *x, double *mrfx);

};

#endif /* PROBLEML2_1B_H_ */
