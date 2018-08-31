/*
 * ProblemSmoothL1.h
 *
 *  Created on: Aug 30, 2018
 *      Author: aivan
 */

#ifndef PROBLEMSMOOTHL1_H_
#define PROBLEMSMOOTHL1_H_

#include "ProblemL2.h"

class ProblemSmoothL1: public ProblemL2 {

private:

	double eps; /* L1 smoothing parameter */

	double Reg_f(const double *x);
	double Reg_fdf(const double *x, double *g);

	void Iterate() {
//		eps /= 2;
	}

public:

	ProblemSmoothL1();
	ProblemSmoothL1(const MSAclass &MSA);
	ProblemSmoothL1(const ProblemL2 &source);
	~ProblemSmoothL1();

	ProblemSmoothL1& operator=(const ProblemL2 &source);

};

#endif /* PROBLEMSMOOTHL1_H_ */
