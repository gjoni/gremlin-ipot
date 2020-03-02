/*
 * ProblemSmoothL1_1b.h
 *
 *  Created on: Aug 30, 2018
 *      Author: aivan
 */

#ifndef PROBLEMSMOOTHL1_1B_H_
#define PROBLEMSMOOTHL1_1B_H_

#include "ProblemL2_1b.h"

class ProblemSmoothL1_1b: public ProblemL2_1b {

private:

	double eps; /* L1 smoothing parameter */

	double Reg_f(const double *x);
	double Reg_fdf(const double *x, double *g);

	void Iterate() {
//		eps /= 2;
	}

public:

	ProblemSmoothL1_1b();
	ProblemSmoothL1_1b(const MSAclass &MSA);
	ProblemSmoothL1_1b(const ProblemL2_1b &source);
	~ProblemSmoothL1_1b();

	ProblemSmoothL1_1b& operator=(const ProblemL2_1b &source);

	void SetSmoothingFactor(double eps);

};

#endif /* PROBLEMSMOOTHL1_1B_H_ */
