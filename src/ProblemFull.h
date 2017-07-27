/*
 * ProblemFull.h
 *
 *  Created on: Jul 17, 2017
 *      Author: ivan
 */

// Standard GREMLIN2 protocol with L2 regularization
// The code for computing the objective function and gradients
// is largely adopted from CCMpred by J.Soeding group
// https://github.com/soedinglab/CCMpred
// V20170720 - first functional class
#ifndef PROBLEMFULL_H_
#define PROBLEMFULL_H_

#include "ProblemBase.h"
#include "MSAclass.h"

class ProblemFull: public ProblemBase {

private:

	/* regularization parameters */
	double lsingle;
	double lpair;

	/* aux array to store asymmetric gradient */
	size_t dim2body;
	double *gaux;

	void Allocate();
	void Free();

public:

	ProblemFull();
	ProblemFull(const MSAclass &MSA);
	ProblemFull(const ProblemFull &source);
	~ProblemFull();

	/* TODO: set reg. parameters */
//	void SetLH();
//	void SetLJ();

	ProblemFull& operator=(const ProblemFull &source);

	/* TODO: include edge masks we[] in the calculations */
	/* TODO: change gsl_vector to double* */
	double f(const gsl_vector *x);
	void df(const gsl_vector *x, gsl_vector *g);
	void fdf(const gsl_vector *x, double *f, gsl_vector *g);

};

#endif /* PROBLEMFULL_H_ */
