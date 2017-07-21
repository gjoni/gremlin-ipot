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
	double *gaux;

	void Allocate();
	void Free();

public:

	ProblemFull(MSAclass &MSA);
	ProblemFull();
	~ProblemFull();

	double f(const gsl_vector *x);
	void df(const gsl_vector *x, gsl_vector *g);
	void fdf(const gsl_vector *x, double *f, gsl_vector *g);

};

#endif /* PROBLEMFULL_H_ */
