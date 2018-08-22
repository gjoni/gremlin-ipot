/*
 * ProblemFull.h
 *
 *  Created on: Jul 17, 2017
 *      Author: ivan
 */

/*
 * Standard GREMLIN2 protocol with L2 regularization
 * The code for computing the objective function and gradients
 * is largely adopted from CCMpred by J.Soeding group
 * https://github.com/soedinglab/CCMpred
 * -----
 * V20170804 - alpha version
 * V20170720 - first functional class
 */

#ifndef PROBLEMFULL_H_
#define PROBLEMFULL_H_

#include "ProblemBase.h"
#include "MSAclass.h"

class ProblemFull: public ProblemBase {

private:

	/* regularization parameters */
	double lsingle;
	double lpair;

	/* Vi and Wij dimensions */
	size_t dim1body;
	size_t dim2body;

public:

	ProblemFull();
	ProblemFull(const MSAclass &MSA);
	ProblemFull(const ProblemFull &source);
	~ProblemFull();

	/* TODO: set reg. parameters */
//	void SetLH();
//	void SetLJ();
	ProblemFull& operator=(const ProblemFull &source);

	double f(const double *x);
	void df(const double *x, double *g);
	void fdf(const double *x, double *f, double *g);

	void GetMRFvector(const double *x, double *mrfx);

//	size_t GetDim();

};

#endif /* PROBLEMFULL_H_ */
