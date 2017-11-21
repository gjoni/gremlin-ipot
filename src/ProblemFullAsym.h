/*
 * ProblemFullAsym.h
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

#ifndef PROBLEMFULLASYM_H_
#define PROBLEMFULLASYM_H_

#include "ProblemBase.h"
#include "MSAclass.h"
#include "RRCE.h"

class ProblemFullAsym: public ProblemBase {

private:

	/* regularization parameters */
	double lsingle;
	double lpair;
	double lskew;

	/* Vi and Wij dimensions */
	size_t dim1body;
	size_t dim2body;

	/* asymmetric priors */
//	RRCE pot;
	double **J;
	double **dx;
	double **dy;
	void SetPriors();
	double fPrior(double x, double mu);
	double dfPrior(double x, double mu);
	double ddfPrior(double x, double mu);

	double Prior(double x, size_t a, size_t b);
	double dPrior(double x, size_t a, size_t b);


	void Allocate();
	void Free();

public:

	ProblemFullAsym();
	ProblemFullAsym(const MSAclass &MSA);
	ProblemFullAsym(const ProblemFullAsym &source);
	~ProblemFullAsym();

	/* TODO: set reg. parameters */
//	void SetLH();
//	void SetLJ();
	ProblemFullAsym& operator=(const ProblemFullAsym &source);

	double f(const double *x);
	void df(const double *x, double *g);
	void fdf(const double *x, double *f, double *g);

	void GetMRFvector(const double *x, double *mrfx);

	void SetLpair(double l);
	void SetLskew(double l);
	void SetL2(double lpair, double lskew);

};

#endif /* PROBLEMFULLASYM_H_ */
