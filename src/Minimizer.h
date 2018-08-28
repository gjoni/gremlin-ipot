/*
 * Minimizer.h
 *
 *  Created on: Jul 17, 2017
 *      Author: ivan
 */

#ifndef MINIMIZER_H_
#define MINIMIZER_H_

#include "ProblemBase.h"
#include "MRFclass.h"

#include "lbfgs.h"

class Minimizer {

private:

	Minimizer();
	virtual ~Minimizer();

	static lbfgsfloatval_t _evaluate(void *instance, const lbfgsfloatval_t *x,
			lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

	static int _progress(void *instance, const lbfgsfloatval_t *x,
			const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
			const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
			const lbfgsfloatval_t step, int n, int k, int ls);

public:

	static void MinimizeLBFGS(ProblemBase &P, int niter, MRFclass &MRF);

	/* TODO: CG minimizer - more memory efficient */
	// static MRFclass MinimizeCG(ProblemBase &P, int niter);
	/* TODO: ADMM minimizer for L1 penalty */

};

#endif /* MINIMIZER_H_ */
