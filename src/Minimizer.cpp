/*
 * Minimizer.cpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ivan
 */

#include <cstring>
#include <cassert>

#include "lbfgs.h"

#include "Minimizer.h"

Minimizer::Minimizer() {
	// TODO Auto-generated constructor stub

}

Minimizer::~Minimizer() {
	// TODO Auto-generated destructor stub
}

void Minimizer::MinimizeLBFGS(ProblemBase &P, int niter, MRFclass &MRF) {

	assert(sizeof(double) == sizeof(lbfgsfloatval_t));

	size_t dim = P.GetDim();

	if (!dim) {
		printf("Error: cannot run Gremlin for an empty MSA\n");
		exit(1);
	}

	lbfgsfloatval_t fx;
	lbfgsfloatval_t *x = MRF.GetX();

	printf("# %-8s%-14s%-14s%-14s%-8s%-12s\n", "iter", "f(x)", "||x||", "||g||",
			"neval", "epsilon");

	lbfgs_parameter_t param;
	lbfgs_parameter_init(&param);

	/* custom params */
	param.max_iterations = niter;
//	param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
//	param.min_step = 1e-40;
	param.m = 3;
//	param.orthantwise_c = 1.;

	int ret = lbfgs(dim, x, &fx, _evaluate, _progress, &P, &param);
	printf("# L-BFGS termination code = %d\n", ret);

}

lbfgsfloatval_t Minimizer::_evaluate(void *instance, const lbfgsfloatval_t *x,
		lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {

	ProblemBase *_P = (ProblemBase*) instance;
	lbfgsfloatval_t f;
	_P->fdf(x, &f, g);
	return f;

}

int Minimizer::_progress(void *instance, const lbfgsfloatval_t *x,
		const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
		const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
		const lbfgsfloatval_t step, int n, int k, int ls) {

	printf("# %-8d%-12.5e  %-12.5e  %-12.5e  %-6d  %-10.5f\n", k, fx, xnorm,
			gnorm, ls, gnorm / xnorm);
	fflush(stdout);

	ProblemBase *_P = (ProblemBase*) instance;
	if (k % 3 == 0) {
		_P->Iterate();
	}

	return 0;
}

