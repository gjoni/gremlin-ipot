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

MRFclass Minimizer::MinimizeLBFGS(ProblemBase &P, int niter) {

	assert(sizeof(double) == sizeof(lbfgsfloatval_t)); /* sse2 disabled */

	size_t dim = P.GetDim();

	if (!dim) {
		printf("Error: cannot run Gremlin for an empty MSA\n");
		exit(1);
	}

	lbfgsfloatval_t fx;
	lbfgsfloatval_t *m_x = lbfgs_malloc(dim);

	if (m_x == NULL) {
		printf("ERROR: Failed to allocate a memory block for variables.\n");
		exit(1);
	}

	/* init the variables */
	memset(m_x, 0, dim * sizeof(lbfgsfloatval_t));

	printf("# %-8s%-14s%-14s%-14s%-8s%-12s\n", "iter", "f(x)", "||x||", "||g||",
			"neval", "epsilon");

	printf("# %-8d%-12.5e\n", 0, P.f(m_x));

	lbfgs_parameter_t param;
	lbfgs_parameter_init(&param);
	param.max_iterations = niter;
	param.m = 3;
//	param.max_linesearch = 3;
//	param.gtol = 0.1;

	int ret = lbfgs(dim, m_x, &fx, _evaluate, _progress, &P, &param);
	printf("# L-BFGS termination code = %d", ret);
	if (ret == -997) {
		printf(" (LBFGSERR_MAXIMUMITERATION)");
	}
	printf("\n");

	/* TODO: a lot of RAM overhead  -- simplify */
	size_t dim_mrf = P.MSA->GetNcol() * MSAclass::NAA
			* (1 + P.MSA->GetNcol() * MSAclass::NAA);
	double *mrfx = (double*) malloc(dim_mrf * sizeof(double));
	P.GetMRFvector(m_x, mrfx);
	MRFclass MRF(mrfx, P.we, P.MSA);

	lbfgs_free(m_x);
	free(mrfx);

	return MRF;

}

void Minimizer::MinimizeLBFGS(ProblemBase &P, int niter, MRFclassNew &MRF) {

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
	param.m = 3;

	int ret = lbfgs(dim, x, &fx, _evaluate, _progress, &P, &param);
	printf("# L-BFGS termination code = %d\n", ret);

}

lbfgsfloatval_t Minimizer::_evaluate(void *instance, const lbfgsfloatval_t *x,
		lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {

	ProblemBase *_P = (ProblemBase*) instance;
	lbfgsfloatval_t f;
	_P->fdf(x, &f, g);
	_P->Iterate();
	return f;

}

int Minimizer::_progress(void *instance, const lbfgsfloatval_t *x,
		const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
		const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
		const lbfgsfloatval_t step, int n, int k, int ls) {

	printf("# %-8d%-12.5e  %-12.5e  %-12.5e  %-6d  %-10.5f\n", k, fx, xnorm,
			gnorm, ls, gnorm / xnorm);
	fflush(stdout);

	return 0;
}

