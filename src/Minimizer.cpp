/*
 * Minimizer.cpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ivan
 */

#include <cstring>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>

#include "Minimizer.h"
#include "ProblemStatic.h"
#include "MRFclass.h"

Minimizer::Minimizer() {
	// TODO Auto-generated constructor stub

}

Minimizer::~Minimizer() {
	// TODO Auto-generated destructor stub
}

MRFclass Minimizer::Minimize(ProblemBase &P) {

	size_t dim = P.GetDim();

	if (!dim) {
		printf("Error: cannot run Gremlin for an empty MSA\n");
		exit(1);
	}

//	double f = 0.0;

	gsl_multimin_function_fdf gsl_func;

	gsl_func.n = dim;
	gsl_func.f = ProblemStatic::f;
	gsl_func.df = ProblemStatic::df;
	gsl_func.fdf = ProblemStatic::fdf;
	gsl_func.params = &P;

	const gsl_multimin_fdfminimizer_type *T;
	T = gsl_multimin_fdfminimizer_vector_bfgs2;

	gsl_multimin_fdfminimizer *s;
	s = gsl_multimin_fdfminimizer_alloc(T, dim);

	gsl_vector *x = gsl_vector_alloc(dim);
	memset(x->data, 0, dim * sizeof(double));

	gsl_multimin_fdfminimizer_set(s, &gsl_func, x, 10.0, 0.2);

	printf("# %-8s%-14s%-12s%-12s%-12s%-12s\n", "iter", "f(x)", "||x||",
			"||dx||", "step", "reltol");

	printf("# %-8d%-10.2f\n", 0, gsl_func.f /* f_full_gpl */(x, &P));

	int status, iter = 0;
	do {
		iter++;

		status = gsl_multimin_fdfminimizer_iterate(s);
		if (status) {
			break;
		}
		status = gsl_multimin_test_gradient(s->gradient, 1e-1);

		double xnorm = gsl_blas_dnrm2(s->x);
		double gnorm = gsl_blas_dnrm2(s->gradient);
		double snorm = gsl_blas_dnrm2(s->dx); /* step */

		printf("# %-8d%-10.2f    %-10.2f  %-10.2f  %-10.2f  %-10.5f\n", iter,
				s->f, xnorm, gnorm, snorm, snorm / xnorm);

		if (status == GSL_SUCCESS) {
			printf("# Minimum found at iteration %d!!!\n", iter);
		}

	} while (status == GSL_CONTINUE && iter < 10);

	/*
	 * minimization protocol
	 */

	/* save results */
	MRFclass MRF(s->x->data, s->x->data + P.MSA->GetNcol() * MSAclass::NAA,
			P.MSA->GetNcol());

	/* clean memory */
	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(x);

	return MRF;

}
