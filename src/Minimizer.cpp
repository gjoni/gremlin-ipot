/*
 * Minimizer.cpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ivan
 */

#include <cstring>
#include <gsl/gsl_multimin.h>

#include "Minimizer.h"
#include "ProblemStatic.h"

Minimizer::Minimizer() {
	// TODO Auto-generated constructor stub

}

Minimizer::~Minimizer() {
	// TODO Auto-generated destructor stub
}

double Minimizer::Minimize(ProblemBase &P) {

	double f = 0.0;

	gsl_multimin_function_fdf gsl_func;

	size_t dim = P.GetDim();
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

	gsl_multimin_fdfminimizer_set(s, &gsl_func, x, 1.0, 0.2);

	/*
	 * minimization protocol
	 */

	/* clean memory */
	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(x);

	return f;

}
