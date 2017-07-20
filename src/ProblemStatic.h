/*
 * ProblemStatic.h
 *
 *  Created on: Jul 17, 2017
 *      Author: ivan
 */

#ifndef PROBLEMSTATIC_H_
#define PROBLEMSTATIC_H_

#include <gsl/gsl_vector.h>

/*
 * wrapper class to enable GSL usage
 */

class ProblemStatic {

private:

	ProblemStatic();
	~ProblemStatic();

public:

	static double f(const gsl_vector *x, void *par);
	static void df(const gsl_vector *x, void *par, gsl_vector *g);
	static void fdf(const gsl_vector *x, void *par, double *f, gsl_vector *g);

};

#endif /* PROBLEMSTATIC_H_ */
