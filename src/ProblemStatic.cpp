/*
 * ProblemStatic.cpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ivan
 */

#include "ProblemStatic.h"
#include "ProblemBase.h"

ProblemStatic::ProblemStatic() {
	// TODO Auto-generated constructor stub
}

ProblemStatic::~ProblemStatic() {
	// TODO Auto-generated destructor stub
}

double ProblemStatic::f(const gsl_vector *x, void *par) {

	ProblemBase *P = (ProblemBase*) par;
	return P->f(x);

}

void ProblemStatic::df(const gsl_vector *x, void *par, gsl_vector *g) {

	ProblemBase *P = (ProblemBase*) par;
	P->df(x, g);

}

void ProblemStatic::fdf(const gsl_vector *x, void *par, double *f,
		gsl_vector *g) {

	ProblemBase *P = (ProblemBase*) par;
	P->fdf(x, f, g);

}
