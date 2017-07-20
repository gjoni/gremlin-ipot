/*
 * ProblemBase.h
 *
 *  Created on: Jul 17, 2017
 *      Author: ivan
 */

#ifndef PROBLEMBASE_H_
#define PROBLEMBASE_H_

#include <vector>
#include <utility>

#include <gsl/gsl_vector.h>

#include "MSAclass.h"

/*
 * base abstract class for handling the Gremlin problem
 */

class ProblemBase {
protected:

	MSAclass *MSA;

	size_t dim; /* problem dimension (number of variables) */

	double *w; /* sequence weights */
	double **we; /* for masking or biasing */

	void AllocateBase();
	void FreeBase();

public:

	ProblemBase();
	ProblemBase(MSAclass &MSA);
	virtual ~ProblemBase();

//	ProblemBase& operator=(const ProblemBase &source);

	virtual double f(const gsl_vector *x) = 0;
	virtual void df(const gsl_vector *x, gsl_vector *g) = 0;
	virtual void fdf(const gsl_vector *x, double *f, gsl_vector *g) = 0;

	void Reweight(double t = 0.8);

	/*
	 * residue indices are as in the REFERENCE sequence
	 */

	/* force e[i,j] to be equal to 0, all others to 1 */
	void MaskEdges(const std::vector<std::pair<int, int> > &e);

	/* force e[i,j] to be equal to 1, all others to 0 */
	void UnmaskEdges(const std::vector<std::pair<int, int> > &e);
	void UnmaskEdges();

	size_t GetDim();

};

#endif /* PROBLEMBASE_H_ */
