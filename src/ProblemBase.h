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

#include "MSAclass.h"

/*
 * base abstract class for handling the Gremlin problem
 */

class ProblemBase {

	friend class Minimizer;

protected:

	const MSAclass *MSA;
	unsigned char *msa; /* MSA converted into indices */

	size_t dim; /* problem dimension (number of variables) */

	bool *we; /* for masking (2D array in 1D representation)*/

	/* regularization parameters */
	double lsingle;
	double lpair;

	void SetLsingle(double);
	void SetLpair(double);

	void AllocateBase();
	void FreeBase();

	/* iteration-specific problem parameters */
	std::vector<double> ipar;
	std::vector<double>::iterator ipar_iter;

	virtual void Iterate() = 0;

public:

	ProblemBase();
	ProblemBase(const MSAclass &MSA);
	ProblemBase(const ProblemBase &source);
	virtual ~ProblemBase();

//	ProblemBase& operator=(const ProblemBase &source);

	virtual double f(const double *x) = 0;
	virtual void df(const double *x, double *g) = 0;
	virtual void fdf(const double *x, double *f, double *g) = 0;

	virtual void GetMRFvector(const double *x, double *mrfx) = 0;

	/*
	 * residue indices are as in the REFERENCE sequence
	 */

	/* set e[i,j] to 0, all others to 1 */
	void MaskEdges(const std::vector<std::pair<size_t, size_t> > &e);

	/* set e[i,j] to 1, all others to 0 */
	void UnmaskEdges(const std::vector<std::pair<size_t, size_t> > &e);

	/* default */
	void UnmaskAllEdges();

	void MaskAllEdges();

	size_t GetDim();

	/*
	 * set up iterations schedule
	 */
	void SetUpIterations(const std::vector<double>&);

};

#endif /* PROBLEMBASE_H_ */
