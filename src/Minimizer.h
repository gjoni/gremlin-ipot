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

class Minimizer {

private:

	Minimizer();
	virtual ~Minimizer();

public:

	static MRFclass Minimize(ProblemBase &P);

};

#endif /* MINIMIZER_H_ */
