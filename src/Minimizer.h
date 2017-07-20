/*
 * Minimizer.h
 *
 *  Created on: Jul 17, 2017
 *      Author: ivan
 */

#ifndef MINIMIZER_H_
#define MINIMIZER_H_

#include "ProblemBase.h"

class Minimizer {

private:

	Minimizer();
	virtual ~Minimizer();

public:

	static double Minimize(ProblemBase &P);

};

#endif /* MINIMIZER_H_ */
