/*
 * EigenRRCE.h
 *
 *  Created on: Jun 1, 2017
 *      Author: ivan
 */

#ifndef EIGENRRCE_H_
#define EIGENRRCE_H_

/*
 * a class for spectral decomposition of contact matrices
 */

class EigenRRCE {

private:

	const int dim = 20; /* operates on 20x20 matrices */

	double **J;

	double *e; /* eigenvalues */
	double **ev; /* eigenvectors - BY ROWS */

	void Decompose();

	void Allocate();
	void Free();

public:

	EigenRRCE(double **J);
	EigenRRCE();

	~EigenRRCE();

	/*
	 * !!! eigenvalues and eigenvectors are indexed 1..20 !!!
	 */

	/* energy for residues a,b from the i-th eigenmode:
	 * E(a,b) = e[i]*ev[i][a]*ev[i][b]
	 * negative i: do the sum over over (i+1)..20 eigenmodes */
	double GetEigenEnergy(const int a, const int b, const int i);

	/* reconstruction (relative) error for first i eigenmodes:
	 * |r1-r2|/|r1+r2|
	 * negative i: do the sum over over (i+1)..20 eigenmodes */
	double GetReconstructionError(const int i);

	/* correlation between original and reconstructed using
	 * first i eigenmodes contact matrices */
	double GetReconstructionCorrel(const int i);

	double GetEigenvalue(const int i);

	/* i-th eigenmode for the interaction matrix: m_i[] = e_i[] ⊗ e_i[]*/
	void GetEigenmatrix(int i, double **em);

};

#endif /* EIGENRRCE_H_ */
