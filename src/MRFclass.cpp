/*
 * MRFclass.cpp
 *
 *  Created on: Jul 20, 2017
 *      Author: ivan
 */

#include <cstring>
#include <cassert>

#include "MRFclass.h"

MRFclass::MRFclass() :
		dim(0), dimh(0), dimJ(0), len_ref(0), mtx_to_a3m(0), a3m_to_mtx(0), h(
		NULL), J(NULL), we(
		NULL) {

	/* */

}

MRFclass::MRFclass(const MRFclass &s) :
		dim(s.dim), dimh(s.dimh), dimJ(s.dimJ), len_ref(s.len_ref), mtx_to_a3m(
				0), a3m_to_mtx(0), h(NULL), J(
		NULL), we(
		NULL) {

	Allocate();

	memcpy(h, s.h, dimh * sizeof(double));
	memcpy(J, s.J, dimJ * sizeof(double));
	memcpy(we, s.we, dim * dim * sizeof(bool));

}

MRFclass::MRFclass(double *x, const MSAclass *MSA) :
		dim(MSA->GetNcol()), len_ref(MSA->GetLen()), mtx_to_a3m(MSA->col_map), a3m_to_mtx(
				MSA->a3m_to_msa), h(NULL), J(NULL), we(NULL) {

	size_t NAA = MSAclass::NAA;

	dimh = dim * NAA;
	dimJ = dimh * dimh;

	double *h_ = x;
	double *J_ = x + dimh;

	Allocate();

	memcpy(h, h_, dimh * sizeof(double));
	memcpy(J, J_, dimJ * sizeof(double));

//	for (size_t i = 0; i < dim; i++) {
//		for (size_t a = 0; a < NAA; a++) {
//			h[i * NAA + a] = h_[a * dim + i];
//		}
//	}
//
//	for (size_t i = 0; i < dim; i++) {
//		for (size_t a = 0; a < NAA; a++) {
//			for (size_t j = 0; j < dim; j++) {
//				for (size_t b = 0; b < NAA; b++) {
//					J[(i * dim + j) * NAA * NAA + a * NAA + b] = J_[((a * dim
//							+ i) * NAA + b) * dim + j];
//				}
//			}
//		}
//	}

	for (size_t i = 0; i < dim * dim; i++) {
		we[i] = true;
	}

}

MRFclass::MRFclass(double *x, bool *we_, const MSAclass *MSA) :
		MRFclass(x, MSA) {

	memcpy(we, we_, dim * dim * sizeof(bool));

}

MRFclass::MRFclass(const std::string &name) :
		dim(0), dimh(0), dimJ(0), len_ref(0), mtx_to_a3m(0), a3m_to_mtx(0), h(
		NULL), J(NULL), we(NULL) {

	FILE *F = fopen(name.c_str(), "r");
	if (F == NULL) {
		printf("Error: cannot open '%s' file to save MRF\n", name.c_str());
		exit(1);
	}

	size_t NAA = MSAclass::NAA;

	/* dimensions (seq length) */
	fscanf(F, "%lu %lu\n", &dim, &len_ref);
	dimh = dim * NAA;
	dimJ = dimh * dimh;

	Allocate();

	/* read maps */
	mtx_to_a3m.resize(dim);
	a3m_to_mtx.resize(len_ref);
	for (size_t i = 0; i < dim; i++) {
		fscanf(F, "%lu ", &(mtx_to_a3m[i]));
	}
	for (size_t i = 0; i < len_ref; i++) {
		fscanf(F, "%lu ", &(a3m_to_mtx[i]));
	}

	/* edge weights */
	{
		size_t SIZE = dim + 2;
		char *buf = (char*) malloc(SIZE * sizeof(char));
		for (size_t i = 0; i < dim; i++) {
			fgets(buf, SIZE, F);
			for (size_t j = 0; j < dim; j++) {
				char c = buf[j];
				if (c == '1') {
					we[i * dim + j] = true;
				} else {
					we[i * dim + j] = false;
				}
			}
		}
		free(buf);
	}

	/* local fields */
	for (size_t i = 0; i < dim; i++) {
		for (size_t a = 0; a < NAA; a++) {
			fscanf(F, "%lf ", h + (i * NAA + a));
		}
	}

	/* couplings */
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			double *Jp = J + (i * dim + j) * NAA * NAA;
			for (size_t aa = 0; aa < NAA * NAA; aa++) {
				fscanf(F, "%lf ", Jp++);
			}
		}
	}

	fclose(F);

}

MRFclass& MRFclass::operator=(const MRFclass &source) {

	assert(this != &source); /* an attempt to assign Residue to itself */

	Free();

	dim = source.dim;
	dimh = source.dimh;
	dimJ = source.dimJ;

	len_ref = source.len_ref;
	mtx_to_a3m = source.mtx_to_a3m;
	a3m_to_mtx = source.a3m_to_mtx;

	Allocate();

	memcpy(h, source.h, dimh * sizeof(double));
	memcpy(J, source.J, dimJ * sizeof(double));
	memcpy(we, source.we, dim * dim * sizeof(bool));

	return *this;

}

MRFclass::~MRFclass() {

	Free();

}

void MRFclass::Allocate() {

	h = (double*) malloc(dimh * sizeof(double));
	J = (double*) malloc(dimJ * sizeof(double));
	we = (bool*) malloc(dim * dim * sizeof(bool));

}

void MRFclass::Free() {

	if (h != NULL) {
		free(h);
	}
	if (J != NULL) {
		free(J);
	}
	if (we != NULL) {
		free(we);
	}

}

size_t MRFclass::GetDim() const {

	return dim;

}

void MRFclass::Save(const std::string &name) const {

	FILE *F = fopen(name.c_str(), "w");
	if (F == NULL) {
		printf("Error: cannot open '%s' file to save MRF\n", name.c_str());
		exit(1);
	}

	size_t NAA = MSAclass::NAA;

	/* dimension (seq length) */
	fprintf(F, "%ld %ld\n", dim, len_ref);

	/* save maps */
	for (size_t i = 0; i < dim; i++) {
		fprintf(F, "%lu ", mtx_to_a3m[i]);
	}
	fprintf(F, "\n");
	for (size_t i = 0; i < len_ref; i++) {
		fprintf(F, "%lu ", a3m_to_mtx[i]);
	}
	fprintf(F, "\n");

	/* edge weights */
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			if (we[i * dim + j] == 1) {
				fputc('1', F);
			} else {
				fputc('0', F);
			}
		}
		fprintf(F, "\n");
	}

	/* local fields */
	for (size_t i = 0; i < dim; i++) {
		for (size_t a = 0; a < NAA; a++) {
			fprintf(F, "%.5e ", h[i * NAA + a]);
		}
		fprintf(F, "\n");
	}

	/* couplings */
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			double *Jp = J + (i * dim + j) * NAA * NAA;
			for (size_t aa = 0; aa < NAA * NAA; aa++) {
				fprintf(F, "%.5e ", *Jp++);
			}
			fprintf(F, "\n");
		}
	}

	fclose(F);

}

double MRFclass::GetPairEnergies(const MSAclass &MSA,
		const std::vector<std::pair<size_t, size_t> > &contacts) const {

	assert(MSA.GetNcol() == dim); /* MRF and cleaned MSA length mismatch */

	double E = 0.0;

	size_t nrow = MSA.GetNrow();
	size_t ncol = MSA.GetNcol();
	size_t NAA = MSAclass::NAA;

	/* prepare cleaned MSA */
	unsigned char *msa = MSA.GetMsa();
	MSAclass::aatoi(msa, nrow * ncol);

	/* loop over all sequences in the msa[] */
	for (size_t i = 0; i < nrow; i++) {

		/* current sequence */
		unsigned char *seq = msa + i * ncol;

		/* loop over all contacts */
		for (const auto& c : contacts) {
			size_t a = seq[c.first];
			size_t b = seq[c.second];

			/* omit pairs with gaps */
			if (a < NAA - 1 && b < NAA - 1) {
				E += J[(c.first * ncol + c.second) * NAA * NAA + a * NAA + b];
			}

		}
	}

	return E;

}

double MRFclass::GetPairEnergies(const unsigned char *msa, size_t nrow,
		const std::vector<std::pair<size_t, size_t> > &contacts) const {

	double E = 0.0;

	size_t NAA = MSAclass::NAA;

	/* loop over all sequences in the msa[] */
	for (size_t i = 0; i < nrow; i++) {

		/* current sequence */
		const unsigned char *seq = msa + i * dim;

		/* loop over all contacts */
		for (const auto& c : contacts) {
			size_t a = seq[c.first];
			size_t b = seq[c.second];

			/* omit pairs with gaps */
			if (a < NAA - 1 && b < NAA - 1) {
				E += J[(c.first * dim + c.second) * NAA * NAA + a * NAA + b];
			}

		}
	}

	return E;

}

std::vector<double> MRFclass::ScoreMSA(const MSAclass &MSA) {

	std::vector<double> v;

	size_t NAA = MSAclass::NAA;

	unsigned char *seq = (unsigned char*) malloc(
			MSA.ncol * sizeof(unsigned char));

	for (size_t i = 0; i < MSA.nrow; i++) {

		/* convert A3M sequence into indices */
		const std::string &a3m_seq = MSA.a3m[MSA.row_map[i]].second;
		for (size_t j = 0; j < MSA.ncol; j++) {
			seq[j] = MSAclass::aatoi(a3m_seq[MSA.col_map[j]]);
		}

		double E = 0.0;

		/* local fields */
		size_t nnongap = 0;
		for (size_t j = 0; j < MSA.ncol; j++) {
			E += h[j * NAA + seq[j]];
			if (seq[j] < NAA - 1) {
				nnongap++;
			}
		}

		/* couplings */
		for (size_t j = 0; j < MSA.ncol; j++) {
			unsigned char a = seq[j];
			for (size_t k = j + 1; k < MSA.ncol; k++) {
				unsigned char b = seq[k];
//				if (a < NAA - 1 && b < NAA - 1) {
				E += J[j * MSA.ncol + k] * NAA * NAA + (a * NAA + b);
//				}
			}
		}

		v.push_back(-E);

	}

	free(seq);

	return v;

}
