/*
 * ProblemBase.cpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ivan
 */

#include <cstring>
#include <cmath>

#include "ProblemBase.h"

ProblemBase::ProblemBase() :
		MSA(NULL), msa(NULL), dim(0), we(NULL), lsingle(0), lpair(0), ipar(), ipar_iter(
				ipar.begin()) {

	/* */

}

ProblemBase::ProblemBase(const MSAclass &MSA_) :
		MSA(&MSA_), msa(NULL), dim(MSA_.GetNcol() * MSA_.GetNrow()), we(
		NULL), lsingle(0), lpair(0), ipar(), ipar_iter(ipar.begin()) {

	AllocateBase();

	msa = MSA->GetMsa();

	MSAclass::aatoi(msa, dim);

	UnmaskAllEdges();

}

ProblemBase::ProblemBase(const ProblemBase &source) :
		MSA(source.MSA), msa(NULL), dim(source.MSA->ncol * source.MSA->nrow), we(
		NULL), lsingle(source.lsingle), lpair(source.lpair), ipar(source.ipar), ipar_iter(
				ipar.begin()) {

	AllocateBase();

	size_t msa_dim = MSA->nrow * MSA->ncol;

	msa = (unsigned char *) malloc(msa_dim * sizeof(unsigned char));

	ipar_iter += source.ipar_iter - source.ipar.begin();

	memcpy(msa, source.msa, msa_dim * sizeof(unsigned char));
	memcpy(we, source.we, MSA->ncol * MSA->ncol * sizeof(bool));

}

ProblemBase::~ProblemBase() {

	FreeBase();

}

void ProblemBase::MaskEdges(const std::vector<std::pair<size_t, size_t> > &e) {

//	/* set all we[][] to 1.0 */
//	UnmaskAllEdges();

	/* mask the specified edges */
	for (auto const& edge : e) {
		size_t i = MSA->GetMsaIdx(edge.first);
		size_t j = MSA->GetMsaIdx(edge.second);
		if (i < SIZE_MAX && j < SIZE_MAX) {
			we[i * MSA->ncol + j] = false;
			we[j * MSA->ncol + i] = false;
		}

	}

}

void ProblemBase::UnmaskEdges(
		const std::vector<std::pair<size_t, size_t> > &e) {

//	/* set all we[][] to false */
//	memset(we, 0, MSA->ncol * MSA->ncol * sizeof(bool));

	/* mask the specified edges */
	for (auto const& edge : e) {
		size_t i = MSA->GetMsaIdx(edge.first);
		size_t j = MSA->GetMsaIdx(edge.second);
		if (i < SIZE_MAX && j < SIZE_MAX) {
			we[i * MSA->ncol + j] = true;
			we[j * MSA->ncol + i] = true;
		}
	}

}

void ProblemBase::UnmaskAllEdges() {

	/* set all we[][] to true */
	for (size_t i = 0; i < MSA->ncol * MSA->ncol; i++) {
		we[i] = true;
	}

}

void ProblemBase::MaskAllEdges() {

	/* set all we[][] to false */
	memset(we, 0, MSA->ncol * MSA->ncol * sizeof(bool));

}

void ProblemBase::AllocateBase() {

	we = (bool*) malloc(MSA->ncol * MSA->ncol * sizeof(bool));

}

void ProblemBase::FreeBase() {

	free(we);
	free(msa);

}

size_t ProblemBase::GetDim() {

	return dim;

}

void ProblemBase::SetLsingle(double l) {

	lsingle = l;

}

void ProblemBase::SetLpair(double l) {

	lpair = l;

}

void ProblemBase::Iterate() {

	if (ipar_iter != ipar.end()) {
		ipar_iter++;
	}

}

void ProblemBase::SetUpIterations(const std::vector<double> &ipar_) {

	ipar = ipar_;
	ipar_iter = ipar.begin();

}
