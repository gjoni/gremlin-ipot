/*
 * MSAclass.h
 *
 *  Created on: Jul 8, 2016
 *      Author: ivan
 */

#ifndef MSACLASS_H_
#define MSACLASS_H_

#include <string>
#include <vector>

/*
 * TODO:
 * ----
 * 1) this class should be able to read A3M files directly
 * 2) cleaning should be incorporated (by masking or weighting ???)
 * 3) weighting - does it belong to this class?
 */

class MSAclass {

private:

	char **msa; // L x N
	size_t nrow; // sequence length
	size_t ncol; // number of sequences in MSA

	void Allocate();
	void Free();

	void TrimRight(char *str);

public:

	static const unsigned char AMINO_INDICES[26];
	static const unsigned char CHAR_INDICES[21];

	MSAclass(const char *name);
	MSAclass(const MSAclass &source);
	MSAclass();

	~MSAclass();

	MSAclass& operator=(const MSAclass &source);

	size_t GetNrow();
	size_t GetNcol();

	char GetResidue(size_t i, size_t j);

	static int aatoi(char aa);

};

#endif /* MSACLASS_H_ */
