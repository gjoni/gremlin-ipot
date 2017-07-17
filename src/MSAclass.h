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
 * 2) cleaning should be incorporated (by masking or weighting ???):
 *     - columns cleaning
 *     - rows cleaning ???
 * 3) weighting - does it belong to this class?
 */

/*
 * Specs:
 *  - read A3M/FASTA, remove lowercase letters
 *  - clean MSA and prepare it for GREMLIN
 *  - append ???
 *  - concat ???
 */

class MSAclass {

private:

	/* raw MSA in A3M/FASTA */
	std::vector<std::pair<std::string, std::string> > a3m;

	/* reference sequence length (1st sequence in the MSA) */
	size_t len_ref;

	size_t CleanLowercase(std::string &str); /* in-place removal of lowercase letters */

	/* cleaned rows and columns */
	std::vector<size_t> row_map;
	std::vector<size_t> col_map;

	size_t CleanRows(double gaps_frac); /* sets row_map */
	size_t CleanCols(double gaps_frac); /* sets col_map */
	void SetMsa();

	/* 1D representation of the CLEANED alignment
	 * msa[row_map.size() * col_map.size()] */
	unsigned char *msa;

	/* size of the CLEANED alignment msa[] */
	uint16_t nrow; // sequence length - 65536 MAX
	uint16_t ncol; // number of sequences in MSA - 65536 MAX

	/* alloc/free msa*/
	void Allocate();
	void Free();

	void TrimRight(char *str);

public:

	/*
	 * conversions between AA letters and indices
	 */
	static const unsigned char AMINO_INDICES[26];
	static const unsigned char CHAR_INDICES[21];
	static unsigned char aatoi(unsigned char aa);
	static unsigned char itoaa(unsigned char aa);
	static void aatoi(unsigned char *str, size_t len);

	MSAclass(const char *name);
	MSAclass(const MSAclass &source);
	MSAclass();

	~MSAclass();

	MSAclass& operator=(const MSAclass &source);

	uint16_t GetNrow();
	uint16_t GetNcol();

	char GetResidue(uint16_t i, uint16_t j);

	/* rows are cleaned first, then columns */
	void CleanMsa(double rgaps, double cgaps);

	/* msa cannot be modified outside the class */
	unsigned char const * GetMsa();

	/* convert msa[] characters into indices {0..20} */
	void CastToIdx();

};

#endif /* MSACLASS_H_ */
