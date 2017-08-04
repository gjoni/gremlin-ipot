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

	friend class ProblemBase;
	friend class ProblemFull;
	friend class MRFclass;

private:

	/* raw MSA in A3M/FASTA */
	std::vector<std::pair<std::string, std::string> > a3m;

	/* reference sequence length (1st sequence in the MSA) */
	size_t len_ref;

	size_t CleanLowercase(std::string &str); /* in-place removal of lowercase letters */

	/* cleaned rows and columns */
	std::vector<size_t> row_map; /* cleaned to reference */
	std::vector<size_t> col_map; /* cleaned to reference */
	std::vector<size_t> a3m_to_msa;

	size_t CleanRows(double gaps_frac); /* sets row_map */
	size_t CleanCols(double gaps_frac); /* sets col_map, ref_to_cleaned */

	/* size of the CLEANED alignment msa[] */
	size_t nrow; // sequence length - 65536 MAX
	size_t ncol; // number of sequences in MSA - 65536 MAX

	/* alloc/free msa*/
	void Allocate();
	void Free();

	void TrimRight(char *str);

public:

	/*
	 * conversions between AA letters and indices
	 */
	static const unsigned char NAA = 21;
	static const unsigned char AMINO_INDICES[26];
	static const unsigned char CHAR_INDICES[NAA];
	static unsigned char aatoi(unsigned char aa);
	static unsigned char itoaa(unsigned char aa);
	static void aatoi(unsigned char *str, size_t len);

	MSAclass(const char *name);
	MSAclass(const MSAclass &source);
	MSAclass();

	~MSAclass();

	MSAclass& operator=(const MSAclass &source);

	/* cleaned MSA dimensions */
	size_t GetNrow() const;
	size_t GetNcol() const;
	size_t GetLen() const;

	/* A3M->MSA mapping*/
	size_t GetMsaIdx(size_t idx) const;

	/* change residue indices in the array of contacts
	 * to match the cleaned msa[]
	 * (contacts with residues missing in the msa[] are omitted) */
	std::vector<std::pair<size_t, size_t> > CastToMsa(
			const std::vector<std::pair<int, int> > &contacts) const;

	/* for a complex */
//	std::vector<std::pair<size_t, size_t> > CastToMsa(
//			const std::vector<std::pair<int, int> > &contacts,
//			size_t len_rec) const;
	char GetA3Mres(size_t i, size_t j) const;
//	char GetMSAres(size_t i, size_t j) const;

	/* rows are cleaned first, then columns */
	void CleanMsa(double rgaps, double cgaps);

	/* 1d array [nrow x ncol] with the cleaned alignment */
	unsigned char * GetMsa() const;

	void SaveMSA(const std::string &name) const;
	void PrintMSA() const;

	/* TODO: entropies */

	/* TODO: list all possible pairs of continuous tuples of length len */
//	std::vector<std::pair<size_t, size_t> > GetTuples(const MSAclass &MSA,
//			size_t len);

};

#endif /* MSACLASS_H_ */
