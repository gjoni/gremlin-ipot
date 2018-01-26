/*
 * MapAlign.h
 *
 *  Created on: Jan 26, 2018
 *      Author: ivan
 */

#ifndef MAPALIGN_H_
#define MAPALIGN_H_

#include <vector>

#include "CMap.h"

class MapAlign {
private:

	MapAlign();
	~MapAlign();

	struct SWDATA {
		unsigned M, N; /* dimensions */
		double **mtx; /* contacts scoring matrix (M x N) */
		double **sco; /* DP scoring matrix (M + 1) x (N + 1) */
		char **label; /* path in the DP matrix (M + 1) x (N + 1) */
		std::vector<double>& gap_a; /* gap opening penalties for A */
		std::vector<double>& gap_b; /* gap opening penalties for B */
	};

	static void Alloc(SWDATA*);
	static void Free(SWDATA*);

	/* Smith-Waterman alignment to calculate score only */
	static double SWsimple(const NListT&, const NListT&);

	/* Smith-Waterman alignment */
	static double SW(SWDATA&, CMap&, CMap&, double gap_e,
			std::vector<double>& a2b);

	static void InitMTX(SWDATA&, CMap&, CMap&, double sep_x, double sep_y);
	static void UpdateMTX(SWDATA&, CMap&, CMap&, int it);
	static void CheckMTX();

public:

	static double Align(const CMap&, const CMap&, std::vector<unsigned>&);

};

#endif /* MAPALIGN_H_ */
