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

using namespace std;

class MapAlign {
private:

	MapAlign();
	~MapAlign();

	struct SWDATA {
		unsigned M, N; /* dimensions */
		double **mtx; /* contacts scoring matrix (M x N) */
		double **sco; /* DP scoring matrix (M + 1) x (N + 1) */
		char **label; /* path in the DP matrix (M + 1) x (N + 1) */
		std::vector<double> gap_a; /* gap opening penalties for A */
		std::vector<double> gap_b; /* gap opening penalties for B */
		std::vector<int> a2b;
		std::vector<int> b2a;
	};

	static void Alloc(SWDATA*);
	static void Free(SWDATA*);

	static double gaussian(double mean, double std, double x);
	static double sepw(double sep);

	/* step1 Smith-Waterman alignment to calculate score only */
	static double SW1(const NListT&, const NListT&, double, double);

	/* step2 Smith-Waterman alignment */
	static double SW2(SWDATA&, double gap_e);

	/* TODO: current implementation is not very efficient */
	static double Intersect(const NListT&, const NListT&,
			const std::vector<int>&, const std::vector<int>&);

	/* a function to assess current alignment
	 * based on contact/gap scores - returned as a vector */
	static std::vector<double> Assess(const SWDATA&, const CMap&, const CMap&,
			double);

	static double GapScore(const std::vector<int>&, const std::vector<double>&,
			double);

	static void InitMTX(SWDATA&, const CMap&, const CMap&, double sep_x,
			double sep_y);

	static void UpdateMTX(SWDATA&, const CMap&, const CMap&, double, int iter);

	static void CheckMTX();

public:

	struct PARAMS {
		double gap_open;
		double gap_ext;
		int kmin;
		int iter;
	};

	static double MaxScore(const CMap&);

	static double Align(const CMap&, const CMap&, const PARAMS&,
			std::vector<int>&);

};

#endif /* MAPALIGN_H_ */
