/*
 * CMap.h
 *
 *  Created on: Jan 26, 2018
 *      Author: ivan
 */

#ifndef CMAP_H_
#define CMAP_H_

#include <vector>
#include <utility>
#include <tuple>
#include <string>

/* neighbors list type: element -> (neighbor; contact_score) */
//typedef std::vector<std::pair<unsigned, double> > NListT;

/* neighbors list type: element -> (neighbor; contact_score; separation) */
typedef std::vector<std::tuple<unsigned, double, unsigned> > NListT;

/* adjacency list type */
typedef std::vector<NListT> AListT;

/*
 * contact map class
 */
class CMap {
private:

	/* dimension (aka sequence length) */
	unsigned size;

	/* neighbors to the left and to the right
	 * of the diagonal */
	AListT left, right;

	/* mappings to positions with nonempty
	 * left/right neighbor lists */
	std::vector<unsigned> mleft, mright;

	/* contact map (2D array in 1D) */
	/* TODO: don't really need this vector */
//	std::vector<double> mtx1d;

public:

	CMap();

	/* from file */
	CMap(const char*);
	CMap(const std::string&);

	/* from adjacency list */
	CMap(const AListT&);

	/* from Chain */
//	CMap(const Chain&);

	~CMap();

	double & operator()(unsigned i, unsigned j);

	const NListT& GetLeftList(unsigned i);
	const NListT& GetRightList(unsigned i);

};

#endif /* CMAP_H_ */
