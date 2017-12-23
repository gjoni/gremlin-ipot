/*
 * LogRegCoeff.h
 *
 *  Created on: Dec 4, 2017
 *      Author: ivan
 */

#ifndef LOGREGCOEFF_H_
#define LOGREGCOEFF_H_

#include "ContactList.h"

/*
 * training set:      2,370 proteins from the Exceptions set
 * distance cut-off:  5A
 * minimal Neff:      50
 * slogreg param A:   -5
 */

// Length (0:100)
const LogRegCoeff coef5A_range1 = { -8.34894, { -4.18859, -0.01950, -7.57406,
		-0.58222, 1.34039, 1.38801 }, { { 0.00000, 0.12934, 5.95199, 0.13022,
		0.19272, 0.92875 }, { 0.12934, 0.00000, -1.94256, 0.07522, 0.02560,
		-0.14146 }, { 5.95199, -1.94256, 0.00000, -0.31318, 0.68424, 0.70356 },
		{ 0.13022, 0.07522, -0.31318, 0.00000, 0.29182, -0.13607 }, { 0.19272,
				0.02560, 0.68424, 0.29182, 0.00000, -0.29263 }, { 0.92875,
				-0.14146, 0.70356, -0.13607, -0.29263, 0.00000 } } };

// Length [100:150)
const LogRegCoeff coef5A_range2 = { -2.83910, { 0.22918, -1.18167, -3.80202,
		-1.18689, 1.18300, 0.25439 }, { { 0.00000, 0.27196, 5.57276, 0.07692,
		0.08387, -0.10410 }, { 0.27196, 0.00000, -1.97408, 0.05480, 0.04575,
		0.08727 }, { 5.57276, -1.97408, 0.00000, -0.23118, 0.65350, 0.00582 }, {
		0.07692, 0.05480, -0.23118, 0.00000, 0.27211, 0.04748 }, { 0.08387,
		0.04575, 0.65350, 0.27211, 0.00000, -0.26797 }, { -0.10410, 0.08727,
		0.00582, 0.04748, -0.26797, 0.00000 } } };

// Length [150:200)
const LogRegCoeff coef5A_range3 = { 5.54992, { -1.88147, 0.28272, -7.51684,
		-1.05569, -0.59154, -1.44362 }, { { 0.00000, 0.34648, 6.22776, -0.06698,
		0.16898, 0.18019 }, { 0.34648, 0.00000, -1.96316, 0.06303, 0.03479,
		-0.19594 }, { 6.22776, -1.96316, 0.00000, -0.30554, 0.71075, 0.63813 },
		{ -0.06698, 0.06303, -0.30554, 0.00000, 0.25732, 0.04130 }, { 0.16898,
				0.03479, 0.71075, 0.25732, 0.00000, 0.09145 }, { 0.18019,
				-0.19594, 0.63813, 0.04130, 0.09145, 0.00000 } } };

// Length [200:250)
const LogRegCoeff coef5A_range4 = { 9.27634, { -7.31468, -2.25108, -4.69865,
		-1.15624, 0.07800, -2.06678 }, { { 0.00000, 0.40524, 6.39193, -0.10068,
		0.08043, 1.28110 }, { 0.40524, 0.00000, -1.92619, 0.06247, 0.05409,
		0.25441 }, { 6.39193, -1.92619, 0.00000, -0.31852, 0.70867, 0.10849 }, {
		-0.10068, 0.06247, -0.31852, 0.00000, 0.24550, 0.07446 }, { 0.08043,
		0.05409, 0.70867, 0.24550, 0.00000, -0.04368 }, { 1.28110, 0.25441,
		0.10849, 0.07446, -0.04368, 0.00000 } } };

// Length [250:300)
const LogRegCoeff coef5A_range5 = { -16.22329, { 3.98167, -1.77373, -7.71974,
		-1.48885, 2.87187, 2.54773 }, { { 0.00000, 0.45589, 6.85143, -0.22942,
		0.04615, -0.73450 }, { 0.45589, 0.00000, -2.12832, 0.08618, 0.06525,
		0.13975 }, { 6.85143, -2.12832, 0.00000, -0.44696, 0.82626, 0.57613 }, {
		-0.22942, 0.08618, -0.44696, 0.00000, 0.23715, 0.13281 }, { 0.04615,
		0.06525, 0.82626, 0.23715, 0.00000, -0.54439 }, { -0.73450, 0.13975,
		0.57613, 0.13281, -0.54439, 0.00000 } } };

// Length [300:400)
const LogRegCoeff coef5A_range6 = { -1.39055, { 2.19482, -0.67260, -10.54219,
		-0.42716, 0.47068, -0.10821 }, { { 0.00000, 0.49827, 6.55125, -0.18701,
		0.02697, -0.41330 }, { 0.49827, 0.00000, -2.15669, 0.06907, 0.06509,
		-0.05106 }, { 6.55125, -2.15669, 0.00000, -0.39296, 0.85304, 1.01932 },
		{ -0.18701, 0.06907, -0.39296, 0.00000, 0.22120, -0.02903 }, { 0.02697,
				0.06509, 0.85304, 0.22120, 0.00000, -0.11742 }, { -0.41330,
				-0.05106, 1.01932, -0.02903, -0.11742, 0.00000 } } };

// Length [400:inf)
const LogRegCoeff coef5A_range7 = { 9.97437, { -0.10094, -1.06926, -8.28565,
		-0.35061, -1.09885, -2.02961 }, { { 0.00000, 0.33470, 5.00737, -0.19642,
		-0.13216, 0.27140 }, { 0.33470, 0.00000, -1.90580, 0.07848, 0.08667,
		-0.00758 }, { 5.00737, -1.90580, 0.00000, -0.35574, 0.75205, 0.65325 },
		{ -0.19642, 0.07848, -0.35574, 0.00000, 0.20668, -0.02567 }, { -0.13216,
				0.08667, 0.75205, 0.20668, 0.00000, 0.14498 }, { 0.27140,
				-0.00758, 0.65325, -0.02567, 0.14498, 0.00000 } } };

/* END 5A PARAMS */

/*
 * training set:      2,370 proteins from the Exceptions set
 * distance cut-off:  8A
 * minimal Neff:      50
 * slogreg param A:   -9
 */

// Length (0:100)
const LogRegCoeff coef8A_range1 = { -1.20781, { -0.68459, -0.58463, -2.58467,
		-1.30928, 0.79473, 0.00286 }, { { 0.00000, -0.09115, 5.38028, 0.36944,
		0.19022, 0.20543 }, { -0.09115, 0.00000, -1.49152, 0.00296, 0.02441,
		0.02097 }, { 5.38028, -1.49152, 0.00000, -0.32777, 0.43708, -0.04944 },
		{ 0.36944, 0.00296, -0.32777, 0.00000, 0.27380, 0.03179 }, { 0.19022,
				0.02441, 0.43708, 0.27380, 0.00000, -0.15943 }, { 0.20543,
				0.02097, -0.04944, 0.03179, -0.15943, 0.00000 } } };

// Length [100:150)
const LogRegCoeff coef8A_range2 = { 0.42317, { -0.58240, -0.95674, -1.53914,
		-1.42116, 0.88362, -0.30439 }, { { 0.00000, -0.00457, 5.37736, 0.42367,
		0.09867, 0.18532 }, { -0.00457, 0.00000, -1.59354, -0.00567, 0.03159,
		0.08460 }, { 5.37736, -1.59354, 0.00000, -0.28225, 0.48357, -0.20202 },
		{ 0.42367, -0.00567, -0.28225, 0.00000, 0.31442, 0.01443 }, { 0.09867,
				0.03159, 0.48357, 0.31442, 0.00000, -0.18129 }, { 0.18532,
				0.08460, -0.20202, 0.01443, -0.18129, 0.00000 } } };

// Length [150:200)
const LogRegCoeff coef8A_range3 = { 6.28020, { 0.15388, -0.26860, -3.93880,
		-1.49475, -0.34326, -1.43642 }, { { 0.00000, 0.01471, 5.49701, 0.33398,
		0.06676, 0.06329 }, { 0.01471, 0.00000, -1.53634, 0.00969, 0.03812,
		-0.06716 }, { 5.49701, -1.53634, 0.00000, -0.26932, 0.49420, 0.25853 },
		{ 0.33398, 0.00969, -0.26932, 0.00000, 0.31748, 0.01805 }, { 0.06676,
				0.03812, 0.49420, 0.31748, 0.00000, 0.06068 }, { 0.06329,
				-0.06716, 0.25853, 0.01805, 0.06068, 0.00000 } } };

// Length [200:250)
const LogRegCoeff coef8A_range4 = { 8.62070, { -8.10800, -1.18919, -6.41102,
		-0.44913, -0.07203, -1.83047 }, { { 0.00000, 0.10326, 5.60631, 0.25804,
		0.01469, 1.61896 }, { 0.10326, 0.00000, -1.53347, 0.01697, 0.04637,
		0.08989 }, { 5.60631, -1.53347, 0.00000, -0.26442, 0.46432, 0.74252 }, {
		0.25804, 0.01697, -0.26442, 0.00000, 0.31060, -0.16545 }, { 0.01469,
		0.04637, 0.46432, 0.31060, 0.00000, 0.00576 }, { 1.61896, 0.08989,
		0.74252, -0.16545, 0.00576, 0.00000 } } };

// Length [250:300)
const LogRegCoeff coef8A_range5 = { -14.28196, { 6.65466, -1.57913, -3.88611,
		-1.58241, 2.72565, 2.32708 }, { { 0.00000, 0.11820, 5.91026, 0.18592,
		0.01557, -1.07084 }, { 0.11820, 0.00000, -1.67135, 0.04120, 0.05164,
		0.14947 }, { 5.91026, -1.67135, 0.00000, -0.30929, 0.54082, 0.23979 }, {
		0.18592, 0.04120, -0.30929, 0.00000, 0.31956, 0.02045 }, { 0.01557,
		0.05164, 0.54082, 0.31956, 0.00000, -0.50142 }, { -1.07084, 0.14947,
		0.23979, 0.02045, -0.50142, 0.00000 } } };

// Length [300:400)
const LogRegCoeff coef8A_range6 = { -1.27629, { 3.95751, -0.64785, -6.89043,
		-0.55736, 0.66091, -0.00651 }, { { 0.00000, 0.18392, 5.01536, 0.10659,
		-0.03807, -0.54616 }, { 0.18392, 0.00000, -1.52856, 0.03475, 0.05644,
		-0.02133 }, { 5.01536, -1.52856, 0.00000, -0.28094, 0.51250, 0.72881 },
		{ 0.10659, 0.03475, -0.28094, 0.00000, 0.29407, -0.11836 }, { -0.03807,
				0.05644, 0.51250, 0.29407, 0.00000, -0.13312 }, { -0.54616,
				-0.02133, 0.72881, -0.11836, -0.13312, 0.00000 } } };

// Length [400:inf)
const LogRegCoeff coef8A_range7 = { 6.69460, { -1.93841, -0.82142, -6.58281,
		-0.40931, -0.30091, -1.36786 }, { { 0.00000, 0.09226, 4.55330, 0.02898,
		-0.08300, 0.55876 }, { 0.09226, 0.00000, -1.47511, 0.03611, 0.06416,
		0.00212 }, { 4.55330, -1.47511, 0.00000, -0.29377, 0.48138, 0.68273 }, {
		0.02898, 0.03611, -0.29377, 0.00000, 0.25020, -0.08216 }, { -0.08300,
		0.06416, 0.48138, 0.25020, 0.00000, 0.02868 }, { 0.55876, 0.00212,
		0.68273, -0.08216, 0.02868, 0.00000 } } };

/* END 8A PARAMS */

#endif /* LOGREGCOEFF_H_ */
