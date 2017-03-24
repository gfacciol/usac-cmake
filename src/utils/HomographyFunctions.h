#ifndef HTOOLS_H
#define HTOOLS_H

#include "MathFunctions.h"
#include <vector>

namespace HTools
{
	void computeDataMatrix(double* data_matrix, unsigned int num_points, double* points);
	void crossprod(double *out, const double *a, const double *b, unsigned int st);
	double computeError(std::vector<double> &gt_pts, double *model);
}
#endif