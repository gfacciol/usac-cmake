#include "HomographyFunctions.h"

namespace HTools
{

	double computeError(std::vector<double> &gt_pts, double *model)
	{
		double acumErr1 = 0, acumErr2 = 0;
		int qtyPts = gt_pts.size() / 6;
		for (unsigned int i=0;i<qtyPts; i++)
		{
			double proj[3]; double pt1[3]; double pt2[3];
			pt1[0] = gt_pts[0+i*6]; pt1[1] = gt_pts[1+i*6]; pt1[2] = gt_pts[2+i*6];
			pt2[0] = gt_pts[3+i*6]; pt2[1] = gt_pts[4+i*6]; pt2[2] = gt_pts[5+i*6];
			//! Project point
			MathTools::vmul(proj, model, pt1, 3);
			//! Normalize
			proj[0] /= proj[2]; proj[1] /= proj[2];
			//! Compute error
			acumErr1 += pow((proj[0] - pt2[0]),2) + pow((proj[1] - pt2[1]),2);
			//! Compute the error by projecting back
			double inv_model[9];
			for (unsigned int j = 0; j < 9; ++j)
				inv_model[j] = model[j];
			MathTools::minv(inv_model, 3);
			for (int k=0;k<9;k++)
				inv_model[k] /= inv_model[8];
			//! Project point
			MathTools::vmul(proj, inv_model, pt2, 3);
			//! Normalize
			proj[0] /= proj[2]; proj[1] /= proj[2];
			//! Compute error
			acumErr2 += pow((proj[0] - pt1[0]),2) + pow((proj[1] - pt1[1]),2);
		}
		return ((acumErr1 / (double) gt_pts.size()) + (acumErr2 / (double) gt_pts.size())) / 2;
	}

	void computeDataMatrix(double* data_matrix, unsigned int num_points, double* points)
	{
		// linearizes corresp. with respect to entries of homography matrix,
		// so that u' = H u -> A h 

		const double *data_ptr;
		double *p;
		unsigned int offset = 2*num_points;

		for (unsigned int i = 0; i < num_points; ++i)
		{
			data_ptr = points + 6*i;
			p = data_matrix + 2*i;

			*p				= 0;
			*(p + offset)	= 0;
			*(p + 2*offset) = 0;
			*(p + 3*offset) = -data_ptr[0];
			*(p + 4*offset) = -data_ptr[1];
			*(p + 5*offset) = -data_ptr[2];
			*(p + 6*offset) = data_ptr[4] * data_ptr[0];
			*(p + 7*offset) = data_ptr[4] * data_ptr[1];
			*(p + 8*offset) = data_ptr[4] * data_ptr[2];

			p = data_matrix + 2*i + 1;

			*p				= data_ptr[0];
			*(p + offset)	= data_ptr[1];
			*(p + 2*offset) = data_ptr[2];
			*(p + 3*offset) = 0;
			*(p + 4*offset) = 0;
			*(p + 5*offset) = 0;
			*(p + 6*offset) = -data_ptr[3] * data_ptr[0];
			*(p + 7*offset) = -data_ptr[3] * data_ptr[1];
			*(p + 8*offset) = -data_ptr[3] * data_ptr[2];
		}
	} // end computeDataMatrix

	void crossprod(double *out, const double *a, const double *b, unsigned int st)
	{
		unsigned int st2 = 2 * st;
		out[0] = a[st]*b[st2] - a[st2]*b[st];
		out[1] = a[st2]*b[0]  - a[0]*b[st2];
		out[2] = a[0]*b[st]   - a[st]*b[0];
	}

}