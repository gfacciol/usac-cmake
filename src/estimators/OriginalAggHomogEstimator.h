//
// Created by rais on 15/03/17.
//

#ifndef USAC_ORIGINALAGGHOMOESTIMATOR_H
#define USAC_ORIGINALAGGHOMOESTIMATOR_H

#include "AggHomogEstimator.h"
#include "../utils/MathFunctions.h"


class OriginalAggHomogEstimator : public AggHomogEstimator
{
public:
    inline void      aggregate_weighted_mean(std::vector<double *> &projected_points, double *res);
    inline void      aggregate_weighted_geometric_median(std::vector<double *> &projected_points, double *wmean, double *res);
    inline void      perform_aggregation();
    OriginalAggHomogEstimator(unsigned int im1w, unsigned int im1h, unsigned int im2w, unsigned int im2h, double pmean, double pwgmed) : p_mean(pmean), p_wgmed(pwgmed)
    {
        srcPoints.clear();
        srcPoints.push_back(new double[3]{0, 0, 1});
        srcPoints.push_back(new double[3]{(double) im1w, 0, 1});
        srcPoints.push_back(new double[3]{0,(double) im1h, 1});
        srcPoints.push_back(new double[3]{(double) im1w,(double) im1h, 1});
    }
    ~OriginalAggHomogEstimator()
    {
        for (size_t i = 0; i < srcPoints.size(); ++i)
        {
            if (srcPoints[i]) { delete[] srcPoints[i]; }
        }
        srcPoints.clear();
    }
protected:
    std::vector<double *> srcPoints;
    double p_mean;
    double p_wgmed;
protected:
    inline void     normalize_weighted_scores(std::vector<unsigned int> &scores, double p, std::vector<double> &weights);

};

void OriginalAggHomogEstimator::normalize_weighted_scores(std::vector<unsigned int> &scores,
                                                          double p, std::vector<double> &weights)
{
    double sum = 0;
    for (int i=0;i<scores.size();i++)
        sum += pow((double) scores[i], p);
    for (int i=0;i<scores.size();i++)
        weights[i] = pow((double) scores[i], p)/sum;

}
void OriginalAggHomogEstimator::aggregate_weighted_geometric_median(std::vector<double *> &projected_points, double *wmean, double *res)
{
    //! First, normalize the scores
    std::vector<double> weights(intermediate_models.size());
    normalize_weighted_scores(intermediate_models_scores, p_wgmed, weights);

    double tol = 1.e-7;
    int max_iter = 500;

    double y[2], yold[2];
    std::copy(wmean, wmean + 2, y);

    //! Check amount of elements
    if (projected_points.size() > 1) {
        double fTiny = 0.0000000001;
        for(int i=0;i<max_iter;i++) {
            double num[2]={0,0}, denom=0;
            for (int j=0;j<projected_points.size();j++)
            {
                double val = MathTools::l2_dist(projected_points[j], y, 2);
                val = (val > fTiny ? val : fTiny);
                denom += weights[j] / val;
                for (int k=0;k<2;k++) num[k] += weights[j] * projected_points[j][k] / val;
            }
            for (int k=0;k<2;k++)
            {
                yold[k] = y[k];
                y[k] = num[k] / denom;

            }
            if (MathTools::l2_dist(yold, y, 2) / MathTools::norm(y, 2) < tol)
                break;
        }
        res[0] = y[0]; res[1] = y[1];
    }
    else {
        if (projected_points.size() == 1) {
            res[0] = projected_points[0][0];
            res[1] = projected_points[0][1];
            res[2] = projected_points[0][2];
        }
    }
}

void OriginalAggHomogEstimator::aggregate_weighted_mean(std::vector<double *> &projected_points, double *res)
{
    //! First, normalize the scores
    std::vector<double> weights(intermediate_models.size());
    normalize_weighted_scores(intermediate_models_scores, p_wgmed, weights);
    res[0] = res[1] = 0;
    res[2] = 1;
    for (int i=0;i<projected_points.size();i++)
    {
        double *pt = projected_points[i];
        res[0] += pt[0] * weights[i];
        res[1] += pt[1] * weights[i];
    }
}

void OriginalAggHomogEstimator::perform_aggregation()
{
    if (intermediate_models.size() > 0) {
        double *wmean_agg_points = new double[6 * srcPoints.size()];
        double *wgmed_agg_points = new double[6 * srcPoints.size()];
        //! For each predefined corner point
        for (int j = 0; j < srcPoints.size(); j++) {
            double *pt = srcPoints[j];
            // Project all hypothesized homographies
            std::vector<double *> projected_points;
            projected_points.clear();
            projected_points.resize(intermediate_models.size());
            for (int i = 0; i < intermediate_models.size(); i++) {
                double *model = intermediate_models[i];
                double h_x[3];
                //! Project
                MathTools::vmul(h_x, model, pt, 3);
                //! Normalize to 2D coordinates
                h_x[0] /= h_x[2];
                h_x[1] /= h_x[2];
                h_x[2] = 1;
                //! Store
                projected_points[i] = new double[2]{h_x[0], h_x[1]};
                //printf("Pt: (%f,%f) -> (%f,%f). Weight: %d\n", pt[0], pt[1], h_x[0], h_x[1], intermediate_models_scores[i]);
            }
            // Now aggregate all projections
            double res_wmean[3];
            for (int i = 0; i < 3; i++) res_wmean[i] /= res_wmean[2];
            aggregate_weighted_mean(projected_points, res_wmean);
            wmean_agg_points[0 + 6 * j] = pt[0] / pt[2];
            wmean_agg_points[1 + 6 * j] = pt[1] / pt[2];
            wmean_agg_points[2 + 6 * j] = 1;
            wmean_agg_points[3 + 6 * j] = res_wmean[0];
            wmean_agg_points[4 + 6 * j] = res_wmean[1];
            wmean_agg_points[5 + 6 * j] = 1;
            double res_wgmed[3];

            aggregate_weighted_geometric_median(projected_points, res_wmean, res_wgmed);
            wgmed_agg_points[0 + 6 * j] = pt[0] / pt[2];
            wgmed_agg_points[1 + 6 * j] = pt[1] / pt[2];
            wgmed_agg_points[2 + 6 * j] = 1;
            wgmed_agg_points[3 + 6 * j] = res_wgmed[0];
            wgmed_agg_points[4 + 6 * j] = res_wgmed[1];
            wgmed_agg_points[5 + 6 * j] = 1;
            //printf("Aggregation results: wgmed: (%f,%f). wmean: (%f,%f)\n", res_wgmed[0], res_wgmed[1], res_wmean[0], res_wmean[1]);

            //! Clear up stuff
            for (int i = 0; i < projected_points.size(); i++)
                delete projected_points[i];

        }
        min_sample_.clear();
        for (unsigned int i = 0; i < 3; ++i) {
            for (unsigned int j = 0; j < 3; ++j) {
                final_model_params_[3 * i + j] /= final_model_params_[3 * 2 + 2];
            }
        }
        for (int i = 0; i < srcPoints.size(); i++)
            min_sample_[i] = i;
        double *wmean_agg_points_normalized = new double[6 * srcPoints.size()];
        double *wgmed_agg_points_normalized = new double[6 * srcPoints.size()];

        if (data_matrix_) {
            delete[] data_matrix_;
            data_matrix_ = NULL;
        }
        unsigned int usac_num_data_points_old = usac_num_data_points_;
        usac_num_data_points_ = srcPoints.size();
        data_matrix_ = new double[18 * usac_num_data_points_]();    // 2 equations per correspondence
        //! Normalized the aggregated points
        FTools::normalizePoints(wmean_agg_points, wmean_agg_points_normalized, srcPoints.size(), m_T1_, m_T2_);
        for (unsigned int i = 0; i < 9; ++i)
            m_T2inv_[i] = m_T2_[i];
        MathTools::minv(m_T2inv_, 3);
        //! Now compute the data matrix and generate the model
        HTools::computeDataMatrix(data_matrix_, srcPoints.size(), wmean_agg_points_normalized);
        generateMinimalSampleModels();
//    MathTools::minv(models_denorm_[0], 3);
        // save the current model as the solution
        for (unsigned int i = 0; i < 9; ++i)
            aggregated_wmean_model_params_[i] = *(models_denorm_[0] + i) / *(models_denorm_[0] + 8);

        // Do the same for the wgmed aggregation
        if (data_matrix_) {
            delete[] data_matrix_;
            data_matrix_ = NULL;
        }
        data_matrix_ = new double[18 * usac_num_data_points_]();    // 2 equations per correspondence
        //! Normalized the aggregated points
        FTools::normalizePoints(wgmed_agg_points, wgmed_agg_points_normalized, srcPoints.size(), m_T1_, m_T2_);
        for (unsigned int i = 0; i < 9; ++i)
            m_T2inv_[i] = m_T2_[i];
        MathTools::minv(m_T2inv_, 3);
        //! Now compute the data matrix and generate the model
        HTools::computeDataMatrix(data_matrix_, srcPoints.size(), wgmed_agg_points_normalized);
        generateMinimalSampleModels();
//    MathTools::minv(models_denorm_[0], 3);
        // save the current model as the solution
        for (unsigned int i = 0; i < 9; ++i)
            aggregated_wgmed_model_params_[i] = *(models_denorm_[0] + i) / *(models_denorm_[0] + 8);
        usac_num_data_points_ = usac_num_data_points_old;

        delete[] wmean_agg_points;
        delete[] wmean_agg_points_normalized;
        delete[] wgmed_agg_points;
        delete[] wgmed_agg_points_normalized;
    } else
    {
        // Copy RANSAC results
        for (unsigned int i = 0; i < 9; ++i) {
            aggregated_wgmed_model_params_[i] = final_model_params_[i];
            aggregated_wmean_model_params_[i] = final_model_params_[i];
        }

    }
}

#endif //USAC_ORIGINALAGGHOMOESTIMATOR_H
