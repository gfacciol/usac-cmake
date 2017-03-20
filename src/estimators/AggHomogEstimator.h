//
// Created by rais on 15/03/17.
//

#ifndef USAC_AGGHOMOGESTIMATOR_H_H
#define USAC_AGGHOMOGESTIMATOR_H_H

#include "HomogEstimator.h"

class AggHomogEstimator: public HomogEstimator
{


public:
    // storage for the intermediate hypotheses
    std::vector<double> aggregated_wgmed_model_params_;
    std::vector<double> aggregated_wmean_model_params_;

public:
    AggHomogEstimator() : HomogEstimator::HomogEstimator()
    {
        aggregated_wgmed_model_params_.clear();
        aggregated_wgmed_model_params_.resize(9);
        aggregated_wmean_model_params_.clear();
        aggregated_wmean_model_params_.resize(9);
    };
    ~AggHomogEstimator()
    {
    };

protected:
    // ------------------------------------------------------------------------
    // problem specific functions
    inline virtual void      perform_aggregation() = 0;

};

#endif //USAC_AGGHOMOGESTIMATOR_H_H
