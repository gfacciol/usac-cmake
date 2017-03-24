//
// Created by rais on 20/03/17.
//

#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>
#endif

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <vector>
#include "estimators/OriginalAggHomogEstimator.h"
#include "config/ConfigParams.h"
#include "estimators/FundmatrixEstimator.h"
#include "estimators/HomogEstimator.h"

extern "C"
    {
    #include "iio/iio.h"
    #include "sift/lib_sift.h"
}

#define AMOUNT_OF_TESTS 10001

// helper functions
bool readGTFromFile(char *inputFilePath, std::vector<double>& pointData, std::vector<double>& H)
{
    // read data from from file
    std::ifstream infile(inputFilePath);
    if (!infile.is_open())
    {
        std::cerr << "Error opening input points file: " << inputFilePath << std::endl;
        return false;
    }
    H.resize(9);
    for (int i=0;i<9;i++)
        infile >> H[i];
    unsigned int numPts;
    infile >> numPts;
    pointData.resize(6*numPts);
    for (unsigned int i = 0; i < numPts; ++i)
    {
        for (unsigned int j=0;j<6;j++)
            infile >> pointData[6*i+j];
    }
    infile.close();
    return true;
}

bool readCorrsFromFile(std::string& inputFilePath, std::vector<double>& pointData, unsigned int& numPts)
{
    // read data from from file
    std::ifstream infile(inputFilePath.c_str());
    if (!infile.is_open())
    {
        std::cerr << "Error opening input points file: " << inputFilePath << std::endl;
        return false;
    }
    infile >> numPts;
    pointData.resize(6*numPts);
    for (unsigned int i = 0; i < numPts; ++i)
    {
        infile >> pointData[6*i] >> pointData[6*i+1] >> pointData[6*i+3] >> pointData[6*i+4];
        pointData[6*i+2] = 1.0;
        pointData[6*i+5] = 1.0;
    }
    infile.close();
    return true;
}

bool readPROSACDataFromFile(std::string& sortedPointsFile, unsigned int numPts, std::vector<unsigned int>& prosacData)
{
    std::ifstream infile(sortedPointsFile.c_str());
    if (!infile.is_open())
    {
        std::cerr << "Error opening sorted indices file: " << sortedPointsFile << std::endl;
        return false;
    }
    for (unsigned int i = 0; i < numPts; ++i)
    {
        infile >> prosacData[i];
    }
    infile.close();
    return true;
}


// ---------------------------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
    // check command line args
    if (argc < 4)
    {
        std::cerr << "Usage: RansaacDemo <file1> <file2> <config file> <gt_file>" << std::endl;
        std::cerr << "\t<config file>: full path to configuration file" << std::endl;
        std::cerr << "\t<gt_file>: full path to ground truth file" << std::endl;
        return(EXIT_FAILURE);
    }
    std::string image1_path = argv[1];
    std::string image2_path = argv[2];
    std::string cfg_file_path = argv[3];

    //! Read GT if availables
    std::vector<double> gt_pts, gt_H;
    if (argc >= 5)
        readGTFromFile(argv[4], gt_pts, gt_H);
    //! Open image files
    // Loading image
    printf("Reading image: %s.\n", argv[1]);
    int w1, h1;
    float *im1 = iio_read_image_float(argv[1], &w1, &h1);
    for(int i=0; i < w1*h1; i++)
        im1[i] /=256.;

    printf("Reading image: %s.\n", argv[2]);
    int w2, h2;
    float *im2 = iio_read_image_float(argv[2], &w2, &h2);
    for(int i=0; i < w2*h2; i++)
        im2[i] /=256.;

    // compute sift keypoints
    struct sift_keypoints *matches_1, *matches_2;
    int qty_of_matches = sift_compute_matches(im1, w1, h1, im2, w2, h2, 0.6, &matches_1, &matches_2, 1);
    std::vector<double> matches(matches_1->size*6);
    for (int i = 0; i < matches_1->size;i++)
    {
        matches[i*6] = (double) matches_1->list[i]->y;
        matches[i*6+1] = (double) matches_1->list[i]->x;
        matches[i*6+2] = 1.0;
        matches[i*6+3] = (double) matches_2->list[i]->y;
        matches[i*6+4] = (double) matches_2->list[i]->x;
        matches[i*6+5] = 1.0;
    }
    free(im1);
    free(im2);
    // seed random number generator
    srand((unsigned int)time(NULL));

    // ------------------------------------------------------------------------
    // initialize the homography estimation problem
    ConfigParamsHomog cfg;
    if ( !cfg.initParamsFromConfigFile(cfg_file_path) )
    {
        std::cerr << "Error during initialization" << std::endl;
        return(EXIT_FAILURE);
    }
//    std::vector<double> point_data;
//    if ( !readCorrsFromFile(cfg.homog.inputFilePath, point_data, cfg.common.numDataPoints) )
//    {
//        return(EXIT_FAILURE);
//    }

    std::vector<unsigned int> prosac_data;
    if (cfg.common.randomSamplingMethod == USACConfig::SAMP_PROSAC)
    {
        prosac_data.resize(cfg.common.numDataPoints);
        if ( !readPROSACDataFromFile(cfg.prosac.sortedPointsFile, cfg.common.numDataPoints, prosac_data) )
        {
            return(EXIT_FAILURE);
        }
        cfg.prosac.sortedPointIndices = &prosac_data[0];
    }
    else
        cfg.prosac.sortedPointIndices = NULL;

//    cfg.common.verifMethod = USACConfig::VERIF_STANDARD;

    std::vector<double> errUSAC(AMOUNT_OF_TESTS); double acumErrUSAC; int valid_tests_usac = 0;
    std::vector<double> errWGMED(AMOUNT_OF_TESTS); double acumErrWGMED; int valid_tests_wgmed = 0;
    std::vector<double> errWMEAN(AMOUNT_OF_TESTS); double acumErrWMEAN; int valid_tests_wmean = 0;

    for (int test_counter = 0;test_counter < AMOUNT_OF_TESTS; test_counter++) {
        //unsigned int current_seed = (unsigned int)time(NULL);
        //current_seed = 1490380499;
//        srand(current_seed);
//        printf("Current seed: %d\n", current_seed);


        OriginalAggHomogEstimator *homog = new OriginalAggHomogEstimator(w1, h1, w2, h2, 12, 6);
        homog->initParamsUSAC(cfg);

        // read in input data points
        cfg.common.numDataPoints = qty_of_matches;

        // set up the homography estimation problem
        homog->initDataUSAC(cfg);
        // homog->initProblem(cfg, &point_data[0]);
        homog->initProblem(cfg, &matches[0]);
        if (!homog->solve())
            return (EXIT_FAILURE);

        // Evaluate the results using the ground truth if it is provided
        if (gt_pts.size() > 0) {
            double errorUSAC = HTools::computeError(gt_pts, &homog->final_model_params_[0]);
            double errorwgmed = HTools::computeError(gt_pts, &homog->aggregated_wgmed_model_params_[0]);
            double errorwmean = HTools::computeError(gt_pts, &homog->aggregated_wmean_model_params_[0]);
            printf("Total error: USAC: %.3f, WGMED: %.3f, WMEAN: %.3f\n", errorUSAC, errorwgmed, errorwmean);
            errUSAC[test_counter] = errorUSAC;
            if (!std::isnan(errorUSAC)) {
                acumErrUSAC += errorUSAC;
                valid_tests_usac++;
            }
            errWGMED[test_counter] = errorwgmed;
            if (!std::isnan(errorwgmed)) {
                acumErrWGMED += errorwgmed;
                valid_tests_wgmed++;
            }
            errWMEAN[test_counter] = errorwmean;
            if (!std::isnan(errorwmean)) {
                acumErrWMEAN += errorwmean;
                valid_tests_wmean++;
            }
        }
        // write out results
//        size_t pos = (cfg.homog.inputFilePath).find_last_of("/\\");
//        std::string working_dir = (cfg.homog.inputFilePath).substr(0, pos + 1);
//        std::ofstream outmodel((working_dir + "H.txt").c_str());
//        for (unsigned int i = 0; i < 3; ++i)
//            for (unsigned int j = 0; j < 3; ++j)
//                outmodel << homog->final_model_params_[3 * i + j] / homog->final_model_params_[3 * 2 + 2] << " ";
//        outmodel.close();
//
//        std::ofstream outmodel2((working_dir + "Hwmean.txt").c_str());
//        for (unsigned int i = 0; i < 3; ++i)
//            for (unsigned int j = 0; j < 3; ++j)
//                outmodel2 << homog->aggregated_wmean_model_params_[3 * i + j] /
//                             homog->aggregated_wmean_model_params_[3 * 2 + 2] << " ";
//        outmodel2.close();
//
//        std::ofstream outmodel3((working_dir + "Hwgmed.txt").c_str());
//        for (unsigned int i = 0; i < 3; ++i)
//            for (unsigned int j = 0; j < 3; ++j)
//                outmodel3 << homog->aggregated_wgmed_model_params_[3 * i + j] /
//                             homog->aggregated_wgmed_model_params_[3 * 2 + 2] << " ";
//        outmodel3.close();
//
//        std::ofstream outinliers((working_dir + "inliers.txt").c_str());
//        for (unsigned int i = 0; i < cfg.common.numDataPoints; ++i)
//            outinliers << homog->usac_results_.inlier_flags_[i] << std::endl;
//        outinliers.close();

        // clean up
    //    point_data.clear();
        homog->cleanupProblem();
        delete homog;
    }
    double meanErrUSAC =  acumErrUSAC / (double)  valid_tests_usac;
    double meanErrWMEAN = acumErrWMEAN/ (double) valid_tests_wmean;
    double meanErrWGMED = acumErrWGMED / (double) valid_tests_wgmed;
    printf("\n\nAvg. errors. USAC: %f. WMEAN: %f. WGMED: %f\n", meanErrUSAC , meanErrWMEAN, meanErrWGMED);
    std::sort(errUSAC.begin(), errUSAC.end());
    std::sort(errWMEAN.begin(), errWMEAN.end());
    std::sort(errWGMED.begin(), errWGMED.end());
    // ASSUME ODD AMOUNT OF TESTS FOR MEDIAN COMPUTATION
    double stddevErrUSAC = 0, stddevErrWGMED = 0, stddevErrWMEAN = 0;
    for (int i=0;i<AMOUNT_OF_TESTS;i++) {
        if (!std::isnan(errUSAC[i]))
            stddevErrUSAC += (errUSAC[i] - meanErrUSAC) * (errUSAC[i] - meanErrUSAC);
        if (!std::isnan(errWGMED[i]))
            stddevErrWGMED += (errWGMED[i] - meanErrWGMED) * (errWGMED[i] - meanErrWGMED);
        if (!std::isnan(errWMEAN[i]))
            stddevErrWMEAN += (errWMEAN[i] - meanErrWMEAN) * (errWMEAN[i] - meanErrWMEAN);
    }
    stddevErrUSAC = sqrt(stddevErrUSAC / valid_tests_usac);
    stddevErrWMEAN = sqrt(stddevErrWMEAN / valid_tests_wmean);
    stddevErrWGMED = sqrt(stddevErrWGMED / valid_tests_wgmed);


    printf("Median errors. USAC: %f. WMEAN: %f. WGMED: %f\n", errUSAC[(AMOUNT_OF_TESTS-1)/2], errWMEAN[(AMOUNT_OF_TESTS-1)/2], errWGMED[(AMOUNT_OF_TESTS-1)/2]);
    printf("Std. dev. USAC: %f. WMEAN: %f. WGMED: %f\n", stddevErrUSAC, stddevErrWMEAN, stddevErrWGMED);
    printf("Amount of failed tests: USAC: %d. WMEAN: %d. WGMED: %d\n", AMOUNT_OF_TESTS - valid_tests_usac,
           AMOUNT_OF_TESTS - valid_tests_wmean, AMOUNT_OF_TESTS - valid_tests_wgmed);

    matches.clear();
    prosac_data.clear();
    return(EXIT_SUCCESS);
}

