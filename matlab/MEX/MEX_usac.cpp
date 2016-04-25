//
// Created by Gabriele Facciolo on 22/04/16.
//
// based on the porting of USAC by Christian Richardt (http://richardt.name).
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

#include "config/ConfigParams.h"
#include "estimators/FundmatrixEstimator.h"
#include "estimators/HomogEstimator.h"

#include "mex.h"

// helper functions
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


int readCorrsFromPointer(double *DATA, int ROWS, std::vector<double> &point_data) {
         for(int i=0;i<ROWS;i++) {
            point_data.push_back(DATA[i+0*ROWS]);
            point_data.push_back(DATA[i+1*ROWS]);
            point_data.push_back(1.0);
            point_data.push_back(DATA[i+2*ROWS]);
            point_data.push_back(DATA[i+3*ROWS]);
            point_data.push_back(1.0);
         }
         return ROWS;
}

int readPROSACDataFromPointer(double *DATA, int ROWS, std::vector<unsigned int> &data) {
         for(int i=0;i<ROWS;i++) {
            data.push_back(DATA[i]);
         }
         return ROWS;
}


void help() {
      mexErrMsgTxt("USAC Matlab wrapper\n\n"
      "[F,inliers] = MEX_usac(ptype, cfgfile, [matches], [sorting], [inlier_thres], [max_hypoth])\n"
      "   \n"
      "   ptype        0:fundamental 1:homography\n"
      "   cfgfile      path to the config file\n"
      "   matches      input 2D matches (4 columns). Overrides the file specified in cfgfile\n"
      "   sorting      PROSAC match sorting, a column vector with the indeces of the matches\n"
      "                ordered by decreasing quality. Use [] to deactivate PROSAC.\n"
      "                Overrides file specified in cfgfile\n"
      "   inlier_thres inlier_threshold parameter. Overrides value specified in cfgfile\n"
      "   max_hypoth   max_hypotheses parameter. Overrides value specified in cfgfile\n"
      "   \n"
      "   F         3x3 result matrix\n"
      "   inliers   column vector with inliers mask\n");
}


/*----------------------------------------------------------------------------*/
/*                                    Matlab Main Function                    */
/*----------------------------------------------------------------------------*/
/** USAC Matlab wrapper 
 * [F,inliers] = MEX_usac(ptype, cfgfile, [verbose], [matches], [sorting], [inlier_thres], [max_hypoth])
 *      ptype        0:fundamental 1:homography 2:essential NOT IMPLEMENTED
 *      cfgfile      path to the config file
 *      verbose      0: silent mode. 1: verbose mode. Overrides the config in cfgfile
 *      matches      input 2D matches (4 columns). Overrides the file specified in cfgfile
 *      sorting      PROSAC match sorting, a column vector with the indeces of the matches
 *                   ordered by decreasing quality. Use [] to deactivate PROSAC.
 *                   Overrides file specified in cfgfile
 *      inlier_thres inlier_threshold parameter. Overrides value specified in cfgfile
 *      max_hypoth   max_hypotheses parameter. Overrides value specified in cfgfile
 *      
 *      F         3x3 result matrix
 *      inliers   column vector with inliers mask
 * */


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (!(nlhs==2))
  {
    help(); mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Two outputs are required.");
  }
  if(nrhs<2)
  {
    help(); mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","At least two inputs are required.");
  }
  if(nrhs>7)
  {
    help(); mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","A maximum of seven inputs are required.");
  }
  if(nrhs==2 || nrhs==3 || nrhs==4 || nrhs==5 || nrhs==6 || nrhs==7)
  {

    if (!(mxIsDouble(prhs[0])))
    {
      help(); mexErrMsgTxt("The first input must be a double.");
    }
    if (!(mxIsChar(prhs[1])))
    {
      help(); mexErrMsgTxt("The second input must be a char array.");
    }
    if (nrhs==3 && (!mxIsDouble(prhs[2])))
    {
      help(); mexErrMsgTxt("The third input must be a double.");
    }
    if (nrhs==4 && (!mxIsDouble(prhs[3])))
    {
      help(); mexErrMsgTxt("The fourth input must be a double matrix.");
    }
    if (nrhs==5 && (!mxIsDouble(prhs[4])))
    {
      help(); mexErrMsgTxt("The fifth input must be a double matrix.");
    }
    if (nrhs==6 && (!mxIsDouble(prhs[5])))
    {
      help(); mexErrMsgTxt("The sixth input must be a double.");
    }
    if (nrhs==7 && (!mxIsDouble(prhs[6])))
    {
      help(); mexErrMsgTxt("The seventh input must be a double.");
    }


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//	int estimation_problem = atoi(argv[1]);
//	std::string cfg_file_path = argv[2];

    double *DATA = NULL, *PROSACDATA = NULL;
    int COLS=-1, ROWS=-1;
    int estimation_problem = mxGetPr(prhs[0])[0];
    char *cfgfile = mxArrayToString(prhs[1]);
    double inlier_threshold=-1;
    int max_hypotheses=-1;
    int activate_prosac = 0;   // -1: deactivate (explicitly []);   0: see cfgfile;   1: activate (explcitly PROSACDATA);
    int verbose = -1; //     -1: see cfgfile (explicitly []);   0: false; 1: true;

	 std::string cfg_file_path ( (const char*) cfgfile );

    if (nrhs >= 3)
    {
      if (!mxIsEmpty(prhs[2]) && mxIsLogical(prhs[2]))
      {
        bool bVerbose = mxGetLogicals(prhs[2])[0];
        if (bVerbose)
            verbose = 1;
        else
            verbose = 0;
      }
    }

    if (nrhs>=4)
    {
      DATA = mxGetPr(prhs[3]);
      COLS = mxGetN(prhs[3]);
      ROWS = mxGetM(prhs[3]);

      if (COLS != 4)
      {
        help(); mexErrMsgTxt("The match matrix must have 4 columns.");
      }
    }

    if (nrhs>=5)
    {
      PROSACDATA = mxGetPr(prhs[4]);
      activate_prosac = 1;

      if (mxGetN(prhs[4]) != 1 || mxGetM(prhs[4]) != ROWS)
      {
         PROSACDATA = NULL;
         activate_prosac = -1;
         if (verbose == 1)
            mexPrintf("Deactivating PROSAC (sorting matrix must have 1 column, and as many rows as matches)\n");
         
      }
    }

    if (nrhs>=6)
    {
      inlier_threshold = mxGetPr(prhs[5])[0];
      if (verbose == 1)
         mexPrintf("Changing inlier_threshold = %f\n", inlier_threshold);
    }
    if (nrhs>=7)
    {
      max_hypotheses = mxGetPr(prhs[6])[0];
      if (verbose == 1)
         mexPrintf("Changing max_hypotheses = %d\n", max_hypotheses);
    }




	// seed random number generator
	srand((unsigned int)time(NULL));

	// initialize the appropriate robust estimation problem
	if (estimation_problem == 0)
	{
		// ------------------------------------------------------------------------
		// initialize the fundamental matrix estimation problem
		ConfigParamsFund cfg;
		if ( !cfg.initParamsFromConfigFile(cfg_file_path) )
		{
         mexErrMsgTxt("Error during initialization");
		}
      // Override PROSAC activate or deactivate
		if (activate_prosac > 0) {
         cfg.common.randomSamplingMethod = USACConfig::SAMP_PROSAC;
      } else if (activate_prosac < 0) {
         cfg.common.randomSamplingMethod = USACConfig::SAMP_UNIFORM;
      }
      // Override inlier_threshold and max_hypotheses
      if (inlier_threshold >= 0) {
         cfg.common.inlierThreshold = inlier_threshold;
      }
      if (max_hypotheses >= 0) {
         cfg.common.maxHypotheses = max_hypotheses;
      }
      if (verbose >= 0) {
         cfg.common.verbose = (verbose == 1);
      }

		FundMatrixEstimator* fund = new FundMatrixEstimator;
		fund->initParamsUSAC(cfg);

		// read in input data points
		std::vector<double> point_data;
      if (DATA) {
         cfg.common.numDataPoints = readCorrsFromPointer(DATA, ROWS, point_data);
      } else {
		   if ( !readCorrsFromFile(cfg.fund.inputFilePath, point_data, cfg.common.numDataPoints) )
		   {
            mexErrMsgTxt("Error reading input data points");
		   }
      }

		// read in prosac data if required
		std::vector<unsigned int> prosac_data;
		if (cfg.common.randomSamplingMethod == USACConfig::SAMP_PROSAC)
		{
         if (PROSACDATA) {
            readPROSACDataFromPointer(PROSACDATA, ROWS, prosac_data);
			   cfg.prosac.sortedPointIndices = &prosac_data[0];
         } else {
			   prosac_data.resize(cfg.common.numDataPoints);
			   if ( !readPROSACDataFromFile(cfg.prosac.sortedPointsFile, cfg.common.numDataPoints, prosac_data) )
			   {
               mexErrMsgTxt("Error reading prosac data");
			   }
			   cfg.prosac.sortedPointIndices = &prosac_data[0];
         }
		} else {
			cfg.prosac.sortedPointIndices = NULL;
		}

		// set up the fundamental matrix estimation problem
		fund->initDataUSAC(cfg);
		fund->initProblem(cfg, &point_data[0]);
		if (!fund->solve())
		{
         mexErrMsgTxt("Error running the algorithm");
		}

		// write out results
      /* The output which is created here will be an array */
      plhs[0]=mxCreateDoubleMatrix(3,3,mxREAL);
      double* Mat=(double*)mxGetPr(plhs[0]);
		for (unsigned int i = 0; i < 3; ++i)
			for (unsigned int j = 0; j < 3; ++j)
				Mat[i+j*3] = fund->final_model_params_[3*i+j];


      plhs[1]=mxCreateDoubleMatrix(cfg.common.numDataPoints,1,mxREAL);
      double* inliers=(double*)mxGetPr(plhs[1]);
		for (unsigned int i = 0; i < cfg.common.numDataPoints; ++i)
			inliers[i] = fund->usac_results_.inlier_flags_[i];

		// clean up
		point_data.clear();
		prosac_data.clear();
		fund->cleanupProblem();
		delete fund;

	} else if (estimation_problem == 1) {
		// ------------------------------------------------------------------------
		// initialize the homography estimation problem
		ConfigParamsHomog cfg;
		if ( !cfg.initParamsFromConfigFile(cfg_file_path) )
		{
         mexErrMsgTxt("Error during initialization");
		}
      // Override PROSAC activate or deactivate
		if (activate_prosac > 0) {
         cfg.common.randomSamplingMethod = USACConfig::SAMP_PROSAC;
      } else if (activate_prosac < 0) {
         cfg.common.randomSamplingMethod = USACConfig::SAMP_UNIFORM;
      }
      // Override inlier_threshold and max_hypotheses
      if (inlier_threshold >= 0) {
         cfg.common.inlierThreshold = inlier_threshold;
      }
      if (max_hypotheses >= 0) {
         cfg.common.maxHypotheses = max_hypotheses;
      }
      if (verbose >= 0) {
         cfg.common.verbose = (verbose == 1);
      }

		HomogEstimator* homog = new HomogEstimator;
		homog->initParamsUSAC(cfg);

		// read in input data points
		std::vector<double> point_data;
      if (DATA) {
         cfg.common.numDataPoints = readCorrsFromPointer(DATA, ROWS, point_data);
      } else {
		   if ( !readCorrsFromFile(cfg.homog.inputFilePath, point_data, cfg.common.numDataPoints) )
		   {
            mexErrMsgTxt("Error reading data points");
		   }
      }

		// read in prosac data if required
		std::vector<unsigned int> prosac_data;
		if (cfg.common.randomSamplingMethod == USACConfig::SAMP_PROSAC)
		{
         if (PROSACDATA) {
            readPROSACDataFromPointer(PROSACDATA, ROWS, prosac_data);
			   cfg.prosac.sortedPointIndices = &prosac_data[0];
         } else {
			   prosac_data.resize(cfg.common.numDataPoints);
			   if ( !readPROSACDataFromFile(cfg.prosac.sortedPointsFile, cfg.common.numDataPoints, prosac_data) )
			   {
               mexErrMsgTxt("Error reading prosac data");
			   }
			   cfg.prosac.sortedPointIndices = &prosac_data[0];
         }
		} else {
			cfg.prosac.sortedPointIndices = NULL;
		}

		// set up the homography estimation problem
		homog->initDataUSAC(cfg);
		homog->initProblem(cfg, &point_data[0]);
		if (!homog->solve())
		{
         mexErrMsgTxt("Error running the algorithm");
		}

		// write out results
      /* The output which is created here will be an array */
      plhs[0]=mxCreateDoubleMatrix(3,3,mxREAL);
      double* Mat=(double*)mxGetPr(plhs[0]);
		for (unsigned int i = 0; i < 3; ++i)
			for (unsigned int j = 0; j < 3; ++j)
				Mat[i+j*3] = homog->final_model_params_[3*i+j];


      plhs[1]=mxCreateDoubleMatrix(ROWS,1,mxREAL);
      double* inliers=(double*)mxGetPr(plhs[1]);
		for (unsigned int i = 0; i < cfg.common.numDataPoints; ++i)
			inliers[i] = homog->usac_results_.inlier_flags_[i];

		// clean up
		point_data.clear();
		prosac_data.clear();
		homog->cleanupProblem();
		delete homog;
	} else {
      mexErrMsgTxt("Estimation problem currently not implemented");
	}



   }

}
