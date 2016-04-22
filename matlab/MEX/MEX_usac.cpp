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
      "[F,inliers] = MEX_usac(ptype, cfgfile, [matches], [sorting])\n"
      "   \n"
      "   ptype     0:fundamental 1:homography\n"
      "   cfgfile   path to the config file\n"
      "   matches   input 2D matches (4 columns). Overrides the file specified in cfgfile\n"
      "   sorting   PROSAC match sorting, a column vector with the indeces of the matches\n"
      "             ordered by decreasing quality. Only used if PROSAC option is activated\n"
      "             in cfgfile. Overrides the file specified in cfgfile\n"
      "   \n"
      "   F         3x3 result matrix\n"
      "   inliers   column vector with inliers mask\n");
}


/*----------------------------------------------------------------------------*/
/*                                    Matlab Main Function                    */
/*----------------------------------------------------------------------------*/
/** USAC Matlab wrapper 
 *  [F,inliers] = usac(ptype, cfgfile, [points(optional)], [prosac sort])
 *
 * ptype     0:fundamental 1:homography 2:essential
 * cfgfile   path to the config file
 * points    override config file points 
 * prosac sort override the sort file specified in the config
 *
 * F         a 3x3 matrix
 * inliers   a column vector 
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
  if(nrhs>4)
  {
    help(); mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","A maximum of two inputs are required.");
  }
  if(nrhs==2 || nrhs==3 || nrhs==4)
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
      help(); mexErrMsgTxt("The third input must be a double matrix.");
    }
    if (nrhs==4 && (!mxIsDouble(prhs[3])))
    {
      help(); mexErrMsgTxt("The third input must be a double matrix.");
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

	 std::string cfg_file_path ( (const char*) cfgfile );

    if (nrhs>=3)
    {
      DATA = mxGetPr(prhs[2]);
      COLS = mxGetN(prhs[2]);
      ROWS = mxGetM(prhs[2]);

      if (COLS != 4)
      {
        help(); mexErrMsgTxt("The match matrix must have 4 columns.");
      }
    }

    if (nrhs==4)
    {
      PROSACDATA = mxGetPr(prhs[3]);

      if (mxGetN(prhs[3]) != 1 || mxGetM(prhs[3]) != ROWS)
      {
        help(); mexErrMsgTxt("The sort matrix must have 1 column, and as many rows as matches");
      }
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
			//return(EXIT_FAILURE);
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
		   	//return(EXIT_FAILURE);
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
			      //return(EXIT_FAILURE);
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
			//return(EXIT_FAILURE);
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

//		// write out results
//		size_t pos = (cfg.fund.inputFilePath).find_last_of("/\\");
//		std::string working_dir = (cfg.fund.inputFilePath).substr(0, pos + 1);
//		std::ofstream outmodel((working_dir + "F.txt").c_str());
//		for (unsigned int i = 0; i < 3; ++i)
//		{
//			for (unsigned int j = 0; j < 3; ++j)
//			{
//				outmodel << fund->final_model_params_[3*i+j] << " ";
//			}
//		}
//		outmodel.close();
//		std::ofstream outinliers((working_dir + "inliers.txt").c_str());
//		for (unsigned int i = 0; i < cfg.common.numDataPoints; ++i)
//		{
//			outinliers << fund->usac_results_.inlier_flags_[i] << std::endl;
//		}
//		outinliers.close();

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
			//return(EXIT_FAILURE);
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
		   	//return(EXIT_FAILURE);
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
			   	//return(EXIT_FAILURE);
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
			//return(EXIT_FAILURE);
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

//		// write out results
//		size_t pos = (cfg.homog.inputFilePath).find_last_of("/\\");
//		std::string working_dir = (cfg.homog.inputFilePath).substr(0, pos + 1);
//		std::ofstream outmodel((working_dir + "H.txt").c_str());
//		for (unsigned int i = 0; i < 3; ++i)
//		{
//			for (unsigned int j = 0; j < 3; ++j)
//			{
//				outmodel << homog->final_model_params_[3*i+j] << " ";
//			}
//		}
//		outmodel.close();
//		std::ofstream outinliers((working_dir + "inliers.txt").c_str());
//		for (unsigned int i = 0; i < cfg.common.numDataPoints; ++i)
//		{
//			outinliers << homog->usac_results_.inlier_flags_[i] << std::endl;
//		}
//		outinliers.close();

		// clean up
		point_data.clear();
		prosac_data.clear();
		homog->cleanupProblem();
		delete homog;
	} else {
      mexErrMsgTxt("Estimation problem currently not implemented");
		//std::cout << "Estimation problem currently not implemented" << std::endl;
	}



   }

}








// ---------------------------------------------------------------------------------------------------------------
//int mainxx(int argc, char **argv)
//{
//	// check command line args
//	if (argc < 3)
//	{
//		std::cerr << "Usage: RunSingleTest <estimation problem> <config file>" << std::endl;
//		std::cerr << "\t<estimation problem>: 0 (fundamental matrix), 1 (homography)" << std::endl;
//		std::cerr << "\t<config file>: full path to configuration file" << std::endl;
//		return(EXIT_FAILURE);
//	}
//	int estimation_problem = atoi(argv[1]);
//	std::string cfg_file_path = argv[2];
//
//	// seed random number generator
//	srand((unsigned int)time(NULL));
//
//	// initialize the appropriate robust estimation problem
//	if (estimation_problem == 0)
//	{
//		// ------------------------------------------------------------------------
//		// initialize the fundamental matrix estimation problem
//		ConfigParamsFund cfg;
//		if ( !cfg.initParamsFromConfigFile(cfg_file_path) )
//		{
//			std::cerr << "Error during initialization" << std::endl;
//			return(EXIT_FAILURE);
//		}
//		FundMatrixEstimator* fund = new FundMatrixEstimator;
//		fund->initParamsUSAC(cfg);
//
//		// read in input data points
//		std::vector<double> point_data;
//		if ( !readCorrsFromFile(cfg.fund.inputFilePath, point_data, cfg.common.numDataPoints) )
//		{
//			return(EXIT_FAILURE);
//		}
//
//		// read in prosac data if required
//		std::vector<unsigned int> prosac_data;
//		if (cfg.common.randomSamplingMethod == USACConfig::SAMP_PROSAC)
//		{
//			prosac_data.resize(cfg.common.numDataPoints);
//			if ( !readPROSACDataFromFile(cfg.prosac.sortedPointsFile, cfg.common.numDataPoints, prosac_data) )
//			{
//				return(EXIT_FAILURE);
//			}
//			cfg.prosac.sortedPointIndices = &prosac_data[0];
//		} else {
//			cfg.prosac.sortedPointIndices = NULL;
//		}
//
//		// set up the fundamental matrix estimation problem
//
//		fund->initDataUSAC(cfg);
//		fund->initProblem(cfg, &point_data[0]);
//		if (!fund->solve())
//		{
//			return(EXIT_FAILURE);
//		}
//
//		// write out results
//		size_t pos = (cfg.fund.inputFilePath).find_last_of("/\\");
//		std::string working_dir = (cfg.fund.inputFilePath).substr(0, pos + 1);
//		std::ofstream outmodel((working_dir + "F.txt").c_str());
//		for (unsigned int i = 0; i < 3; ++i)
//		{
//			for (unsigned int j = 0; j < 3; ++j)
//			{
//				outmodel << fund->final_model_params_[3*i+j] << " ";
//			}
//		}
//		outmodel.close();
//		std::ofstream outinliers((working_dir + "inliers.txt").c_str());
//		for (unsigned int i = 0; i < cfg.common.numDataPoints; ++i)
//		{
//			outinliers << fund->usac_results_.inlier_flags_[i] << std::endl;
//		}
//		outinliers.close();
//
//		// clean up
//		point_data.clear();
//		prosac_data.clear();
//		fund->cleanupProblem();
//		delete fund;
//
//	} else if (estimation_problem == 1) {
//		// ------------------------------------------------------------------------
//		// initialize the homography estimation problem
//		ConfigParamsHomog cfg;
//		if ( !cfg.initParamsFromConfigFile(cfg_file_path) )
//		{
//			std::cerr << "Error during initialization" << std::endl;
//			return(EXIT_FAILURE);
//		}
//
//		HomogEstimator* homog = new HomogEstimator;
//		homog->initParamsUSAC(cfg);
//
//		// read in input data points
//		std::vector<double> point_data;
//		if ( !readCorrsFromFile(cfg.homog.inputFilePath, point_data, cfg.common.numDataPoints) )
//		{
//			return(EXIT_FAILURE);
//		}
//
//		// read in prosac data if required
//		std::vector<unsigned int> prosac_data;
//		if (cfg.common.randomSamplingMethod == USACConfig::SAMP_PROSAC)
//		{
//			prosac_data.resize(cfg.common.numDataPoints);
//			if ( !readPROSACDataFromFile(cfg.prosac.sortedPointsFile, cfg.common.numDataPoints, prosac_data) )
//			{
//				return(EXIT_FAILURE);
//			}
//			cfg.prosac.sortedPointIndices = &prosac_data[0];
//		} else {
//			cfg.prosac.sortedPointIndices = NULL;
//		}
//
//		// set up the homography estimation problem
//		homog->initDataUSAC(cfg);
//		homog->initProblem(cfg, &point_data[0]);
//		if (!homog->solve())
//		{
//         mexErrMsgTxt("Error running the algorithm");
//			//return(EXIT_FAILURE);
//		}
//
//		// write out results
//		size_t pos = (cfg.homog.inputFilePath).find_last_of("/\\");
//		std::string working_dir = (cfg.homog.inputFilePath).substr(0, pos + 1);
//		std::ofstream outmodel((working_dir + "H.txt").c_str());
//		for (unsigned int i = 0; i < 3; ++i)
//		{
//			for (unsigned int j = 0; j < 3; ++j)
//			{
//				outmodel << homog->final_model_params_[3*i+j] << " ";
//			}
//		}
//		outmodel.close();
//		std::ofstream outinliers((working_dir + "inliers.txt").c_str());
//		for (unsigned int i = 0; i < cfg.common.numDataPoints; ++i)
//		{
//			outinliers << homog->usac_results_.inlier_flags_[i] << std::endl;
//		}
//		outinliers.close();
//
//		// clean up
//		point_data.clear();
//		prosac_data.clear();
//		homog->cleanupProblem();
//		delete homog;
//	} else {
//		std::cout << "Estimation problem currently not implemented" << std::endl;
//	}
//
//	return(EXIT_SUCCESS);
//}
//
