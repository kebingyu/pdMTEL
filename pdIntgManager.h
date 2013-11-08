/*

The pdIntgManager class takes care of the integration over volume calculation.

*/

#ifndef PDINTGMANAGER_H
#define PDINTGMANAGER_H

#include <limits>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "pdBond.h"
#include "pdNode.h"
#include "pdVector.h"
#include "pdMatrix.h"
#include "pdIntgPoint.h"
#include "pdBdryConditionDispl.h"
#include "pdBdryConditionDisplGrad.h"
using namespace std;

const double PI = acos(-1.0);
const double MEGA = numeric_limits<double>::max();
const double TINY = numeric_limits<double>::min();

// here JMAX limits the total number of steps
const double TOL = 1.0e-6; // tolerance for error control
const int JMAX = 20;
const int JMAXP = JMAX + 1;
const int SMAX = 10;

// MLS 
const int POLYBASISDIM = 10; // polynomial basis vector dimension (order 2)
const double ILLCONDITION = 1.0e10; // if the condition number of the moment matrix is higher than this, it is considered ill-conditioned


class pdIntgManager;
// function pointers
//typedef double (pdIntgManager::*FUNC_T)(double, double, double);
//typedef double (pdIntgManager::*FUNC_D)(double, double);
//typedef double (pdIntgManager::*FUNC_DV)(double, vector<string>);
//typedef double (pdIntgManager::*FUNC_S)(double);

typedef void (pdIntgManager::*FUNC_T)(double, double, double, double*);
typedef double (pdIntgManager::*FUNC_D)(double, double);
typedef void (pdIntgManager::*FUNC_DV)(double, vector<string>, double*);
typedef double (pdIntgManager::*FUNC_S)(double);

class pdIntgManager
{
public:
	// Constructor	
	pdIntgManager(const double, const double, pdBond*, pdNode*, pdNode*);
	pdIntgManager(const double, const double, pdBond*, const double*, pdNode*, pdNode*);
	
	// Destructor
	~pdIntgManager();

	// Function	
	void CalcIntg3D(string, double, double, vector<string>, double*);	
	void AdaptiveIntegration(string, double*);
	void CalcConfigurationType3(string, double*);
	void CalcBondBasedForce(const double, const double, const double, double*);
	//void CalcVolumeToPointForceDensity(const double, const double, const double, double*);
	int GetEndIntgPoint(int);

	void SetTrapFlag(const bool);
	void SetErrorControlFlag(const bool);
	void SetErrorControlEPS(const double);
	void SetGaussianFlag(const bool);
	void SetGaussianAbscis(const vector<double>&);
	void SetGaussianWeight(const vector<double>&);
	void SetMLSFlag(const bool);
	void SetIntgPoint(const int);
	void SetInitMLSCoeffFlag(const bool);
	void SetCBNodeDispl(double**);
	void SetCurrentRunTime(const double);

private:
	double qgaus(FUNC_S, double, double);
	void sobseq(int *n, double x[]);

	// Copy Control
	pdIntgManager(const pdIntgManager&);
	pdIntgManager& operator=(const pdIntgManager&);
	
	// Function	
	void CalcTrapzdQuad(FUNC_DV, double, double, vector<string>, int, int, double*);
	void CalcGaussianQuad(FUNC_DV, double, double, vector<string>, int, int, double*);
	void CalcIntgX1(double, vector<string>, double*);
	void CalcIntgX2(double, vector<string>, double*);
	void CalcIntgX3(double, vector<string>, double*);
	double SelectYLimitFuncValue(double, string);
	double SelectZLimitFuncValue(double, double, string);	
	void CalcConfigurationType1(string, double*);
	void CalcConfigurationType2(string, double*);
	int FindIntersectionPointOnEightEdges(vector<string>, double*);
	int FindIntersectionPointOnFourEdges(vector<string>, double*);
	void CalcConfigurationType1Sub1(string, string, string, string, double, double*);
	void CalcConfigurationType1Sub2(string, string, string, string, double, double*);
	void CalcConfigurationType1Sub3(string, string, string, string, double*);
	void CalcConfigurationType1Sub4(string, string, string, string, double, double*);
	void CalcConfigurationType1Sub5(string, string, string, string, double, double*);
	void CalcConfigurationType2Sub1(string, string, string, string, double*);
	void CalcConfigurationType2Sub2(string, string, string, string, double*);
	double GetCoordAtDirection(string);
	double SolveBinaryQuadric(double, double, double, int, double);
	void GetAppDispl(double, double, double, vector<pdNode*>, double*);
	pdVector* GetBasis(double, double, double, int);
	double GetWeight(pdNode*, double, double, double, double);
	bool ConsistencyCheck(double, double, double, vector<pdNode*>, pdVector*);
	double CalcBetaFactor(double dist, double delta, double dx);
	void ArrayAddition(double* source, double* target, const int length);
	bool IsConverge(const double, const double, const int, int&);
	void InitMLSCoeff(const double, const double, const double, pdNode*);

	// Table of integrand function
	void UnitValue(double, double, double, double*);
	void GetStrainEnergyDensity(double, double, double, double*);
	

	// Table of integration limit function for Y
	double YLimitXLow(double);
	double YLimitXHigh(double);
	double YLimitYLow(double);
	double YLimitYHigh(double);
	double YLimitZLow(double);
	double YLimitZHigh(double);
	double YLimitHighXYPos(double);
	double YLimitHighXYNeg(double);
	double YLimitHighXZPos(double);
	double YLimitHighXZNeg(double);
	double YLimitHighYXPos(double);
	double YLimitHighYXNeg(double);
	double YLimitHighYZPos(double);
	double YLimitHighYZNeg(double);
	double YLimitHighZXPos(double);
	double YLimitHighZXNeg(double);
	double YLimitHighZYPos(double);
	double YLimitHighZYNeg(double);
	double YLimitLowXYPos(double);
	double YLimitLowXYNeg(double);
	double YLimitLowXZPos(double);
	double YLimitLowXZNeg(double);
	double YLimitLowYXPos(double);
	double YLimitLowYXNeg(double);
	double YLimitLowYZPos(double);
	double YLimitLowYZNeg(double);
	double YLimitLowZXPos(double);
	double YLimitLowZXNeg(double);
	double YLimitLowZYPos(double);
	double YLimitLowZYNeg(double);	

	// Table of integration limit function for Z
	double ZLimitZSpherePos(double, double);
	double ZLimitZSphereNeg(double, double);
	double ZLimitZLow(double, double);
	double ZLimitZHigh(double, double);
	double ZLimitYSpherePos(double, double);
	double ZLimitYSphereNeg(double, double);
	double ZLimitYLow(double, double);
	double ZLimitYHigh(double, double);
	double ZLimitXSpherePos(double, double);
	double ZLimitXSphereNeg(double, double);
	double ZLimitXLow(double, double);
	double ZLimitXHigh(double, double);


	// Function pointer selector
	FUNC_T SelectIntegrand(string);
	FUNC_S SelectYLimitFunc(string);
	FUNC_D SelectZLimitFunc(string);
	FUNC_S SelectSingleArguFunc(string);

	// Data
	double _xSave, _ySave;	
	double _delta; // horizon
	double _dx; // grid spacing	
	double _x1i, _x2i, _x3i; // lcoal coordinate of node i which is (0, 0, 0) since it is located at the center
	double _u1i, _u2i, _u3i;
	double _x1j, _x2j, _x3j; // local coordinate of node j (shift from its global coordinate)
	double _eps; // desired error control accuracy (used in the adaptive integration)
	int _intgPoint; 
	/*
	  1. In adaptive integration, _intgPoint means the starting trapezoidal index.
	  2. In fixed Gaussian integration, it means the number of Gaussian points in each direction.
	*/
	int _intgPointEnd[4];
	/* 
	  1. _intgPointEnd[1]~[3] are used in trapezoidal quadrature to 
		store the number of the endding trapezoidal (integration) points in three directions.
	  2. _intgPointEnd[0] is used in Gaussian quadrature to store the total number of Gaussian (integration) points.
	*/							  
	int _intgPointCounter; // integration point counter
	/*
	  1. This index only works for the Gaussian integration since only one integration equation is used for all 
	     possible configurations.
      2. Two reasons to make it malfunctioning in trapezoidal integration
	     1) In some configurations, the integration is done by two or three parts. Therefore the index is resetted 
		    once one part is done.
		 2) In some configurations, the lower and upper integration limits can be the same. Therefore the function 
		    call is returned in CalcIntgX1 or CalcIntgX2. CalcIntgX3 is never called and the index is not self-added. 
	*/
	bool _aiOn; // trapezoidal quadrature method flag
	bool _errorControlOn; // error control flag for adaptive integration
	bool _fgiOn; // Gaussian quadrature method flag
	bool _mlsOn; // moving least square approximation method flag
	bool _initMLSCoeffOn;
	vector<double> _abscis; // store abscissas of Gaussian points
	vector<double> _weights; // store weights of Gaussian points
	pdBond* _bond;
	pdNode* _node_i; // the source node
	pdNode* _node_j; // the family node
	FUNC_T _nrfunc; // function pointer to the integrand function
	// Adaptive integration with MLS (error control (adaptivity) is manually controlled in this mode to save computation time)
	double** _cbNodeDispl; // contain trial displacements of contributing nodes of node_j from RK4 steps
	double _curRunTime; // the current run time. It is used to determine the boundary condition value
};

#endif