/*

The pdDataManager class contains all the data for problem setup.

*/

#ifndef PDDATAMANAGER_H
#define PDDATAMANAGER_H

#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include "math.h"
#include "pdNode.h"
#include "pdBond.h"
#include "pdMaterialBond.h"
#include "pdMaterialState.h"
#include "pdSpace.h"
#include "pdSpaceMaterial.h"
#include "pdSpaceBdryCondition.h"
#include "pdSpaceDeletion.h"
#include "pdSpacePlane.h"
#include "pdBdryConditionDispl.h"
#include "pdBdryConditionVeloc.h"
#include "pdBdryConditionDisplGrad.h"
#include "pdBdryConditionTract.h"
#include "pdHandle.cpp"
#include "pdIntgManager.h"
#include "pdVector.h"
#include "pdMatrix.h"
//#include "mpi.h"
using namespace std;
using std::string;

const int BASISDIM = 10; // MLS polynomial basis vector dimension (order of 2)

class pdDataManager
{
public:
	// Constructor
	pdDataManager();

	// Destructor
	~pdDataManager();

	// Major functions (in pdDataManager.cpp)
	void ReadInputFile();
	void InitNode();
	void InitGridBoundary();
	void InitMaterial();
	void InitBdryEffectFactor();
	void CalcStateBasedMotion(const double);
	void CalcBondBasedMotion(const double, const double);
	void CalcBondBasedMotionRK4(const double, const double, const vector<double>&, vector<double>&);
	void RK4(const double, const double);
	void UpdateNodeStatus(const double, const double);
	void CalcIntgPointDisplRK4(const double, const vector<double>&);
	void UpdateNodeStatusWithBC(const double);
	void UpdateIntgPointStatusWithBC(const double);
	void CalcApproxDispl(pdNode*);
	void InitMPI();
	void InitLocalBonds();

	// Basic functions (in pdDataManager2.cpp)	
	void InitBond();
	void InitIntegrationPoint();
	void SetGaussianPoint();
	void InitProcBdry();
	void InitDynamicSolve();
	double CalcStableTimeStep();
	void CalcStateBasedForce(pdNode*, const double);
	void CalcBondForceIntg(pdBond*, pdNode*, pdNode*, double*, int&);
	void CalcBondForceIntg(pdBond*, const double*, pdNode*, double**, pdNode*, const double, double*, int&);
	void CalcBondForceCCI(pdBond*, pdNode*, pdNode*, double*, int&);
	void CalcBondForceCCI(const double*, const double*, const double, const double, const double, const double, double*, double&);
	void CalcArealForceDensity(pdNode*, double*, int&, ofstream&);
	void CalcBdryEffFactor(pdNode*);
	double CalcBdryEffFactorTypeOne(double, pdMaterial*);
	//void CalcForceThroughPlane(pdBond*, const double [], const double);
	template <typename T>
	void ReadKeyword(const string, const int, T&);
	bool FindKeyword(const string, const int);
	pdVector* GetBasis(double, double, double, int);
	double GetWeight(pdNode*, double, double, double, double);
	bool ConsistencyCheck(double, double, double, vector<pdNode*>, pdVector*);
	double CalcDampingCoeff(const double, pdVector*);
	void CalcVolumeToPointForceDensity(const double*, const double, double*);
	void CalcVolumeToPointForceDensity(const pdNode*, double*);

	// other functions (in pdDataManager3.cpp)
	void ReadExternalGridFile(ifstream&);
	void ReadEMUPlotFile(ifstream&, ofstream&) const;
	void PostFormatData(ifstream&, ofstream&) const;
	void PostAFD(ifstream&, ofstream&);
	void CalcDistanceToCube(pdNode*, pdNode*, double&, double&);
	void CalcDistanceToFace(int, double, pdNode*, pdNode*, double&, double&);
	double GetShortestBetween(double x_in, double x_lo, double x_hi) const;
	double GetLongestBetween(double x_in, double x_lo, double x_hi) const;
	void CalcForceOfPlusSpace(pdNode*, pdNode*, const int, double*, int&, ofstream&);
	void TestbdryEffFactor(pdNode*, ofstream&);
	void SetGaussianPointsValue();
	string KeywordAppend(const string, const int) const;
	int CalcKroneckerDelta(const int, const int);
	double CalcPermutationSymbol(const int, const int, const int);
	string GetProblemTitle() const;
	int GetReadKeywordError() const;
	int GetNumNode() const;
	int GetNumBond() const;
	int GetNumMaterial() const;
	int GetNumBC() const;
	int GetNumMatSpace() const;
	int GetNumDelSpace() const;
	int GetNumBdSpace() const;
	int GetStepDumpFrequency() const;
	int GetGridDimensionX1() const;
	int GetGridDimensionX2() const;
	int GetGridDimensionX3() const;
	int GetTotalStep() const;
	int GetTrapIndex() const;
	double GetHorizon() const;
	double GetGridCenterX1() const;
	double GetGridCenterX2() const;
	double GetGridCenterX3() const;
	double GetGridSpacing() const;
	double GetBoundingX1Low() const;
	double GetBoundingX1High() const;
	double GetBoundingX2Low() const;
	double GetBoundingX2High() const;
	double GetBoundingX3Low() const;
	double GetBoundingX3High() const;
	double GetTotalTime() const;
	double GetElastTotal() const;
	double GetKinetTotal() const;
	double GetWorkBdryTotal() const;
	double GetElastStep() const;
	double GetWorkBdryStep() const;
	double GetX1force() const;
	double GetX2force() const;
	double GetX3force() const;
	double GetForce1() const;
	double GetForce2() const;
	double GetForce3() const;
	double GetErroControlValue() const;
	bool IsStateOn() const;
	bool IsExtGridFileOn() const;
	bool IsFixedTimeStepOn() const;
	bool IsErrorControlOn() const;
	bool IsMLSOn() const;
	pdNode* GetNode(const int) const;
	pdBond* GetBond(const int) const;
	pdMaterial* GetMaterial(const int) const;
	pdBdryCondition* GetBC(const int) const;
	pdSpaceMaterial* GetMaterialSpace(const int) const;
	pdSpaceDeletion* GetDelSpace(const int) const;
	pdSpaceBdryCondition* GetBdrySpace(const int) const;
	pdSpacePlane* GetPlane(const int) const;
	void Debug();
	void PrintForceThroughPlane(ofstream&, const int) const;
	void PrintData(ofstream&) const;
	void SetWkDir();
	string GetWkDir() const;
	int GetRank() const;
	double CalcBetaFactor(double dist, double delta, double dx);

private:
	// Data
	vector<pdNode*> _nodes;	// global data containers
    vector<pdBond*> _bonds;
    vector<pdSpacePlane*> _planes;
    vector<pdSpaceDeletion*> _delSpaces;
    vector<pdSpaceMaterial*> _matSpaces;
    vector<pdSpaceBdryCondition*> _bdrySpaces;
	vector<vector<double>> _tractSpaces; // store the space limis and strain field of each traction boundary condition (07/19/2011)
    vector<pdHandle<pdMaterial>> _mats;
    vector<pdHandle<pdBdryCondition>> _BCs;
	vector<double> _abscis; // Gaussion points data
	vector<double> _weights;    
    string _title;	// title of the problem
    string _workDir; // the working folder where input and output files are located
    string _inFilePath, _outFilePath;
    int _rderr;			// flag for error type when reading keywd from input file:
    int _numNodes;	// total number of nodes
	int _numMaxFamily; // total number of family nodes in a full horizon sphere
    int _numBonds;	// total number of bonds
    int _numMats;		// total number of maertials
    int _numBCs;		// total number of boundary conditions
    int _numMatsSpaces;	// number of material Spaces
    int _numDelSpaces;	// number of deletion Spaces
    int _numBdSpaces;	// number of boundary Spaces
	int _gridDimensionX1, _gridDimensionX2, _gridDimensionX3;	// grid dimension
	int _intgPoint; // starting trapezoidal index for adaptive integration
	int _totalStep;	// total time step
    int _stepDumpFrequency;	// frequency to write nodal data
    double _horizon;	// horizon for entire grid
    double _gridSpacing;	// global grid spacing
    double _gridCenterX1, _gridCenterX2, _gridCenterX3;	// coord of grid center	
    double _gridBoundaryX1Low, _gridBoundaryX1High;	// grid boundary limits
    double _gridBoundaryX2Low, _gridBoundaryX2High;
    double _gridBoundaryX3Low, _gridBoundaryX3High;
    double _minMaterialDensity;	// minimum material density
    double _maxBulkSoundSpeed;	// maximum material bulk sound speed
    double _minBondLength;	// minimum bond length
    double _forcePlaneX1, _forcePlaneX2, _forcePlaneX3;	// total force plane position
    double _forceX1, _forceX2, _forceX3;	// force on total force plane
    double _elasticEnergyTotal;	// total elastic energy of all nodes
    double _elasticEnergyStep;	// total elastic energy of all nodes at each time step 
    double _kineticEnergyStep;	// total kinetic energy of all nodes at each time step
    double _externalWorkTotal;	// total work done through the boundary
    double _externalWorkStep;
	double _totalTime;	// total running time
    double _fixedTimeStep;	// fixed time step to replace the default stable time step
    double _safeFactor; // safety factor applied to stable time step
    double _tolerance;	// default tolerance
	double _dampCoeff; // linear damping coefficient
	double _eps; // desired accuracy for adaptive integration
	double _rw; // support domain size, used to get the approxmated displacement at arbitrary integration points
	bool _stateMatOn;			// flag: true if state-based material is used 
    bool _extGridFileOn;			// flag: true if external grid file is used
    bool _fixedTimeStepOn;	// flag: true if fixed time step is used
	bool _bdryEffOn; // flag: true if boundary effect compensation factor is turned on	
	bool _cciOn; // flag: true if cubic-cell integration method is used
	bool _cciMdf; // flag: 1 for modified CCI beta (same as 1*1*1 Gaussian integration) 0 for original CCI beta
	bool _fgiOn; // flag: true if integration with fixed Gaussian points method is used
	bool _mlsOn; // flag: true if moving least square approximation is used
	bool _aiOn; // flag: true if adaptive integration method is used
	bool _errorControlOn; // flag: true if error control is used (adaptive integration)		
    enum SpaceGeometryType { CUBOID = 1, CYLINDER };
    enum MaterialType { BOND = 1, STATE };
    enum BCType { DISPL = 1, VELOC, DISPLGRAD, TRACTION };
	enum KeywordType { OPT = 0, REQ }; // OPT=optional, REQ=required
	vector<int> _rowOfNode; // store the row number for 1st dof of each node
	bool _rk4On; // flaG: true if explicit fourth-order Runge-Kutta method is used.

	// MPI data
	int _proc1, _proc2, _proc3; 	// number of processores allocated in three direcions
	int _numProcs; // total number of processors
	int _myRank; // rank of current processor
	vector<pdBond*> _myBonds;
	double _myBdryX1Low, _myBdryX1High; // boundary limits owned by this processor
	double _myBdryX2Low, _myBdryX2High;
	double _myBdryX3Low, _myBdryX3High; 
	double _myNeedBdryX1Low, _myNeedBdryX1High; // boundary limits needed by this processor (add buffer area)
	double _myNeedBdryX2Low, _myNeedBdryX2High;
	double _myNeedBdryX3Low, _myNeedBdryX3High;

	// Function
	template <typename T>
	void ClearVectorList(vector<T*>& vec);

	// copy control
	pdDataManager(const pdDataManager&);
	pdDataManager& operator=(const pdDataManager&);

	//// state-based parameter ??
	//int num_nodes_in_state, nbond_not_broken;
	//double ident[3][3], def_grad[3][3], def_gra_inv[3][3], shape[3][3], shape_inv[3][3];
	//double sum_disp_grad[3][3], disp_grad[3][3], vel_grad[3][3], sum_vel_grad[3][3];
	//double damage_state;

	//// temporary store ?? (better removed)
	//int node_state, node_type_state, node_matmdl_state;
	//double refpos_state[3], disp_state[3], vel_state[3], deform_state[3], force_state[3], vel_family[3]; 
	//double stretch_state, exten_state, vol_state;
	//bool broken_state_old, broken_state_new;
	//double t_shape[3][3], t_disp_grad[3][3], t_vel_grad[3][3];
};

template <typename T>
inline void pdDataManager::ReadKeyword(const string keywd, const int position, T& value)
{
	/* 
	   This function reads value under the keyword at desired position from the input file

	   Input:
		keywd: desired keyword
		position: the value index of the desired keyword (positive integer)

	   Output:
		value: the value at desired position

	   Read error type:
		0 for position is not a positive integer
	    1 for fail to find keyword in the input file
	    2 for missing value at desired position
	    3 for keywd length error
	*/

	ifstream infile(_inFilePath.c_str(), ios::in);
	ofstream outfile(_outFilePath.c_str(), ios::app);
	try
	{
		// check if keyword length is larger than 1
		if (keywd.empty())
		{
			_rderr = 3;
			throw "ReadKeyword: keyword empty!\n";
		}
		// check if the position is not a positive integer
		if (position <= 0)
		{
			_rderr = 0;
			throw "ReadKeyword: position is not a positive integer!, keyword = " + keywd + "\n";
		}

		// start the search
		int find = 0;
		string line;
		do 
		{
			getline(infile, line);
			// erase the space at two ends of the read-in string
			line.erase(line.find_last_not_of(' ')+1, string::npos);
			line.erase(0, line.find_first_not_of(' '));
			if (line==keywd)  // keyword found
			{
				++find;
				// read the line next to the keyword line
				getline(infile, line); 
				// erase the space at two ends of the string
				line.erase(line.find_last_not_of(' ')+1, string::npos);
				line.erase(0, line.find_first_not_of(' '));
				// convert string as input stream
				istringstream iss(line, istringstream::in);
				// output value at defined position
				for (int i=1;!iss.eof();++i) 
				{
					iss >> value;
					if (i<position && iss.eof()) // missing value at desired position
					{
						_rderr = 2;
						throw "ReadKeyword: value missing at desired position, keyword = " + keywd + "\n";
					}
					if (i==position) // value found at desired position
					{
						infile.close();
						outfile.close();
						return;
					}
				}	
			}
		}
		while (!infile.eof());

		// keyword not in the input file
		if (find==0)
		{
			_rderr = 1;
		}
		infile.close();
		outfile.close();
		return;
	}
	catch (char* str)
	{
		if (_myRank==0)
		{
			cerr << str;
			outfile << str;
		}
		exit(0);
	}
}

inline bool pdDataManager::FindKeyword(const string keywd, const int keywdType)
{
	// This function checks if the desired keyword is shown in the input file and returns a bool value
	// If the keyword is required (keywdType=REQ), fail to find the keyword will stop the program
	// If the keyword is optional (keywdType=OPT), just return the search result

	ifstream infile(_inFilePath.c_str(), ios::in);
	ofstream outfile(_outFilePath.c_str(), ios::app);
	try 
	{
		// check if keyword length is correct
		if (keywd.empty())
		{
			throw "FindKeyword: keyword empty!\n";
		}
		else
		{
			// start the search
			string line;
			do 
			{
				getline(infile, line);
				// erase the space at two ends of the read-in string
				line.erase(line.find_last_not_of(' ')+1, string::npos);
				line.erase(0, line.find_first_not_of(' '));
				if (line==keywd)  // keyword found
				{
					infile.close();
					outfile.close();
					return true;
				}
			}
			while (!infile.eof());
			// keyword not found
			if (keywdType==REQ)
			{
				throw "FindKeyword: required keyword " + keywd + " not found!\n";
			}
			else
			{
				infile.close();
				outfile.close();
				return false;
			}
		}
	}
	catch (char* str)
	{
		if (_myRank==0)
		{
			cerr << str;
			outfile << str;
		}
		exit(0);
	}
}

#endif