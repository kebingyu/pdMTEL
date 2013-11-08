/*

The pdSolve class manages the solving procedure and prints out results.

*/

#ifndef PDSOLVE_H
#define PDSOLVE_H

#include <fstream>
#include <vector>
#include "pdIntgManager.h"
using std::ofstream;
using std::vector;

class pdDataManager;

class pdSolve
{
public:
	// Function
	static void InitialProblem(pdDataManager&);
	static void SolveDynamic(pdDataManager&);
	static void SolveStatic(pdDataManager&);
	static void Debug(pdDataManager&);
	static void WriteDetail(const pdDataManager&, const int, const double);
	static void WriteTecplot(const pdDataManager&, const int, const double);
	static void WriteEnergy(const pdDataManager&, const int, const double);
	static void WriteEachStep(const pdDataManager&, const int, const double, ofstream&);
	static void PostFormatData(const pdDataManager&);
	static void PostAFD(pdDataManager&);

private:
	// Constructor
	pdSolve();

	// copy control
	pdSolve(const pdSolve&);
	pdSolve& operator=(const pdSolve&);

	// Destructor
	~pdSolve();

	// Function
	static void AssembleEquation(pdDataManager&, const double time, vector<vector<double>>& matrix, vector<double>& load);
	static void Gauss(pdDataManager&, vector<vector<double>> m, vector<double> rhs, vector<double>& sol);
};

#endif
