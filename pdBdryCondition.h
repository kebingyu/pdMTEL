/*

This class is the base class for boundary condtion type.

Boundary condition type:
  1: Prescribed displacement for t<total_time
  2: Prescribed velocity for t<end_time
  3: Prescribed displacement gradient (prescribed strain) for t<end_time
  4: Traction boundary condition

Boundary condition loading type:
  1: Ramp loading
  2: Step loading

*/

#ifndef PDBDRYCONDITION_H
#define PDBDRYCONDITION_H

#include <iostream>
#include <vector>
#include "math.h"
using namespace std;

class pdBdryCondition
{
public:
	enum BCType { DISPL=1, VELOC, DISPLGRAD, TRACTION };
	//enum LoadType { RAMP=1, STEP };

	// Constructor
	pdBdryCondition(int id, int f1, int f2, int f3, vector<double> value, double tend, double loadType);

	// Destructor
	virtual ~pdBdryCondition();

	// Interface
	int GetID() const;
	BCType GetType() const;
	int GetDirFlagX1() const;
	int GetDirFlagX2() const;
	int GetDirFlagX3() const;	
	double GetValue(int) const;
	double GetValue(int, int) const;
	double GetEndTime() const; 
	virtual double GetTimeFactor(double) const = 0;  // return a time-dependant factor

	// Overload operator
	friend ostream& operator << (ostream&, const pdBdryCondition&);	

protected:
	// Data
	int _id;   // boundary condition id, number starts from 0
	BCType _type; // boundary conditon type
	double _loadType; // the type boundary condition is loaded.
	int _dirx1, _dirx2, _dirx3; // direction flag, either 0 or 1
	vector<double> _value; // boundary condition value
	double _endTime; // the end time when this boundary condition does not apply

	// Never used constructor and copy control
	pdBdryCondition();
	pdBdryCondition(const pdBdryCondition&);
	pdBdryCondition& operator=(const pdBdryCondition&);

	// Function
	virtual void Print(ostream&) const;
};

#endif