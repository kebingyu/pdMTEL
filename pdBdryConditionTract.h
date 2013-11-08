/*

This class implement a traction force boundary condition to boundary region.
It follows the idea of 'general areal force density' (07/19/2011).

Current implementation is a traction force on the boundary surface without tangential components. 
  1. The affected boundary region will have constant traction force (_f1,_f2,_f3).
  2. This condition is applied only in one normal directions xJ for which _dirxJ is Set to 1, where J=1,2,3.
  3. The surface applied with traction force is the whole boundary surface.

keyword format in the input file:
	boundary_condition_(id)
	4 (flag_x1) (flag_x2) (flag_x3) (value_x1) (value_x2) (value_x3) (en_time) (load_type)

	note: id = boundary condition id, number starts from 1

*/

#ifndef PDBDRYCONDITIONTRACT_H
#define PDBDRYCONDITIONTRACT_H

#include "pdBdryCondition.h"

class pdBdryConditionTract : public pdBdryCondition
{
public:
	// Consturctor
	pdBdryConditionTract(int id, int x1, int x2, int x3, vector<double> value, double tend, double loadType);
	
	// Destructor
	~pdBdryConditionTract();
	
	// Function
	double GetTimeFactor(double) const;

private:
	// Data
	int _numnd; // number of nodes in this boundary condition region	
	
	// Never used constructor and copy control
	pdBdryConditionTract();
	pdBdryConditionTract(const pdBdryConditionTract&);	
	pdBdryConditionTract& operator=(const pdBdryConditionTract&);

	// Function
	void Print(ostream&) const;
};

#endif