/*

This class implement a traction force boundary condition to boundary region.
The affected boundary region will have constant force (_f1,_f2,_f3).
This condition is applied only in the directions xJ for which _dirxJ is Set to 1, where J=1,2,3.
If _dirxJ is Set to 0, then the node motion in the xJ direction is determined by the equation of motion. 

keyword format in p_input.txt:
	boundary_condition_(id)
	3 (flag_x1) (flag_x2) (flag_x3) (value_x1) (value_x2) (value_x3) (en_time)

	note: id = boundary condition id, number starts from 1 

*/

#pragma once
#include "pdBdryCondition.h"

class pdBCTract : public pdBdryCondition
{
public:
	// Consturctor
	pdBCTract(int id, int x1, int x2, int x3, vector<double> value, double tend);
	
	// Destructor
	virtual ~pdBCTract();
	
	// Function


private:
	// Data
	int _numnd; // number of nodes in this boundary condition region	
	
	// Never used constructor and copy control
	pdBCTract();
	pdBCTract(const pdBCTract&);	
	pdBCTract& operator=(const pdBCTract&);

	// Function
	void Print(ostream&) const;
};