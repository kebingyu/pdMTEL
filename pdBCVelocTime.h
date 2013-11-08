/*

This class implement a time-dependant velocity boundary condition to boundary region.
The affected boundary region has constant velocity (_v1,_v2,_v3) for time<tend; and 
velocity = 0 thereafter.  This condition is applied only in the directions xJ for which _dirxJ is Set to 
1, where J=1,2,3.  If _dirxJ is Set to 0, then the node motion in the xJ direction is determined by the 
equation of motion. 

keyword format in p_input.txt:
	boundary_condition_(id)
	4 (shape) (flag_x1) (flag_x2) (flag_x3) (value_x1) (value_x2) (value_x3) (t1) (t2) (tend)

	note: id = boundary condition id, number starts from 1

*/

#pragma once
#include "pdBdryCondition.h"

class pdBCVelocTime : public pdBdryCondition
{
public:
	// Consturctor
	pdBCVelocTime(int id, int shape, int x1, int x2, int x3, vector<double> value, double t1, double t2, double tend);
	
	// Destructor
	virtual ~pdBCVelocTime();
	
	//function
	double GetFactor(double time) const;	

private:
	// Data
	int _shape; // 1 for a positive square pulse of v from 0<t<t1 followed by a negative square pulse of -v from t1<t<t2 
	double _t1, _t2; // time set points
	
	// Never used constructor and copy control
	pdBCVelocTime();
	pdBCVelocTime(const pdBCVelocTime&);
	pdBCVelocTime& operator=(const pdBCVelocTime&);

	// Function
	void Print(ostream&) const;
};
