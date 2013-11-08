/*

This class implements an prescribed displacement gradient boundary condition to boundary space.
The affected boundary space will have constant displacement gradient (_e11 to _e33, 9 components).
This condition is applied only in the directions xJ for which _dirxJ is Set to 1, where J=1,2,3.
If _dirxJ is Set to 0, then the node motion in the xJ direction is determined by the equation of motion. 

keyword format in the input file:
	boundary_condition_(id)
	3 (flag_x1) (flag_x2) (flag_x3) (value_x11) (value_x12) (value_x13) (value_x21) (value_x22) (value_x23)(value_x31) (value_x32) (value_x33) (end_time) (load_type)

	note: id = boundary condition id, number starts from 1 

*/

#ifndef PDBDRYCONDITIONDISPLGRAD_H
#define PDBDRYCONDITIONDISPLGRAD_H

#include "pdBdryCondition.h"

class pdBdryConditionDisplGrad : public pdBdryCondition
{
public:
	// Consturctor
	pdBdryConditionDisplGrad(int id, int x1, int x2, int x3, vector<double> v, double tend, double loadType);
	
	// Destructor
	~pdBdryConditionDisplGrad();
	
	// Function
	double GetTimeFactor(double) const;

private:
	// Never used constructor and copy control
	pdBdryConditionDisplGrad();
	pdBdryConditionDisplGrad(const pdBdryConditionDisplGrad&);
	pdBdryConditionDisplGrad& operator=(const pdBdryConditionDisplGrad&);
	
	// Function
	void Print(ostream&) const;
};

#endif