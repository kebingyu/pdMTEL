/*

This class implement a prescribed velocity boundary condition to boundary space.
The affected boundary space has constant velocity (_v1,_v2,_v3) for time<tend; and 
velocity = 0 thereafter.  This condition is applied only in the directions xJ for which _dirxJ is Set to 
1, where J=1,2,3.  If _dirxJ is Set to 0, then the node motion in the xJ direction is determined by the 
equation of motion. 

keyword format in the input file:
	boundary_condition_(id)
	2 (flag_x1) (flag_x2) (flag_x3) (value_x1) (value_x2) (value_x3) (end_time) (load_type)

	note: id = boundary condition id, number starts from 1 

*/

#ifndef PDBDRYCONDITIONVELOC_H
#define PDBDRYCONDITIONVELOC_H

#include "pdBdryCondition.h"

class pdBdryConditionVeloc : public pdBdryCondition
{
public:
	// Consturctor
	pdBdryConditionVeloc(int id, int x1, int x2, int x3, vector<double> value, double tend, double loadType);

	// Destructor
	~pdBdryConditionVeloc();

	// Function
	double GetTimeFactor(double) const;

private:
	// Never used constructor and copy control
	pdBdryConditionVeloc();
	pdBdryConditionVeloc(const pdBdryConditionVeloc&);
	pdBdryConditionVeloc& operator=(const pdBdryConditionVeloc&);

	// Function
	void Print(ostream&) const;
};

#endif