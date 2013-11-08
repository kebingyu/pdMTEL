/*

This class implement a prescribed displacment boundary condition to boundary space.
The affected boundary space will have a total displacement of _uJ from zero to end time.
This condition is applied only in the directions xJ for which _dirxJ is Set to 1, where J=1,2,3.
If _dirxJ is Set to 0, then the node motion in the xJ direction is determined by the equation of motion. 

keyword format in the input file:
	boundary_condition_(id)
	1 (flag_x1) (flag_x2) (flag_x3) (value_x1) (value_x2) (value_x3) (end_time) (load_type)

	note: id = boundary condition id, number starts from 1 

*/

#ifndef PDBDRYCONDITIONDISPL_H
#define PDBDRYCONDITIONDISPL_H

#include "pdBdryCondition.h"

class pdBdryConditionDispl : public pdBdryCondition
{
public:
	// Consturctor
	pdBdryConditionDispl(int id, int x1, int x2, int x3, vector<double> value, double tend, double loadType);
	
	// Destructor
	~pdBdryConditionDispl();
	
	// Function
	double GetTimeFactor(double) const;

private:
	// Never used constructor and copy control
	pdBdryConditionDispl();
	pdBdryConditionDispl(const pdBdryConditionDispl&);	
	pdBdryConditionDispl& operator=(const pdBdryConditionDispl&);

	// Function
	void Print(ostream&) const;
};

#endif