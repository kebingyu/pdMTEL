/*

This class implements a space with boundary condition.

keyword format in the input file:
	boundary_space_(id)
	(geoType) (bcid) (x1_lo) (x1_hi) (x2_lo) (x2_hi) (x3_lo) (x3_hi)

note: 
id = boundary condition space id, number starts from 1 
bcid = boundary condition id associated with this space, number starts from 1

*/

#ifndef PDSPACEBDRYCONDITION_H
#define PDSPACEBDRYCONDITION_H

#include "pdSpace.h"

class pdBdryCondition;

class pdSpaceBdryCondition : public pdSpace
{
public:
	// Consturctor
	pdSpaceBdryCondition(int, int, pdBdryCondition*);

	// Destructor
	~pdSpaceBdryCondition();

	// Function
	pdBdryCondition* GetBC() const;

private:
	// Data
	pdBdryCondition* _bc; //boundary condition associated with this boundary space

	// Never used constructor and copy control
	pdSpaceBdryCondition();
	pdSpaceBdryCondition(const pdSpaceBdryCondition&);
	pdSpaceBdryCondition& operator=(const pdSpaceBdryCondition&);

	// Function
	void Print(ostream&) const;
};

#endif