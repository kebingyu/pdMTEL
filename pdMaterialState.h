/* 

This class implements a State-Based peridynamic material.

keyword format in the input file:
	material_(id)
	2 (young's modulus) (yield strength) (poisson's ratio) (critical stretch) (density)

	note: id = material id, number starts from 1 

*/

#ifndef PDMATERIALSTATE_H
#define PDMATERIALSTATE_H

#include "pdMaterial.h"

class pdMaterialState : public pdMaterial
{
public:
	// Constructors		
	pdMaterialState(int id, int type, double ymod, double sigmay, double nu, double ecrit, double denst);

	// Destructor
	~pdMaterialState();
	
	// Function

private:
	// Data

	// Never used constructor and copy control
	pdMaterialState();
	pdMaterialState(const pdMaterialState&);	
	pdMaterialState& operator=(const pdMaterialState&);

	// Function
	void Print(ostream&) const;
};

#endif