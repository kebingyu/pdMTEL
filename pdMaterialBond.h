/* 

This class implements an bond-based peridyanmic material (Prototype Microelastic Brittal material).

keyword format in the input file:
	material_(id)
	1 (young's modulus) (yield strength) (critical stretch) (density)

	note: id = material id, number starts from 1 

*/

#ifndef PDMATERIALBOND_H
#define PDMATERIALBOND_H

#include "pdMaterial.h"

class pdMaterialBond : public pdMaterial
{
public:
	// Constructor		
	pdMaterialBond(int id, int type, double ymod, double sigmay, double ecrit, double denst);
	
	// Destructor
	~pdMaterialBond();

	// Function

private:
	// Data

	// Never used constructor and copy control
	pdMaterialBond();
	pdMaterialBond(const pdMaterialBond&);	
	pdMaterialBond& operator=(const pdMaterialBond&);	

	// Function
	void Print(ostream&) const;
};

#endif