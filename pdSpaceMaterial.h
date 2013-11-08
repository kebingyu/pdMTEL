/*

This class implements a space with material properties.

keyword format in the input file:
	material_space_(id)
	(geometry_type) (mid) (x1_lo) (x1_hi) (x2_lo) (x2_hi) (x3_lo) (x3_hi) 
	
note: 
id = material space id, number starts from 1 
mid = material id associated with this space, number starts from 1

*/

#ifndef PDSPACEMATERIAL_H
#define PDSPACEMATERIAL_H

#include "pdSpace.h"

class pdMaterial;

class pdSpaceMaterial : public pdSpace
{
public:
	// Consturctor
	pdSpaceMaterial(int, int, pdMaterial*);

	// Destructor
	~pdSpaceMaterial();

	// Function
	pdMaterial* GetMaterial() const;		

private:
	// Data
	pdMaterial* _mat; //material associated with this material space

	// Never used constructor and copy control
	pdSpaceMaterial();
	pdSpaceMaterial(const pdSpaceMaterial&);	
	pdSpaceMaterial& operator=(const pdSpaceMaterial&);

	// Function
	void Print(ostream&) const;
};

#endif