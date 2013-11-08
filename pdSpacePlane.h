/*

This is a 2D subclass of pdSpace. It can be rectangle (type1) or circle (type2). 
Right now it only allows for plane who is normal to one of the three axes. 

*/

#ifndef PDSPACEPLANE_H
#define PDSPACEPLANE_H
#include "pdSpace.h"

class pdSpacePlane : public pdSpace
{
public:
	// Constructor
	pdSpacePlane(int, int);

	// Destructor
	~pdSpacePlane();

	// Function
	void SetForce(double);
	void AddForce(double);
	int GetNorm();
	double GetForce() const;
	double GetArea() const;

private:
	// Data
	int _norm; // norm direction (1, 2, or 3)
	double _force; // the total force across the plane in the norm direction
	
	// Never used constructors and copy control
	pdSpacePlane();
	pdSpacePlane(const pdSpacePlane&);
	pdSpacePlane& operator=(const pdSpacePlane&);

	// Function
	void Print(ostream&) const;
};

#endif