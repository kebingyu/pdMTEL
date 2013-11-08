/*

This is the base class for all types of spaces.

*/

#ifndef PDSPACE_H
#define PDSPACE_H

#include <iostream>
using std::ostream;
using std::endl;

class pdSpace
{
public:
	// Constructor
	pdSpace(int, int);

	// Destructor
	virtual ~pdSpace();

	// Interface
	void SetCuboidSpaceBdry(double x1lo, double x1hi, double x2lo, double x2hi, double x3lo, double x3hi);
	void SetCylinderSpaceBdry(double rad, double x1cen, double x2cen, double x3lo, double x3hi);
	int GetID() const;
	int GetGeomType() const;
	double GetX1Low() const;
	double GetX1High() const;
	double GetX2Low() const;
	double GetX2High() const;
	double GetRadius() const;
	double GetX1Center() const;
	double GetX2Center() const;
	double GetX3Low() const;
	double GetX3High() const;

	// Overload operator
	friend ostream& operator << (ostream&, const pdSpace&);

protected:
	// Data
	int _id; // space id
	int _geomType;	// space geometry type (3D): 1 for cuboid, 2 for cylinder
								// for 2D case (slit and plane): 1 for rectangle, 2 for circle
	double _x1Low, _x1High;	// space boundary for cuboid
	double _x2Low, _x2High;
	double _x3Low, _x3High;	// (share with cylinder type)
	double _radius;	//space boundary for cylinder
	double _x1Center, _x2Center; 

	// Never used constructor and copy control
	pdSpace();
	pdSpace(const pdSpace&);
	pdSpace& operator=(const pdSpace&);

	// Function
	virtual void Print(ostream&) const;
};

#endif