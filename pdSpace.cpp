#include "pdSpace.h"

pdSpace::pdSpace() {}

pdSpace::pdSpace(const pdSpace &reg) {}

pdSpace::pdSpace(int id,  int type) : _id(id), _geomType(type) {}

pdSpace::~pdSpace() {}

void pdSpace::SetCuboidSpaceBdry(double x1lo, double x1hi, 
		double x2lo, double x2hi,double x3lo, double x3hi)
{
	_x1Low = x1lo;
	_x1High = x1hi;
	_x2Low = x2lo;
	_x2High = x2hi;
	_x3Low = x3lo;
	_x3High = x3hi;
}

void pdSpace::SetCylinderSpaceBdry(double rad, double x1cen, 
		double x2cen, double x3lo, double x3hi)
{
	_radius = rad;
	_x1Center = x1cen;
	_x2Center = x2cen;
	_x3Low = x3lo;
	_x3High = x3hi;
}

int pdSpace::GetID() const
{
	return _id;
}

int pdSpace::GetGeomType() const
{
	return _geomType;
}

double pdSpace::GetX1Low() const
{
	return _x1Low;
}

double pdSpace::GetX1High() const
{
	return _x1High;
}

double pdSpace::GetX2Low() const
{
	return _x2Low;
}

double pdSpace::GetX2High() const
{
	return _x2High;
}

double pdSpace::GetRadius() const
{
	return _radius;
}

double pdSpace::GetX1Center() const
{
	return _x1Center;
}

double pdSpace::GetX2Center() const
{
	return _x2Center;
}

double pdSpace::GetX3Low() const
{
	return _x3Low;
}

double pdSpace::GetX3High() const
{
	return _x3High;
}

void pdSpace::Print(ostream& os) const
{
	os << "Space " << _id+1 << endl;
	if (_geomType==1)
	{
		os << "  x1= " << _x1Low << "\t" << _x1High << "\t" 
			<< " x2= " << _x2Low << "\t" << _x2High << "\t" << " x3= " << _x3Low << "\t" << _x3High;
	}
	else 
	{
		os << "  center= " << _x1Center << "\t" << _x2Center << "\t" << " radius= " << _radius 
			<< "\t" << " x3= " << _x3Low << "\t" << _x3High;
	}
}

ostream& operator << (ostream& os, const pdSpace& s)
{
	s.Print(os);
	return os;
}