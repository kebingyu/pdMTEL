#include "pdSpacePlane.h"
#include "math.h"

enum PlaneGeometryType {RECTANGLE=1, CIRCLE};

pdSpacePlane::pdSpacePlane() {}

pdSpacePlane::pdSpacePlane(const pdSpacePlane &p) {}

pdSpacePlane::pdSpacePlane(int id, int type) : pdSpace(id, type) {}

pdSpacePlane::~pdSpacePlane() {}

void pdSpacePlane::SetForce(double f)
{
	_force = f;
}

void pdSpacePlane::AddForce(double f)
{
	_force += f;
}

double pdSpacePlane::GetForce() const
{
	return _force;
}

int pdSpacePlane::GetNorm()
{
	if (_geomType = RECTANGLE)
	{
		if (_x1Low==_x1High)
		{
			_norm = 1;
		}
		else if (_x2Low==_x2High)
		{
			_norm = 2;
		}
		else if (_x3Low==_x3High)
		{
			_norm = 3;
		}
	}
	else if (_geomType = CIRCLE)
	{
		_norm = 3;
	}
	return _norm;
}

double pdSpacePlane::GetArea() const
{
	if (abs(_norm)==1)
	{
		return (_x2High-_x2Low)*(_x3High-_x3Low);
	}
	else if (abs(_norm)==2)
	{
		return (_x1High-_x1Low)*(_x3High-_x3Low);
	}
	else
	{
		return (_x1High-_x1Low)*(_x2High-_x2Low);
	}
}

void pdSpacePlane::Print(ostream& os) const
{
	pdSpace::Print(os);
	os << "\t" << "Type: 2D plane, norm=" << "\t" << _norm << "\t" << "force_in_norm=" << "\t" << _force;
	os << "\t" << "areal_force_density=" << "\t" << _force/GetArea() << endl;
}
