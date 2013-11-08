#include "pdIntgPoint.h"

pdIntgPoint::pdIntgPoint() {}

pdIntgPoint::pdIntgPoint(const pdIntgPoint &) {}

pdIntgPoint::~pdIntgPoint() {}

/*
 (07/04/2011) The reason I am not using void SetX1(const double) is because this function is called
 in pdIntgManager::CalcIntgX3 to set up trapezoidal point coordinates and the arguments can be changed.
*/
void pdIntgPoint::SetX1(const double x)
{
	_x1 = x;
}

void pdIntgPoint::SetX2(double x)
{
	_x2 = x;
}

void pdIntgPoint::SetX3(double x)
{
	_x3 = x;
}

void pdIntgPoint::SetU1(const double u)
{
	_u1 = u;
}

void pdIntgPoint::SetU2(const double u)
{
	_u2 = u;
}

void pdIntgPoint::SetU3(const double u)
{
	_u3 = u;
}

double pdIntgPoint::GetX1() const
{
	return _x1;
}

double pdIntgPoint::GetX2() const
{
	return _x2;
}

double pdIntgPoint::GetX3() const
{
	return _x3;
}

double pdIntgPoint::GetU1() const
{
	return _u1;
}

double pdIntgPoint::GetU2() const
{
	return _u2;
}

double pdIntgPoint::GetU3() const
{
	return _u3;
}

