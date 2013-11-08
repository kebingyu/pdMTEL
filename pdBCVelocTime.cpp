#include "pdBCVelocTime.h"

pdBCVelocTime::pdBCVelocTime() {};

pdBCVelocTime::pdBCVelocTime(const pdBCVelocTime &vbc) {};

pdBCVelocTime::pdBCVelocTime(int id, int shape, int x1, int x2, int x3, vector<double> v, double t1, double t2, double tend)
: pdBdryCondition(id, x1, x2, x3, v, tend), _shape(shape), _t1(t1), _t2(t2)
{
	_type = VELOC_TIME;
}

pdBCVelocTime::~pdBCVelocTime() {}

double pdBCVelocTime::GetFactor(double time) const
{
	if (_shape==1)
	{
		if (time<=_t1)
		{
			return 1.;
		}
		else if (time<=_t2)
		{
			return -1.;
		}
		else
		{
			return 0.;
		}
	}
	else
	{
		return 1.;
	}
}

void pdBCVelocTime::Print(ostream & os) const
{
	pdBdryCondition::Print(os);
	os << "  Set Point = " << _t1 << "\t" << _t2 << endl;
	os << "  Type: Prescribed Time-dependant Velocity" << endl;
}