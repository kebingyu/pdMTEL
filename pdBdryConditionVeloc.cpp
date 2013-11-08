#include "pdBdryConditionVeloc.h"

pdBdryConditionVeloc::pdBdryConditionVeloc() {};

pdBdryConditionVeloc::pdBdryConditionVeloc(const pdBdryConditionVeloc &vbc) {};

pdBdryConditionVeloc::pdBdryConditionVeloc(int id, int x1, int x2, int x3, vector<double> v, double tend, double loadType)
: pdBdryCondition(id, x1, x2, x3, v, tend, loadType) 
{
	_type = VELOC;
}

pdBdryConditionVeloc::~pdBdryConditionVeloc() {}

double pdBdryConditionVeloc::GetTimeFactor(double time) const
{	
	return (time <= _endTime) ? 1.0: 0.0;
}

void pdBdryConditionVeloc::Print(ostream & os) const
{
	pdBdryCondition::Print(os);
	os << "  Type: Prescribed Velocity" << endl;
}