#include "pdBdryConditionDispl.h"

pdBdryConditionDispl::pdBdryConditionDispl() {};

pdBdryConditionDispl::pdBdryConditionDispl(const pdBdryConditionDispl &dbc) {};

pdBdryConditionDispl::pdBdryConditionDispl(int id, int x1, int x2, int x3, vector<double> value, double tend, double loadType)
: pdBdryCondition(id, x1, x2, x3, value, tend, loadType)
{
	_type = DISPL;
}

pdBdryConditionDispl::~pdBdryConditionDispl() {}

double pdBdryConditionDispl::GetTimeFactor(double currentTime) const
{
	//(1) ramp loading
	if (_loadType==1.0)
	{
		return (currentTime <= _endTime) ? currentTime / _endTime: 1.0;
	}
	//(2) step loading
	else
	{
		return 1.0; 
	}
}

void pdBdryConditionDispl::Print(ostream & os) const
{
	pdBdryCondition::Print(os);
	os << "  Type: Prescribed Displacement" << endl;
}