#include "pdBdryConditionDisplGrad.h"
using std::vector;

pdBdryConditionDisplGrad::pdBdryConditionDisplGrad() {};

pdBdryConditionDisplGrad::pdBdryConditionDisplGrad(const pdBdryConditionDisplGrad &idg) {};

pdBdryConditionDisplGrad::pdBdryConditionDisplGrad(int id, int x1, int x2, int x3, vector<double> v, double tend, double loadType)
: pdBdryCondition(id, x1, x2, x3, v, tend, loadType)
{
	_type = DISPLGRAD;
}

pdBdryConditionDisplGrad::~pdBdryConditionDisplGrad() {}

double pdBdryConditionDisplGrad::GetTimeFactor(double currentTime) const
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

void pdBdryConditionDisplGrad::Print(ostream & os) const
{
	pdBdryCondition::Print(os);
	os << "  Type: Prescribed Displacement Gradient " << endl;
}