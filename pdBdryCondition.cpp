#include "pdBdryCondition.h"

pdBdryCondition::pdBdryCondition(int id, int f1, int f2, int f3, vector<double> v, double tend, double loadType)
	: _id(id), _dirx1(f1), _dirx2(f2), _dirx3(f3), _value(v), _endTime(tend), _loadType(loadType) {}

pdBdryCondition::pdBdryCondition() {}

pdBdryCondition::pdBdryCondition(const pdBdryCondition &bc) {}

pdBdryCondition::~pdBdryCondition() {}

int pdBdryCondition::GetID() const
{
	return _id;
}

pdBdryCondition::BCType pdBdryCondition::GetType() const
{
	return _type;
}

 // return direction flag in three directions
int pdBdryCondition::GetDirFlagX1() const
{
	return _dirx1;
}

int pdBdryCondition::GetDirFlagX2() const
{
	return _dirx2;
}

int pdBdryCondition::GetDirFlagX3() const
{
	return _dirx3;
}

double pdBdryCondition::GetValue(int i) const
{
	// this function returns one of the three components for this boundary condition value (i=1, 2, or 3)
	return _value[i - 1];
}

double pdBdryCondition::GetValue(int i, int j) const
{
	// this function returns one of the nine components for this boundary condition value (i, j=1, 2, or 3)
	return _value[3 * (i - 1) + (j - 1)];
}

double pdBdryCondition::GetEndTime() const
{
	return _endTime;
}

void pdBdryCondition::Print(ostream & os) const
{
	os << "Boundary Condition " << _id+1 << endl;
	os << "  Direction flag = " << _dirx1 << "\t" << _dirx2 << "\t" << _dirx3 << endl;
	os << "  Value = ";
	for (vector<double>::const_iterator itor=_value.begin();itor!=_value.end();++itor)
	{
		os << (*itor) << "\t";
	}
	os << endl << "  End time = " << _endTime << "\t";
	if (_loadType==1.0)
	{
		os << "Loading type = Ramp loading";
	}
	else
	{
		os << "Loading type = Step loading";
	}
	os << endl;
}

ostream& operator << (ostream& os, const pdBdryCondition& bc)
{
	bc.Print(os);
	return os;
}





