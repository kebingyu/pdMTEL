#include "pdBCTract.h"

pdBCTract::pdBCTract() {};

pdBCTract::pdBCTract(const pdBCTract &dbc) {};

pdBCTract::pdBCTract(int id, int x1, int x2, int x3, vector<double> v, double tend)
: pdBdryCondition(id, x1, x2, x3, v, tend) 
{
	_type = TRACTION;
}

pdBCTract::~pdBCTract() {}

void pdBCTract::Print(ostream & os) const
{
	pdBdryCondition::Print(os);
	os << "  Type: Traction Force" << endl;
}