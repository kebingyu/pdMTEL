#include "pdSpaceBdryCondition.h"
#include "pdBdryCondition.h"

pdSpaceBdryCondition::pdSpaceBdryCondition() {}

pdSpaceBdryCondition::pdSpaceBdryCondition(const pdSpaceBdryCondition &bspace) {}

pdSpaceBdryCondition::pdSpaceBdryCondition(int id, int type, pdBdryCondition* bc) : pdSpace(id, type), _bc(bc) {}

pdSpaceBdryCondition::~pdSpaceBdryCondition() 
{
	//if (_bc!=NULL)
	//	delete _bc;
}

pdBdryCondition* pdSpaceBdryCondition::GetBC() const
{
	return _bc;
}

void pdSpaceBdryCondition::Print(ostream& os) const
{
	pdSpace::Print(os);
	os << endl << "  Type: boundary condition space associates with BC id = " << _bc->GetID() << endl;
}