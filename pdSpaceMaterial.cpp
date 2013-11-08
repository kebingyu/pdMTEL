#include "pdSpaceMaterial.h"
#include "pdMaterial.h"

pdSpaceMaterial::pdSpaceMaterial() {}

pdSpaceMaterial::pdSpaceMaterial(const pdSpaceMaterial &mspace) {}

pdSpaceMaterial::pdSpaceMaterial(int id, int type, pdMaterial* mat) : pdSpace(id, type), _mat(mat) {}

pdSpaceMaterial::~pdSpaceMaterial() 
{
	//if (_mat!=NULL)
	//	delete _mat;
}

pdMaterial* pdSpaceMaterial::GetMaterial() const
{
	return _mat;
}

void pdSpaceMaterial::Print(ostream& os) const
{
	pdSpace::Print(os);
	os << endl << "  Type: material space associates with material id = " << _mat->GetID() << endl;
}