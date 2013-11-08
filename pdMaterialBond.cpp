#include "pdMaterialBond.h"

pdMaterialBond::pdMaterialBond() {}

pdMaterialBond::pdMaterialBond(const pdMaterialBond& mat) {}

pdMaterialBond::pdMaterialBond(int id, int type, double ymod, double sigmay, double ecrit, double denst)
: pdMaterial(id, type, ymod, sigmay, ecrit, denst)
{
	_nu = 0.25;
}

pdMaterialBond::~pdMaterialBond() {}

void pdMaterialBond::Print(ostream& os) const
{
	pdMaterial::Print(os);
	os << "  Type: Bond-Based Material (PMB) " << endl;
}

