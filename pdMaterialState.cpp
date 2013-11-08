#include "pdMaterialState.h"
#include "math.h"

pdMaterialState::pdMaterialState() {}

pdMaterialState::pdMaterialState(const pdMaterialState& mat) {}

pdMaterialState::pdMaterialState(int id, int type, double ymod, double sigmay, double nu, double ecrit, double denst)
	:pdMaterial(id, type, ymod, sigmay, ecrit, denst)
{
	_nu = nu;
}

pdMaterialState::~pdMaterialState() {}

void pdMaterialState::Print(ostream& os) const
{
	pdMaterial::Print(os);
	os << "  Type: State-Based Material" << endl;
}

