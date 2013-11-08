#include "pdMaterial.h"
#include "math.h"

pdMaterial::pdMaterial() {}

pdMaterial::pdMaterial(const pdMaterial& mat) {}

pdMaterial::pdMaterial(int id, int type, double ymod, double sigmay, double ecrit, double denst)
: _id(id), _type(type), _ymod(ymod), _sigmay(sigmay), _critStretch(ecrit), _density(denst) {}

pdMaterial::~pdMaterial() {}

void pdMaterial::SetGravity(double g1, double g2, double g3)
{
	_grav1 = g1;
	_grav2 = g2;
	_grav3 = g3;
}

int pdMaterial::GetID() const
{
	return _id;
}

int pdMaterial::GetType() const
{
	return _type;
}

double pdMaterial::GetNu() const
{
	return _nu;
}

double pdMaterial::GetCriticalStretch() const
{
	return _critStretch;
}

double pdMaterial::GetDensity() const
{
	return _density;
}

double pdMaterial::GetYoungModulus() const
{
	return _ymod;
}

double pdMaterial::GetGrav1() const
{
	return _grav1;
}

double pdMaterial::GetGrav2() const
{
	return _grav2;
}

double pdMaterial::GetGrav3() const
{
	return _grav3;
}

void pdMaterial::SetSpringConstant(double spconst)
{
	_spconst = spconst;
}

double pdMaterial::GetYieldStrength() const
{
	return _sigmay;
}

double pdMaterial::GetBulkModulus() const
{
	return _ymod / (3.0 * (1.0 - 2.0 * _nu));
}

double pdMaterial::GetShearModulus() const
{
	return _ymod / ((1.0 + _nu) * 2.0);
}

double pdMaterial::GetBulkSoundSpeed() const
{
	return sqrt(this->GetBulkModulus() / _density);
}

double pdMaterial::GetYieldStretch() const 
{
	return _sigmay / (2.0 * _ymod);
}

double pdMaterial::GetSpringConstant() const
{
	return _spconst;
}

void pdMaterial::Print(ostream& os) const
{
	os << "Material " << _id+1 << endl
		<< "  Young's modulus= " << _ymod << "\t" << "Yield strength= " << _sigmay << endl
		<< "  Poisson's ratio= " << _nu << "\t" << "Density= " << _density << endl
		<< "  Critical stretch= " << _critStretch << "\t" << "Spring constant= " << _spconst << endl;
}

ostream& operator << (ostream& os, const pdMaterial& m)
{
	m.Print(os);
	return os;
}

