#include "pdGauss.h"

pdGauss::pdGauss() {}

pdGauss::pdGauss(const pdGauss &) {}

pdGauss::~pdGauss() {}

void pdGauss::SetAbscissa(const double x, const double y, const double z)
{
	_abscissaX1 = x;
	_abscissaX2 = y;
	_abscissaX3 = z;
}

void pdGauss::SetWeight(const double x, const double y, const double z)
{
	_weightX1 = x;
	_weightX2 = y;
	_weightX3 = z;
}

double pdGauss::GetAbscissaX1() const
{
	return _abscissaX1;
}

double pdGauss::GetAbscissaX2() const
{
	return _abscissaX2;
}

double pdGauss::GetAbscissaX3() const
{
	return _abscissaX3;
}

double pdGauss::GetWeightX1() const
{
	return _weightX1;
}

double pdGauss::GetWeightX2() const
{
	return _weightX2;
}

double pdGauss::GetWeightX3() const
{
	return _weightX3;
}