#include "pdStress.h"

pdStress::pdStress(const pdStress& sigma) {}

pdStress::pdStress()
{
	s_rot_tens = new double* [3];
	s_left_stretch = new double* [3];
	s_stress_cauchy = new double* [3];
	s_stress_piola = new double* [3];
	s_strain_cauchy = new double* [3];
	for (int i=0;i!=3;++i)
	{
		s_rot_tens[i] = new double [3];
		s_left_stretch[i] = new double [3];
		s_stress_cauchy[i] = new double [3];
		s_stress_piola[i] = new double [3];
		s_strain_cauchy[i] = new double [3];
	}
}

pdStress::~pdStress()
{
	for (int i=0;i!=3;++i)
	{
		delete s_rot_tens[i];
		delete s_left_stretch[i];
		delete s_stress_cauchy[i];
		delete s_stress_piola[i];
		delete s_strain_cauchy[i];
	}
	delete [] s_rot_tens;
	delete [] s_left_stretch;
	delete [] s_stress_cauchy;
	delete [] s_stress_piola;
	delete [] s_strain_cauchy;
}

void pdStress::Initial()
{
	// These tensors are initialized at the first cycle in "InitNode" function.
	for (int i=0;i!=3;++i)
	{
		for (int j=0;j!=3;++j)
		{
			s_rot_tens[i][j] = (i==j);
			s_left_stretch[i][j] = (i==j);
			s_stress_cauchy[i][j] = 0.;
			s_stress_piola[i][j] = 0.;
			s_strain_cauchy[i][j] = 0.;
		}
	}
}

void pdStress::Set_rot_tens(int i, int j, double s)
{
	s_rot_tens[i][j] = s;
}

void pdStress::Set_left_stretch(int i, int j, double s)
{
	s_left_stretch[i][j] = s;
}

void pdStress::Set_stress_cauchy(int i, int j, double s)
{
	s_stress_cauchy[i][j] = s;
}

void pdStress::Set_stress_piola(int i, int j, double s)
{
	s_stress_piola[i][j] = s;
}

void pdStress::Set_strain_cauchy(int i, int j, double s)
{
	s_strain_cauchy[i][j] = s;
}

double pdStress::Get_rot_tens(int i, int j)
{
	return s_rot_tens[i][j];
}

double pdStress::Get_left_stretch(int i, int j)
{
	return s_left_stretch[i][j];
}

double pdStress::Get_stress_cauchy(int i, int j)
{
	return s_stress_cauchy[i][j];
}

double pdStress::Get_stress_piola(int i, int j)
{
	return s_stress_piola[i][j];
}

double pdStress::Get_strain_cauchy(int i, int j)
{
	return s_strain_cauchy[i][j];
}

void pdStress::Print(ostream& os) const
{
	os << "stress cauchy:" << endl;
	for (int k=0;k!=3;++k)
	{
		for (int l=0;l!=3;++l)
		{
			os << "\t" << s_stress_cauchy[k][l];
		}
		os << endl;
	}
	os << "stress piola:" << endl;
	for (int k=0;k!=3;++k)
	{
		for (int l=0;l!=3;++l)
		{
			os << " \t" << s_stress_piola[k][l];
		}
		os << endl;
	}
	os << "strain cauchy:" << endl;
	for (int k=0;k!=3;++k)
	{
		for (int l=0;l!=3;++l)
		{
			os << "\t" << s_strain_cauchy[k][l];
		}
		os << endl;
	}
}

ostream& operator << (ostream& os, const pdStress& s)
{
	s.Print(os);
	return os;
}
