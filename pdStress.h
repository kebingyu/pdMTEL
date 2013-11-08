/*

This class stores stress and related value for state-based peridynamic force calculation.

*/

#pragma once
#include <iostream>
using namespace std;

class pdStress
{	
public:
	// Constructors
	pdStress();

	// Destructor
	virtual ~pdStress();

	// Functions
	void Initial();
	void Set_rot_tens(int, int, double);
	void Set_left_stretch(int, int, double);
	void Set_stress_cauchy(int, int, double);
	void Set_stress_piola(int, int, double);
	void Set_strain_cauchy(int, int, double);
	double Get_rot_tens(int, int);
	double Get_left_stretch(int, int);
	double Get_stress_cauchy(int, int);
	double Get_stress_piola(int, int);
	double Get_strain_cauchy(int, int);

	// Overload operator
	friend ostream& operator << (ostream&, const pdStress&);

private:
	// Data
	double** s_rot_tens;
	double** s_left_stretch;
	double** s_stress_cauchy;
	double** s_stress_piola;
	double** s_strain_cauchy;

	// Never used constructors and copy control
	pdStress(const pdStress&);
	pdStress& operator=(const pdStress&);

	// Function
	void Print(ostream&) const;
};