/*

This is the base material class.

*/

#ifndef PDMATERIAL_H
#define PDMATERIAL_H

#include <iostream>
using namespace std;

class pdMaterial
{
public:
	// Constructor
	pdMaterial(int id, int type, double ymod, double sigmay, double ecrit, double denst);

	// Destructor
	virtual ~pdMaterial();

	// Functions
	void SetSpringConstant(double);
	void SetGravity(double, double, double);
	int GetID() const;
	int GetType() const;
	double GetNu() const;
	double GetCriticalStretch() const;
	double GetYieldStrength() const;
	double GetDensity() const;
	double GetYoungModulus() const;
	double GetGrav1() const;
	double GetGrav2() const;
	double GetGrav3() const;
	virtual double GetBulkModulus() const;
	virtual double GetShearModulus() const;
	virtual double GetBulkSoundSpeed() const;
	virtual double GetSpringConstant() const;
	virtual double GetYieldStretch() const;

	// Overload operator
	friend ostream& operator << (ostream&, const pdMaterial&);

protected:
	// Data
	int _id; // Material id
	int _type; // Material type:
					// 1 for PMB material, 2 for state-based material
	double _nu;     // poisson's ratio
	double _ymod;   // Young's modulus
	double _sigmay; // yield strength
	double _critStretch;  // Critical stretch
	double _density;  // Material density
	double _spconst;   // spring constant		
	double _grav1, _grav2, _grav3; // material gravity, default=0 in three directions

	// Never used constructors and copy operator
	pdMaterial();
	pdMaterial(const pdMaterial&);
	pdMaterial& operator=(const pdMaterial&);

	// Function
	virtual void Print(ostream&) const;
};

#endif



