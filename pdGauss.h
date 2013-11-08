/* 

   This class implments a Gaussian-Legendre integration point. 
   It is a sub class of pdIntgPoint.
   It also contains the abscissas and weights in three directions.
*/

#ifndef PDGAUSS_H
#define PDGAUSS_H

#include "pdIntgPoint.h"

class pdGauss : public pdIntgPoint
{
public:
	// Constructor
	pdGauss();

	// Destructor
	~pdGauss();

	// Function
	void SetAbscissa(const double, const double, const double);
	void SetWeight(const double, const double, const double);
	double GetAbscissaX1() const;
	double GetAbscissaX2() const;
	double GetAbscissaX3() const;
	double GetWeightX1() const;
	double GetWeightX2() const;
	double GetWeightX3() const;

private:
	// Data
	double _abscissaX1, _abscissaX2, _abscissaX3;
	double _weightX1, _weightX2, _weightX3;	

	// Never used copy control
	pdGauss(const pdGauss&);
	pdGauss& operator=(const pdGauss&);
};

#endif