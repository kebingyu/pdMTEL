/* 

The pdVector class stores a column vector.

*/

#ifndef PDVECTOR_H
#define PDVECTOR_H

#include <vector>
#include <iostream>
#include "math.h"
using namespace std;
class pdMatrix;

class pdVector
{
	// This class implements a column vector.

public:
	// Constructors
	pdVector(const int);

    // Destructor
    ~pdVector();

    // Functions	
	void Zero();
	int GetNumRows() const;
	void SetCoeff(const int, const double);
	void AddCoeff(const int, const double);
	double GetCoeff(const int) const;
	pdMatrix* DyadicProd(pdVector*);
	double Mult(pdVector*);
	void Mult(const double);
	double Norm() const;	

	// Overload operator
	pdVector* operator +(const pdVector&);
	pdVector* operator -(const pdVector&);
	friend ostream& operator << (ostream&, const pdVector&);

private:
	// Data
    vector<double> _coeff;		// Coefficients of the vector

	// Never to be used constructor
	pdVector();
	pdVector(const pdVector&);
	void Print(ostream&) const;
};

#endif