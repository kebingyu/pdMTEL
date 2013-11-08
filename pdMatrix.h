/* pdMatrix

The pdMatrix class is for a full matrix. It stores the numbers 
of rows and columns and all the coefficients of the matrix. It also has 
member functions that perform matrix and vector multiplications, and
solution of linear system equations using the Gaussian elimination.
The matrix implemented can be unsymmetric and full.

*/

#ifndef PDMATRIX_H
#define PDMATRIX_H

#include <iostream>
#include <cstdio>
#include "math.h"
#include "pdVector.h"
using namespace std;

class pdMatrix
{		
public:
	// Constructors	
    pdMatrix(const int myRow, const int myCol);	//construct a myRow by myCol matrix

    // Destructor
    ~pdMatrix();

    // Functions	
	void SetCoeff(const int, const int, const double);
	void AddCoeff(const int, const int, const double);
	double GetCoeff(const int, const int) const;
	int GetNumRows() const;	
	int GetNumCols() const;
	void AddColumn(pdVector*);
	void Zero();
	pdMatrix* Mult(const pdMatrix*);
	pdVector* Mult(const pdVector*);
	void Mult(const double);
    void Trans(pdMatrix*);
	void Gauss(pdVector*, pdVector*);
	/*
	This function solves a linear system of equations [_coeff]*{x} = {b}.  
	The stiffness matrix is the matrix data member _coeff of this class, the vector data member of 
	class b is the RHS vector, and the vector data member of class x will store 
	the solution. The Gauss() function uses a Gauss elimination with full pivoting;
	it works for any full, unsymmetric matrix.
	*/	
	double RowSumNorm();
	double FBnorm();

	// Overload operator
	pdMatrix* operator +(const pdMatrix&);
	pdMatrix* operator -(const pdMatrix&);
	friend ostream& operator << (ostream&, const pdMatrix&);	

 private:
	 // Data
	 vector<pdVector*> _coeff;

	// Never to be used constructors
	 pdMatrix();
	 pdMatrix(const pdMatrix&);
	 void Print(ostream&) const;
};

#endif
