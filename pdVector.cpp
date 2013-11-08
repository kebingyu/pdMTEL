#include "pdVector.h"
#include "pdMatrix.h"

pdVector::pdVector(const pdVector &vec) {}

pdVector::pdVector() {}

pdVector::pdVector(const int myRow)
{
	_coeff.resize(myRow);
}

pdVector::~pdVector() 
{
	_coeff.clear(); 
}

void pdVector::Zero()
{
	// set all coefficient of the vector to be zero
	for (vector<double>::iterator itor=_coeff.begin();itor!=_coeff.end();++itor)
	{
		(*itor) = 0.0;
	}
}

int pdVector::GetNumRows() const
{
	return (int)_coeff.size();
}

void pdVector::SetCoeff(const int i, const double value)
{
	_coeff[i] = value;
}

void pdVector::AddCoeff(const int i, const double value)
{
	_coeff[i] += value;
}

double pdVector::GetCoeff(const int i) const
{
	return _coeff[i];
}

pdVector* pdVector::operator +(const pdVector& b)
{
	try 
	{
		int myRow = this->GetNumRows();
		int bRow = b.GetNumRows();
		// rank check
		if (myRow != bRow)
		{
			throw "pdVector::operator +: ranks do not match!\n";
		}
		// else continue
		pdVector* temp = new pdVector(myRow);
		for (int i=0; i!=myRow; ++i)
		{
			temp->SetCoeff(i, this->GetCoeff(i) + b.GetCoeff(i));
		}
		return temp;
	}
	catch (char* str)
	{
		cerr << str;
		exit(0);
	}
}

pdVector* pdVector::operator -(const pdVector& b)
{
	try 
	{
		int myRow = this->GetNumRows();
		int bRow = b.GetNumRows();
		// rank check
		if (myRow != bRow)
		{
			throw "pdVector::operator -: ranks do not match!\n";
		}
		// else continue
		pdVector* temp = new pdVector(myRow);
		for (int i=0; i!=myRow; ++i)
		{
			temp->SetCoeff(i, this->GetCoeff(i) - b.GetCoeff(i));
		}
		return temp;
	}
	catch (char* str)
	{
		cerr << str;
		exit(0);
	}
}

pdMatrix* pdVector::DyadicProd(pdVector* b)
{
	// c[i][j] = _coeff[i] * (b[j])^T, dyadic product
	int myRow = this->GetNumRows();
	int bRow = b->GetNumRows();
	pdMatrix* c = new pdMatrix(myRow, bRow);
	for (int i = 0; i != myRow; ++i)
	{
		for (int j = 0; j != bRow; ++j)
		{
			c->SetCoeff(i, j, _coeff[i] * b->GetCoeff(j));
		}
	}
	return c; 
}

double pdVector::Mult(pdVector* b)
{
	// c[i][j] = (_coeff[i])^T * b[j], dot product of two vectors
	try 
	{
		int myRow = this->GetNumRows();
		int bRow = b->GetNumRows();
		// rank check
		if (myRow != bRow)
		{
			throw "pdVector::Mult(pdVector*): ranks do not match!\n";
		}
		// else continue		
		double sum = 0.0;
		for (int i = 0; i != myRow; ++i)
		{
			sum += _coeff[i] * b->GetCoeff(i);
		}
		return sum;
	}
	catch (char* str)
	{
		cerr << str;
		exit(0);
	}
 }

void pdVector::Mult(double b)
{
	// vector multiply by a constant b
	int myRow = this->GetNumRows();
	for (int i=0; i!=myRow; ++i)
	{
		this->SetCoeff(i, _coeff[i] * b);
	}
}

double pdVector::Norm() const
{
	int myRow = this->GetNumRows();
    double sum = 0.0;
    for (int i = 0; i != myRow; ++i)
	{
        sum += _coeff[i] * _coeff[i];
	}
    return sqrt(sum);
}

void pdVector::Print(ostream &os) const
{
	int myRow = this->GetNumRows();
	for (int i=0; i!=myRow; ++i)
	{
		os << _coeff[i] << "\t";
	}
	os << endl;
}

ostream& operator << (ostream& os, const pdVector& v)
{
	v.Print(os);
	return os;
}

