#include "pdMatrix.h"


pdMatrix::pdMatrix() {}

pdMatrix::pdMatrix(const pdMatrix &m) {}

pdMatrix::pdMatrix(const int myRow, const int myCol) 
{
	for (int i=0; i!=myCol; ++i)
	{
		pdVector* column = new pdVector(myRow);
		//_coeff.resize(myCol);
		_coeff.push_back(column);
	}
}

pdMatrix::~pdMatrix() 
{
	for (vector<pdVector*>::iterator itor=_coeff.begin();itor!=_coeff.end();++itor)
	{
		if (*itor!=NULL)
			delete *itor;
	}
	_coeff.clear();
}

void pdMatrix::SetCoeff(const int i, const int j, const double value)
{
	//set value to the compnent at ith row, jth column
	_coeff[j]->SetCoeff(i, value);
}

void pdMatrix::AddCoeff(const int i, const int j, const double value) 
{
	//add value to the compnent at ith row, jth column
	_coeff[j]->AddCoeff(i, value);
}

double pdMatrix::GetCoeff(const int i, const int j) const
{
	//return the compnent at ith row, jth column
	return _coeff[j]->GetCoeff(i);
}

void pdMatrix::AddColumn(pdVector* col)
{
	// need to be used carefully 
	// after add new column, the number of matrix row and column increase
	_coeff.push_back(col);
}

int pdMatrix::GetNumRows() const 
{
	return _coeff[0]->GetNumRows();
}

int pdMatrix::GetNumCols() const 
{
	return (int)_coeff.size();
}

pdMatrix* pdMatrix::operator +(const pdMatrix& b)
{
	try 
	{
		int myRow = this->GetNumRows();
		int myCol = this->GetNumCols();
		int bRow = b.GetNumRows();
		int bCol = b.GetNumCols();
		// rank check
		if (myRow != bRow)
		{
			throw "pdMatrix::operator +: row ranks do not match!\n";
		}
		if (myCol != bCol)
		{
			throw "pdMatrix::operator +: column ranks do not match!\n";
		}
		// else continue
		pdMatrix* temp = new pdMatrix(myRow, myCol);
		for (int i = 0; i != myRow; ++i)
		{
			for (int j = 0; j != myCol; ++j)
			{
				temp->SetCoeff(i, j, this->GetCoeff(i, j) + b.GetCoeff(i, j));
			}
		}
		return temp;
	}
	catch (char* str)
	{
		cerr << str;
		exit(0);
	}
}

pdMatrix* pdMatrix::operator -(const pdMatrix& b)
{
    try 
	{
		int myRow = this->GetNumRows();
		int myCol = this->GetNumCols();
		int bRow = b.GetNumRows();
		int bCol = b.GetNumCols();
		// rank check
		if (myRow != bRow)
		{
			throw "pdMatrix::operator -: row ranks do not match!\n";
		}
		if (myCol != bCol)
		{
			throw "pdMatrix::operator -: column ranks do not match!\n";
		}
		// else continue
		pdMatrix* temp = new pdMatrix(myRow, myCol);
		for (int i = 0; i != myRow; ++i)
		{
			for (int j = 0; j != myCol; ++j)
			{
				temp->SetCoeff(i, j, this->GetCoeff(i, j) - b.GetCoeff(i, j));
			}
		}
		return temp;
	}
	catch (char* str)
	{
		cerr << str;
		exit(0);
	}
}

pdMatrix* pdMatrix::Mult(const pdMatrix *b)  
{
	// c[i][j] = (sum) _coeff[i][k]*b[k][j]   right multiplication of another matrix
	try 
	{
		int myRow = this->GetNumRows();
		int myCol = this->GetNumCols();
		int bRow = b->GetNumRows();
		int bCol = b->GetNumCols();
		// rank check
		if (myCol != bRow)
		{
			throw "pdMatrix::Mult(pdMatrix*): matrix ranks do not match!\n";
		}
		// else continue
		pdMatrix* c = new pdMatrix(myRow, bCol);
		for (int i=0; i!=myRow; ++i)
		{
			for (int j=0; j!=bCol; ++j) 
			{
				double sum = 0.0;
				for (int k=0; k!=myCol; ++k)
				{
					sum += this->GetCoeff(i, k) * b->GetCoeff(k, j);
				}
				c->SetCoeff(i, j, sum);
			}
		}
		return c;
	}
	catch (char* str)
	{
		cerr << str;
		exit(0);
	}
}

pdVector* pdMatrix::Mult(const pdVector *b)  
{
	// c[i] = (sum) _coeff[i][j]*b[j]   right multiplication of a vector
	try 
	{
		int myRow = this->GetNumRows();
		int myCol = this->GetNumCols();
		int bRow = b->GetNumRows();
		// rank check
		if (myCol != bRow)
		{
			throw "pdMatrix::Mult(pdVector*): matrix and vector ranks do not match!\n";
		}
		// else continue
		pdVector* c = new pdVector(myRow);
		for (int i=0; i!=myRow; ++i) 
		{
			double sum = 0.0;
			for (int j=0; j!=bRow; ++j)
			{
				sum += this->GetCoeff(i, j) * b->GetCoeff(j);
			}
			c->SetCoeff(i, sum); 
		}
		return c;
	}
	catch (char* str)
	{
		cerr << str;
		exit(0);
	}
}

void pdMatrix::Mult(const double b)
{
	// matrix multiply by a constant b
	int myRow = this->GetNumRows();
	int myCol = this->GetNumCols();
	for (int i=0; i!=myRow; ++i)
	{
		for (int j=0; j!=myCol; ++j)
		{
			this->SetCoeff(i, j, b * this->GetCoeff(i, j));
		}
	}
}

void pdMatrix::Trans(pdMatrix *b)
{
	 //transpose itself and store the tranposed matrix in b 
	try 
	{
		int myRow = this->GetNumRows();
		int myCol = this->GetNumCols();
		int bRow = b->GetNumRows();
		int bCol = b->GetNumCols();
		// rank check
		if (myCol != bRow)
		{
			throw "pdMatrix::Trans(pdMatrix*): matrix ranks do not match!\n";
		}
		if (myRow != bCol)
		{
			throw "pdMatrix::Trans(pdMatrix*): matrix ranks do not match!\n";
		}
		// else continue
		for (int i=0; i!=myCol; ++i)
		{
			for (int j=0; j!=myRow; ++j)
			{
				b->SetCoeff(i, j, this->GetCoeff(j, i));
			}
		}
	}
	catch (char* str)
	{
		cerr << str;
		exit(0);
	}
}
		
void pdMatrix::Gauss(pdVector *b, pdVector *x)
{
	//Here is the routine for Gauss-Jordan elimination with full pivoting:
	//Never use Gauss-Jordan elimination without pivoting!
	//It solves equation [a]{x}={b} where matrix [a] is n by n
	//Once done, the solution {x}=inverse([a])*{b} is stored in {b} vector
	//and the inverse([a]) in [a] matrix

	#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

	int n, i, icol, irow, j, k, l, ll;
	double big, dum, pivinv, temp, hold;
	//	dynamially create 3 integer arrays. The integer arrays 
	// ipiv, indxr, and indxc are used for bookkeeping on the pivoting. 

	n=this->GetNumRows();
	int *indxc=new int[n]; 
	int *indxr=new int[n];
	int *ipiv=new int[n];
	
	for (j=0;j<n;j++)
		ipiv[j]=0;
	for (i=0;i<n;i++) 
	{ 
		// This is the main loop over the columns to be reduced. 
		big=0.0;
		for (j=0;j<n;j++) 
		{
			// This is the outer loop of the search for a pivot element. 
			if (ipiv[j] != 1)
			{	
				for (k=0;k<n;k++) 
				{	
					if (ipiv[k] == 0) 
					//{	if (fabs(a[j][k]) >= big) 
						//{	big=fabs(a[j][k]);
					{
						if (fabs(this->GetCoeff(j, k)) >= big) 
						{
							big=fabs(this->GetCoeff(j, k));
							irow=j;
							icol=k;
						}
					}
					else if (ipiv[k] > 1) 
					{
						cerr << "pdMatrix::Gauss: Singular Matrix-1" << endl;
						exit(0);
					}
				}
			}
		}
		++(ipiv[icol]);

	/* We now have the pivot element, so we interchange rows, if needed, to put the pivot
	element on the diagonal. The columns are not physically interchanged, only relabeled:
	indxc[i], the column of the ith pivot element, is the ith column that is reduced, while
	indxr[i] is the row in which that pivot element was originally located. If indxr[i]
	6 = indxc[i] there is an implied column interchange. With this form of bookkeeping, the
	solution b's will end up in the correct order, and the inverse matrix will be scrambled
	by columns. */
		if (irow != icol)
		//{	for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l]);
		//	SWAP(b[irow],b[icol]);
		{	
			for (l=0;l<n;l++) 
			{
				double rl = this->GetCoeff(irow, l);
				double cl = this->GetCoeff(icol, l);
				SWAP(rl, cl);
				this->SetCoeff(irow, l, rl);
				this->SetCoeff(icol, l, cl);
			}
			hold=b->GetCoeff(irow);
			b->SetCoeff(irow,b->GetCoeff(icol));
			b->SetCoeff(icol,hold);
		}

	// We are now ready to divide the pivot row by the pivot element, located at irow and icol.	
		indxr[i]=irow; 
		indxc[i]=icol;
		//if (a[icol][icol] == 0.0) 
		if (this->GetCoeff(icol, icol) == 0.0) 
		{
			cerr << "pdMatrix::Gauss: Singular Matrix-2" << endl;
			exit(0);
		}
		//pivinv=1.0/a[icol][icol];
		//a[icol][icol]=1.0;
		//for (l=0;l<n;l++) {a[icol][l] *= pivinv;}
		//b[icol] *= pivinv;
		pivinv=1.0/this->GetCoeff(icol, icol);
		this->SetCoeff(icol, icol, 1.0);
		for (l=0;l<n;l++) 
		{
			this->SetCoeff(icol, l, this->GetCoeff(icol, l)*pivinv);
		}
		b->SetCoeff(icol, b->GetCoeff(icol)*pivinv);
		for (ll=0;ll<n;ll++) //Next, we reduce the rows... except for the pivot one, of course.
		{	
			if (ll != icol) 
			//{	dum=a[ll][icol];
			//	a[ll][icol]=0.0;
			//	for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
				//b[ll] -= b[icol]*dum;
			{	
				dum=this->GetCoeff(ll, icol);
				this->SetCoeff(ll, icol, 0.0);
				for (l=0;l<n;l++) 
					this->SetCoeff(ll, l, this->GetCoeff(ll, l) - this->GetCoeff(icol, l)*dum);
				b->AddCoeff(ll, -dum*b->GetCoeff(icol));
			}
		}
	}

	/* This is the end of the main loop over columns of the reduction. It only remains to unscram-
	ble the solution in view of the column interchanges. We do this by interchanging pairs of
	columns in the reverse order that the permutation was built up. */
	for (l=n-1;l>=0;l--)
	{	
		if (indxr[l] != indxc[l])
		{	
			for (k=0;k<n;k++)
			//{	SWAP(a[k][indxr[l]],a[k][indxc[l]]);
			{	
				double kr = this->GetCoeff(k, indxr[l]);
				double kc = this->GetCoeff(k, indxc[l]);
				SWAP(kr, kc);
				this->SetCoeff(k, indxr[l], kr);
				this->SetCoeff(k, indxc[l], kc);
			}
		}
	} // And we are done.
	delete [] indxc; //delete dynamically created arrays
	delete [] indxr;
	delete [] ipiv;

	for(i=0; i<n; i++){	//copy the solution vector to x object
		x->SetCoeff(i, b->GetCoeff(i));
	}
}

void pdMatrix::Zero() 
{
	// set all matrix components to zero
	for (vector<pdVector*>::iterator itor=_coeff.begin(); itor!=_coeff.end(); ++itor)
	{
		(*itor)->Zero();
	}
}

double pdMatrix::RowSumNorm()
{
    // return the uniform-matrix norm (row-sum norm)
	int myRow = this->GetNumRows();
	int myCol = this->GetNumCols();
    double maxSum = 0.0;
    for (int i = 0; i != myRow; ++i)
    {
        double rowSum = 0.0;
        for (int j = 0; j != myCol; ++j)
        {
            rowSum += abs(this->GetCoeff(i, j));
        }
        maxSum = max(maxSum, rowSum);
    }
    return maxSum;
}

double pdMatrix::FBnorm()
{
    // return the Frobenius norm of the matrix
	int myRow = this->GetNumRows();
	int myCol = this->GetNumCols();
    double sum = 0.0;
    for (int i = 0; i != myRow; ++i)
	{
        for (int j = 0; j != myCol; ++j)
		{
            sum += pow(this->GetCoeff(i, j), 2.0);
		}
	}
    return sqrt(sum);
}

void pdMatrix::Print(std::ostream &os) const
{
	int myRow = this->GetNumRows();
	int myCol = this->GetNumCols();
	for (int i = 0; i != myRow; ++i)
	{
        for (int j = 0; j != myCol; ++j)
		{
			os << this->GetCoeff(i, j) << "\t";
		}
		os << endl;
	}
}

ostream& operator << (ostream& os, const pdMatrix& v)
{
	v.Print(os);
	return os;
}
