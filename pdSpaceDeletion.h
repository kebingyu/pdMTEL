/*

This class implements a space which is void.

deletion space keyword format in the input file:
	deletion_space_(id)
	(geometry_type) (x1_lo) (x1_hi) (x2_lo) (x2_hi) (x3_lo) (x3_hi)
	
note: 
id = deletion space id, number starts from 1 

*/

#ifndef PDSPACEDELETION_H
#define PDSPACEDELETION_H

#include "pdSpace.h"

class pdDeletion;

class pdSpaceDeletion : public pdSpace
{
public:
	// Consturctor
	pdSpaceDeletion(int, int);

	// Destructor
	~pdSpaceDeletion();

	// Function

private:
	// Data

	// Never used constructor and copy control
	pdSpaceDeletion();
	pdSpaceDeletion(const pdSpaceDeletion&);
	pdSpaceDeletion& operator=(const pdSpaceDeletion&);

	// Function
	void Print(ostream&) const;
};

#endif