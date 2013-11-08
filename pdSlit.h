/*

This class stores data for slit. Slit is a 2D subclass of pdSpace used for simulation the crack.
Only rectangular slit is available right now. The rectangle slit has to be normal to one of the three coordinate axises. 

keyword format in the input file:
	slit_(id)
	1 (x1_lo) (x1_hi) (x2_lo) (x2_hi) (x3_lo) (x3_hi)

note: id = slit id, number starts from 1 

*/

#pragma once
#include "pdSpace.h"

class pdSlit : public pdSpace
{
public:
	// Consturctor
	pdSlit(int id, int type);
	
	// Destructor
	virtual ~pdSlit();
	
	// Function
	void SetNorm(int);
	void SetHeadID(int);
	void SetTailID(int);
	void SetExt_head(double);
	void SetExt_tail(double);
	int GetHeadID() const;
	int GetTailID() const;
	int GetNorm() const;
	double GetExt_head() const;
	double GetExt_tail() const;

private:
	// Data
	int _norm;   // slit normal direction: i for xi axis, i=1, 2, 3
	int _hea_id, _tail_id; // node id closest to slit head and tail (depend on the normal direction)
	double _ext_head, _ext_tail; // crack extension of slit head and tail (depend on the normal direction)
	
	// Never used constructor and copy control
	pdSlit();
	pdSlit(const pdSlit&);
	pdSlit& operator=(const pdSlit&);

	// Function
	void Print(ostream&) const;
};