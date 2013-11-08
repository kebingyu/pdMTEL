/* 

   This class implments an integration point for calculation of the peridynamic force integration over volume.
   It has global coordinates and displacement associated with it. 

*/

#ifndef PDINTGPOINT_H
#define PDINTGPOINT_H

class pdIntgPoint
{
public:
	// Constructor
	pdIntgPoint();

	// Destructor
	virtual ~pdIntgPoint();

	// Function
	void SetX1(const double);
	void SetX2(double);
	void SetX3(double);
	void SetU1(const double);
	void SetU2(const double);
	void SetU3(const double);
	double GetX1() const;
	double GetX2() const;
	double GetX3() const;
	double GetU1() const;
	double GetU2() const;
	double GetU3() const;

private:
	// Data
	double _x1, _x2, _x3;
	double _u1, _u2, _u3;

	// Never used copy control
	pdIntgPoint(const pdIntgPoint&);
	pdIntgPoint& operator=(const pdIntgPoint&);
};

#endif