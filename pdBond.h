/*

The pdBond class implement a peridynamic bond.

*/

#ifndef PDBOND_H
#define PDBOND_H

#include <vector>
#include <iostream>
using std::ostream;
using std::vector;

class pdSpacePlane;
class pdNode;

typedef void (*new_handler)();

class pdBond
{
public:
	// Constuctor
	//pdBond(int, pdNode*, pdNode*);
	pdBond(pdNode*, pdNode*);

	// Destructor
	~pdBond();

	// Function
	void SetBreak(bool);
	//void SetStretch(double);
	void SetMinDistance(double);
	void SetMaxDistance(double);
	//int GetID() const;
	pdNode* GetNodeI() const;
	pdNode* GetNodeJ() const;
	//double GetStretch() const;
	double GetMinDistance() const;
	double GetMaxDistance() const;
	bool IsBreak() const;
	// force-through-plane related function (10-13-2008)
	//void AddPlane(pdSpacePlane*);
	//int GetNumberPlane();
	//pdSpacePlane* GetPlane(int);
	
	// Overload operator
	friend ostream& operator << (ostream&, const pdBond&);

	//// customered operator new and delete
	//static new_handler set_new_handler(new_handler p);
	//static void* operator new(size_t size);
	//void outOfMemory();
	//static void operator delete(void* rawMemory, size_t size);

private:
	// Data
	//int _id;	// bond id
	pdNode* _node_i;	// pointers to nodes on two ends. (the node id of node_i is always smaller than node_j)
	pdNode* _node_j;
	//double _stretch;   //bond stretch
	bool _break;       // flag: true if bond between i and j is break (default=false)
	double _minDistance, _maxDistance; // minimum and maximum distance from one node to the other node cubic
	// force-through-plane related data (10-13-2008)
	//vector<pdSpacePlane*> _planes; // store the pointers of planes this bond passes through
	//int _influence;
	//double _force_state_1, _force_state_2, _force_state_3;
	//bool _break_state;
	//static new_handler currentHandler;

	// Never used constructor and copy control
	pdBond();
	pdBond(const pdBond&);
	pdBond& operator=(const pdBond&);

	// Function
	void Print(ostream&) const;
};

#endif
