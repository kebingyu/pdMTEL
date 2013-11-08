/*

This class implements a peridynamic node.

mid-point method is used in the time integration so the nodal displcaement is calculated at full step, 
and the nodal velocity is calculated at half step.			
	v_(n+0.5) = v_(n-0.5) + dt*f(u_n, t_n) or vmid_(n+1) = vmid_n + dt*f(u_n, t_n)
	u_(n+1) = u_n + v_(n+0.5)*dt or u_(n+1) = u_n + vmid_(n+1)*dt

*/

#ifndef PDNODE_H
#define PDNODE_H

#include <vector>
#include <iostream>
#include <algorithm>
#include "pdBond.h"
#include "pdMaterial.h"
#include "pdBdryCondition.h"
#include "pdGauss.h"
#include "pdVector.h"
using std::ostream;
using std::endl;
using std::vector;

class pdNode
{
public:
	// Constructor
	pdNode(int id);

	// Destructor
	~pdNode();

	// Functions
	void SetX1(double);
	void SetX2(double);
	void SetX3(double);	
	void SetU1(double);
	void SetU2(double);
	void SetU3(double);	
	void SetV1(double);
	void SetV2(double);
	void SetV3(double);
	void SetV1Buffer(double);
	void SetV2Buffer(double);
	void SetV3Buffer(double);
	void SetMaterial(pdMaterial*);
	void SetBC(pdBdryCondition*);
	void SetNumDmgBond(int);
	void SetNumYldBond(int);
	void SetElastEnergyDensity(double);	
	void SetKinetEnergy(double);	
	void SetDissiEnergy(double);	
	void SetNodeVolume(double);
	void SetSupportDomainSize(double);
	void AddFamilyBond(pdBond*);
	void AddNumDmgBond(int);
	void AddNumYldBond(int);
	void AddElastEnergyDensity(double);
	void AddKinetEnergy(double);
	void AddDissiEnergy(double);
	void AddV1Buffer(double);
	void AddV2Buffer(double);
	void AddV3Buffer(double);
	int GetID() const;
	pdMaterial* GetMaterial() const;
	pdBdryCondition* GetBC() const;
	int GetNumDmgBond() const;	// return the current total number of broken bonds
	int GetNumYldBond() const; // return the current total number of yielded bonds
	pdBond* GetFamilyBond(int n) const; // return the nth family bond
	pdNode* GetFamilyNode(int n) const; // return node_j on the other end of the nth family bond
	int GetNumFamily() const; // return the total number of family at initialization
	double GetX1() const;  
	double GetX2() const;	
	double GetX3() const;
	double GetU1() const;
	double GetU2() const;
	double GetU3() const;
	double GetV1() const;
	double GetV2() const;
	double GetV3() const;
	double GetV1Buffer() const;
	double GetV2Buffer() const;
	double GetV3Buffer() const;
	double GetDmgFraction() const;  //return total damage
	double GetElastEnergyDensity() const; //return local strain energy density
	double GetKinetEnergy() const; //return kinetic energy
	double GetDissiEnergy() const;
	double GetYldFraction() const; //return yield fraction
	double GetNodeVolume() const;
	double GetSupportDomainSize() const;
	void AddContribNode(pdNode*);
	void ClearContribNodeList();
	void GetContribNodeList(vector<pdNode*>&);
	pdNode* GetContribNode(const int) const;
	int GetContribNodeListSize() const;
	void AddGaussianPoint(pdGauss*);
	void GetGaussianPointList(vector<pdGauss*>&);
	int GetGaussianPointListSize() const;
	pdGauss* GetGaussianPoint(const int);
	void ClearGaussianPointList();
	void SetBdryEffFactor(double);
	double GetBdryEffFactor() const;
	void AddMLSCoeff(pdVector*);
	pdVector* GetMLSCoeff(const int);
	int GetNumTrapPoint() const;

	// Overload operator
	friend ostream& operator << (ostream&, const pdNode&);

private:
	// Data	
	int _id;	         // node id, start from 0
	int _dmgBond;       // total number of damaged bond
	int _yldBond;       // total number od yielded bond
	double _x1, _x2, _x3; // 3D coordinates
	double _u1, _u2, _u3; // displacement at current full time step
	double _v1, _v2, _v3; // velocity at current half time step
	double _v1Buff, _v2Buff, _v3Buff; // updated velocity at each time step
	double _elasticEnergyDensity;        // local strain energy density (per unit volume)
	double _kineticEnergy;       // node kinetic energy
	double _dissiEnergy;		// dissipated energy density	
	double _volume;	// node volume
	double _rw; // support domain size. right now it is equal to all of the supporting nodes
	double _bdryEffFactor; // boundary effect compensation factor
	pdMaterial* _mat; // pointer to material associated with this node, 0 for not in any material space(fatal error) 
	pdBdryCondition* _bc; // pointer to boundary condition associated with this node, 0 for no initial boundary condition
	vector<pdBond*> _family; // container for bonds connecting to this node	(family bonds)	
    vector<pdNode*> _cbList; // container of contributing nodes within the support domain (use to get the approximated displacement at the integration point)
	vector<pdGauss*> _gaussPointList; // container of fixed Gaussian integration points in this node's cell
	vector<pdVector*> _mlsCoeffList; // container of MLS coefficients for the trapezoidal points in this node

	// Never used constructors and copy control
	pdNode();
	pdNode(const pdNode&);
	pdNode& operator=(const pdNode&);

	// Function
	void Print(ostream&) const;
};

#endif
