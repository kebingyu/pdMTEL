#include "pdNode.h"
#include "math.h"
#include <iostream>
using namespace std;

pdNode::pdNode() {}

pdNode::pdNode(const pdNode& node) {}

pdNode::pdNode(int id) : _id(id) {}

pdNode::~pdNode()
{
	_family.clear();
	_cbList.clear();
	this->ClearGaussianPointList();
	for (vector<pdVector*>::iterator itor=_mlsCoeffList.begin();itor!=_mlsCoeffList.end();++itor)
		delete (*itor);
	_mlsCoeffList.clear();
}

int pdNode::GetID() const
{
	return _id;
}

void pdNode::SetX1(double x)
{
	_x1 = x;
}

void pdNode::SetX2(double x)
{
	_x2 = x;
}
void pdNode::SetX3(double x)
{
	_x3 = x;
}

void pdNode::SetU1(double u)
{
	_u1 = u;
}

void pdNode::SetU2(double u)
{
	_u2 = u;
}
void pdNode::SetU3(double u)
{
	_u3 = u;
}

void pdNode::SetV1(double v)
{
	_v1 = v;
}

void pdNode::SetV2(double v)
{
	_v2 = v;
}

void pdNode::SetV3(double v)
{
	_v3 = v;
}

void pdNode::SetV1Buffer(double v)
{
	_v1Buff = v;
}

void pdNode::SetV2Buffer(double v)
{
	_v2Buff = v;
}

void pdNode::SetV3Buffer(double v)
{
	_v3Buff = v;
}

void pdNode::SetMaterial(pdMaterial* mat)
{
	_mat = mat;
}

void pdNode::SetBC(pdBdryCondition* bc)
{
	_bc = bc;
}

void pdNode::SetNumDmgBond(int bond)
{
	_dmgBond = bond;
}

void pdNode::SetNumYldBond(int bond)
{
	_yldBond = bond;
}

void pdNode::SetDissiEnergy(double d)
{
	_dissiEnergy = d;
}

void pdNode::AddDissiEnergy(double d)
{
	_dissiEnergy += d;
}

double pdNode::GetDissiEnergy() const
{
	return _dissiEnergy;
}

void pdNode::SetElastEnergyDensity(double e)
{
	_elasticEnergyDensity = e;
}

void pdNode::AddElastEnergyDensity(double e)
{
	_elasticEnergyDensity += e;
}

double pdNode::GetElastEnergyDensity()  const
{
	return _elasticEnergyDensity;
}

void pdNode::SetKinetEnergy(double k)
{
	_kineticEnergy = k;
}

void pdNode::AddKinetEnergy(double k)
{
	_kineticEnergy += k;
}

double pdNode::GetKinetEnergy()  const
{
	return _kineticEnergy;
}

void pdNode::SetNodeVolume(double volnod)
{
	_volume = volnod;
}

double pdNode::GetX1() const
{
	return _x1;
}

double pdNode::GetX2() const
{
	return _x2;
}

double pdNode::GetX3() const
{
	return _x3;
}

double pdNode::GetU1() const
{
	return _u1;
}

double pdNode::GetU2() const
{
	return _u2;
}

double pdNode::GetU3() const
{
	return _u3;
}

double pdNode::GetV1() const
{
	return _v1;
}

double pdNode::GetV2() const
{
	return _v2;
}

double pdNode::GetV3() const
{
	return _v3;
}

double pdNode::GetV1Buffer() const
{
	return _v1Buff;
}

double pdNode::GetV2Buffer() const
{
	return _v2Buff;
}

double pdNode::GetV3Buffer() const
{
	return _v3Buff;
}

pdMaterial* pdNode::GetMaterial() const
{
	return _mat;
}

pdBdryCondition* pdNode::GetBC() const
{
	return _bc;
}

double pdNode::GetDmgFraction()  const
{
	int numFam = this->GetNumFamily();
	if (numFam != 0)
		return _dmgBond / numFam;
	else
		return 0.0;
}

void pdNode::AddNumDmgBond(int b)
{
	_dmgBond += b;
}

void pdNode::AddNumYldBond(int b)
{
	_yldBond += b;
}

int pdNode::GetNumDmgBond() const
{
	return _dmgBond;
}

int pdNode::GetNumYldBond() const
{
	return _yldBond;
}

double pdNode::GetYldFraction()  const
{
	int numFam = this->GetNumFamily();
	if (numFam != 0)
		return _yldBond / numFam;
	else
		return 0.0;
}

double pdNode::GetNodeVolume()  const
{
	return _volume;
}

void pdNode::AddV1Buffer(double v)
{
	_v1Buff += v;
}

void pdNode::AddV2Buffer(double v)
{
	_v2Buff += v;
}

void pdNode::AddV3Buffer(double v)
{
	_v3Buff += v;
}

void pdNode::AddFamilyBond(pdBond* bd)
{
	_family.push_back(bd);
}

pdBond* pdNode::GetFamilyBond(int n) const
{
	return _family[n];
}

pdNode* pdNode::GetFamilyNode(int n) const
{
	pdBond* bd = GetFamilyBond(n);
	if ( bd->GetNodeI() == this)
	{
		return bd->GetNodeJ();
	}
	else
	{
		return bd->GetNodeI();
	}
}

int pdNode::GetNumFamily() const
{
	return int(_family.size());
}

void pdNode::SetSupportDomainSize(double rw)
{
	_rw = rw;
}

double pdNode::GetSupportDomainSize() const
{
	return _rw;
}

void pdNode::AddContribNode(pdNode* n)
{
	_cbList.push_back(n);
}

void pdNode::ClearContribNodeList()
{
	_cbList.clear();
}

void pdNode::GetContribNodeList(vector<pdNode*>& target)
{
	// copy this node's contributing node list to target container
	target.resize(_cbList.size());
	copy(_cbList.begin(), _cbList.end(), target.begin());
}

pdNode* pdNode::GetContribNode(const int idx) const
{
	return _cbList[idx];
}

int pdNode::GetContribNodeListSize() const
{
	return (int)_cbList.size();
}

void pdNode::AddGaussianPoint(pdGauss* g)
{
	_gaussPointList.push_back(g);
}

void pdNode::GetGaussianPointList(vector<pdGauss*>& target)
{
	// copy this node's Gaussian point list to target container
	target.resize(_gaussPointList.size());
	copy(_gaussPointList.begin(), _gaussPointList.end(), target.begin());
}

int pdNode::GetGaussianPointListSize() const
{
	return _gaussPointList.size();
}

pdGauss* pdNode::GetGaussianPoint(const int idx)
{
	return _gaussPointList[idx];
}

void pdNode::ClearGaussianPointList()
{
	for (vector<pdGauss*>::iterator itor=_gaussPointList.begin();itor!=_gaussPointList.end();++itor)
		delete *itor;
	_gaussPointList.clear();
}

void pdNode::SetBdryEffFactor(double f)
{
	_bdryEffFactor = f;
}

double pdNode::GetBdryEffFactor() const
{
	return _bdryEffFactor;
}

void pdNode::AddMLSCoeff(pdVector* mlsList)
{
	_mlsCoeffList.push_back(mlsList);
}

pdVector* pdNode::GetMLSCoeff(const int idx)
{
	try
	{
		if (idx > int(_mlsCoeffList.size()))
		{
			throw "pdNode::GetMLSCoeff: index out of vector scope!";
		}
		else
		{
			return _mlsCoeffList[idx];
		}
	}
	catch (char* str)
	{
		cerr << str;
		exit(0);
	}

}

int pdNode::GetNumTrapPoint() const
{
	return int(_mlsCoeffList.size());
}

void pdNode::Print(ostream &os) const
{
	os << _id 
		<< "\t" << _x1 << "\t" << _x2 << "\t" << _x3
		<< "\t" << _u1 << "\t" << _u2 << "\t" << _u3 
		<< "\t" << _v1 << "\t" << _v2 << "\t" << _v3 << "\t" << _mat->GetID();
	if (_bc!=0)
	{
		os << "\t" << _bc->GetID();
	}
	else
	{
		os << "\t" << -1; // indicate no boundary condition is applied to this node
	}
	os << "\t" << _bdryEffFactor << "\t" << GetDmgFraction();
}

//void pdNode::Print(ostream &os) const
//{
//	os << "Node (id = " << _id << "):" << endl
//		<< "coord: " << _x1 << "\t" << _x2 << "\t" << _x3 << endl
//		<< "displ: " << _u1 << "\t" << _u2 << "\t" << _u3 << endl
//		<< "veloc: " << _v1 << "\t" << _v2 << "\t" << _v3 << endl
//		<< "volume= " << _volume << "\t" << "material id= " << _mat->GetID();
//	if (_bc!=0)
//	{
//		os << "\t" << "boundary condition id= " << _bc->GetID() << endl;
//	}
//	else
//	{
//		os << "\t" << "boundary condition id= NULL " << endl;
//	}
//	os << "elastic energy density= " << _elasticEnergyDensity << "\t" << "kinetic energy= " << _kineticEnergy << endl
//		<< "damage fraction= " << GetDmgFraction() << "\t" << "yield fraction= " << GetYldFraction() 
//		<< "\t" << "bdryEff= " << _bdryEffFactor << endl;
//}

ostream& operator << (ostream& os, const pdNode& n)
{
	n.Print(os);
	return os;
}