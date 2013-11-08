#include "pdDataManager.h"
//#include "random.h"
using namespace std;

pdDataManager::pdDataManager() {}

pdDataManager::~pdDataManager()
{
	_mats.clear();
	_BCs.clear();
	ClearVectorList(_planes);
	ClearVectorList(_nodes);
	ClearVectorList(_bonds);
	ClearVectorList(_myBonds);
	ClearVectorList(_delSpaces);
	ClearVectorList(_matSpaces);
	ClearVectorList(_bdrySpaces);
}

template <typename T>
void pdDataManager::ClearVectorList(std::vector<T*>& vec)
{
	for (vector<T*>::iterator itor=vec.begin();itor!=vec.end();++itor)
	{
		if (*itor!=0)
		{
			delete *itor;
		}
	}
	vec.clear();
}

void pdDataManager::Debug()
{
	string ofpath = _workDir + "pds_debug.dat";
	ofstream fout(ofpath.c_str(), ios::out);

	//cout << " grid spacing = ";
	//cin >> _gridSpacing;
	//_tolerance = 1e-4*_gridSpacing;
	//cout << " horizon = ";
	//cin >> _horizon;	
	////cout << " Gaussian point = ";
	////cin >> _intgPoint;	
	//// initial node list
	//cout << " grid dimension = ";
	//int dimen;
	//cin >> dimen;
	//_gridDimensionX1 = dimen;
	//_gridDimensionX2 = dimen;
	//_gridDimensionX3 = dimen;
	//_gridCenterX1 = 0;
	//_gridCenterX2 = 0;
	//_gridCenterX3 = 0;
	//int id = 0;
	//for (int n=1;n!=_gridDimensionX3+1;++n) // z direction
	//{
	//	for (int m=1;m!=_gridDimensionX2+1;++m) // y direction
	//	{
	//		for (int l=1;l!=_gridDimensionX1+1;++l) // x direction
	//		{
	//			double x1i = (l-(_gridDimensionX1+1)*0.5)*_gridSpacing + _gridCenterX1;
	//			double x2i = (m-(_gridDimensionX2+1)*0.5)*_gridSpacing + _gridCenterX2;
	//			double x3i = (n-(_gridDimensionX3+1)*0.5)*_gridSpacing + _gridCenterX3;
	//			pdNode* node_i = new pdNode(id);
	//			node_i->SetX1(x1i);
	//			node_i->SetX2(x2i);
	//			node_i->SetX3(x3i);
	//			node_i->SetU1(0.0);
	//			node_i->SetU2(0.0);
	//			node_i->SetU3(0.0);
	//			_nodes.push_back(node_i);
	//			++id;
	//		}
	//	}
	//}
	//_numNodes = int(_nodes.size());
	
	//// find the center node and only set up bonds connecting to this node
	//pdNode* node_i;
	//for (vector<pdNode*>::const_iterator itor_i=_nodes.begin();itor_i!=_nodes.end();++itor_i)
	//{
	//	pdNode* node = *itor_i;
	//	if (node->GetX1() == 0.5*_gridSpacing && node->GetX2() == 0.5*_gridSpacing 
	//		&& node->GetX3() == 0.5*_gridSpacing) // find the center node
	//	{
	//		node_i = node;
	//		double x1i = node_i->GetX1();
	//		double x2i = node_i->GetX2();
	//		double x3i = node_i->GetX3();
	//		// loop over all nodes to set up bonds just for the center node
	//		for (vector<pdNode*>::const_iterator itor_j=_nodes.begin();itor_j!=_nodes.end();++itor_j)
	//		{
	//			pdNode* node_j = *itor_j;
	//			if (node_j != node_i)
	//			{
	//				double x1j = node_j->GetX1();
	//				double x2j = node_j->GetX2();
	//				double x3j = node_j->GetX3();
	//				double dist_ij = sqrt(pow(x1j - x1i, 2.0) + pow(x2j - x2i, 2.0) + pow(x3j - x3i, 2.0));
	//				if (_horizon < dist_ij - 1.5*_gridSpacing)
	//				{
	//					// a rough judge to skip calculation for two far away nodes
	//					continue;
	//				}
	//				else
	//				{
	//					// find the shortest and longest distance from node i to the cubic volume of node j
	//					double dist_min, dist_max;
	//					CalcDistanceToCube(node_i, node_j, dist_min, dist_max);
	//					if (dist_min < _horizon + _tolerance)
	//					{
	//						// add new bond to the container
	//						pdBond* bond_ij = new pdBond(node_i, node_j); // here _numBonds serves as bond id
	//						// save the min and max distance, set the initial bond breakage status to unbroken
	//						bond_ij->SetMinDistance(dist_min);
	//						bond_ij->SetMaxDistance(dist_max);
	//						bond_ij->SetBreak(false);
	//						_bonds.push_back(bond_ij);
	//						// find the shortest bond length
	//						_minBondLength = min(_minBondLength, dist_ij);
	//						// store ptr to this bond to node i and j family container
	//						node_i->AddFamilyBond(bond_ij);
	//						node_j->AddFamilyBond(bond_ij);
	//					}
	//					//// horizon as a cube (right now the length of the horizon cube always equals n*_gridSpacing)
	//					//// check if the family node cell is fully inside the horizon cube
	//					//if (x1j+0.5*_gridSpacing <= x1i+_horizon && x1j-0.5*_gridSpacing >= x1i-_horizon
	//					//	&& x2j+0.5*_gridSpacing <= x2i+_horizon && x2j-0.5*_gridSpacing >= x2i-_horizon
	//					//	&& x3j+0.5*_gridSpacing <= x3i+_horizon && x3j-0.5*_gridSpacing >= x3i-_horizon)
	//					//{
	//					//	// add new bond to the container
	//					//	pdBond* bond_ij = new pdBond(id, node_i, node_j); // here _numBonds serves as bond id
	//					//	// save the min and max distance, set the initial bond breakage status to unbroken
	//					//	bond_ij->SetBreak(false);
	//					//	_bonds.push_back(bond_ij);
	//					//	// find the shortest bond length
	//					//	_minBondLength = min(_minBondLength, dist_ij);
	//					//	// store ptr to this bond to node i and j family container
	//					//	node_i->AddFamilyBond(bond_ij);
	//					//	node_j->AddFamilyBond(bond_ij);
	//					//	++id;
	//					//}
	//				}
	//			}
	//		}
	//		break;
	//	}
	//}

	//// find the center node (bonds are already set up in IniNode subroutine)
	//pdNode* node_i;
	//for (vector<pdNode*>::const_iterator itor_i=_nodes.begin();itor_i!=_nodes.end();++itor_i)
	//{
	//	pdNode* node = *itor_i;
	//	if (node->GetX1() == 0.5*_gridSpacing && node->GetX2() == 0.5*_gridSpacing
	//		&& node->GetX3() == 0.5*_gridSpacing) // find the center node
	//	{
	//		node_i = node;
	//		break;
	//	}
	//}

	//// volume of the horizon sphere calculation
	//fout << "grid spacing = " << _gridSpacing << "\n"
	//	<< "horizon = " << _horizon << "\n\n";
	////fout << " Gaussian_point vol_code vol_accurate error(%)" << endl;
	////int GPoint [7] = {1, 2, 3, 4, 5, 6, 8};
	////for (int g=0;g!=7;++g)
	////{
	//	//_intgPoint = GPoint[g];
	//	//SetGaussianPointsValue();
	//	//SetGaussianPoints();

	//	// calculate intersection volume for each family node of the center node
	//	double intgVol = 0.0;
	//	int totalIntgPoint = 0;
	//	for (int f = 0; f != node_i->GetNumFamily(); ++f)
	//	{
	//		pdBond* bond = node_i->GetFamilyBond(f);
	//		pdNode* node_j = node_i->GetFamilyNode(f);
	//		int idj = node_j->GetID();
	//		double vol;
	//		if (_cciOn)// use CCI method to count the family nodes and calculate the intersection volume
	//		{
	//			double xi1 = node_j->GetX1() - node_i->GetX1();
	//			double xi2 = node_j->GetX2() - node_i->GetX2();
	//			double xi3 = node_j->GetX3() - node_i->GetX3();
	//			double r = sqrt(xi1*xi1 + xi2*xi2 + xi3*xi3);
	//			double beta = this->CalcBetaFactor(r, _horizon, _gridSpacing);
	//			vol = pow(_gridSpacing, 3.0)*beta;
	//			if (beta>0.0)
	//				++totalIntgPoint;

	//			//// refined grid size
	//			//// the problem is you have to refine every family node cell at different level ??
	//			//vol = 0.0;
	//			//int x1 = 8;
	//			//int x2 = 2;
	//			//int x3 = 2;
	//			//for (int n=1;n!=x3+1;++n) // z direction
	//			//{
	//			//	for (int m=1;m!=x2+1;++m) // y direction
	//			//	{
	//			//		for (int l=1;l!=x1+1;++l) // x direction
	//			//		{
	//			//			double x1j = (l - (x1 + 1) * 0.5) * (_gridSpacing/x1) + node_j->GetX1();
	//			//			double x2j = (m - (x2 + 1) * 0.5) * (_gridSpacing/x2) + node_j->GetX2();
	//			//			double x3j = (n - (x3 + 1) * 0.5) * (_gridSpacing/x3) + node_j->GetX3();
	//			//			double xi1 = x1j - node_i->GetX1();
	//			//			double xi2 = x2j - node_i->GetX2();
	//			//			double xi3 = x3j - node_i->GetX3();
	//			//			double r = sqrt(xi1*xi1 + xi2*xi2 + xi3*xi3);
	//			//			double effDX = pow(pow(_gridSpacing, 3)/(x1*x2*x3), 1.0/3.0); //effective refined grid size
	//			//			double beta = this->CalcBetaFactor(r, _horizon, effDX);
	//			//			vol += pow(_gridSpacing, 3.0)*beta/(x1*x2*x3);	
	//			//		}
	//			//	}
	//			//}

	//			// print data
	//			/*fout << idj	<< "\t" << node_j->GetX1() <<"\t" << node_j->GetX2() 
	//				<< "\t" << node_j->GetX3() << "\t" << vol << "\n";*/
	//		}
	//		else
	//		{
	//			pdIntgManager* im = new pdIntgManager(_horizon, _gridSpacing, bond, node_i, node_j);
	//			im->SetIntgPoint(_intgPoint);
	//			im->SetTrapFlag(_aiOn);
	//			im->SetErrorControlFlag(_errorControlOn);
	//			im->SetErrorControlEPS(_eps);
	//			im->SetGaussianFlag(_fgiOn);
	//			im->SetGaussianAbscis(_abscis);
	//			im->SetGaussianWeight(_weights);
	//			im->SetMLSFlag(_mlsOn);
	//			double vTemp[3];
	//			if (_horizon <= bond->GetMinDistance()) // too far away (should never get here)
	//			{
	//				vol = 0.0;
	//			}
	//			else if (_horizon >= bond->GetMaxDistance()) // cubic of node j is fully in the horizon of node i
	//			{					
	//				im->CalcConfigurationType3("UnitValue", vTemp);
	//				vol = vTemp[0];
	//			}
	//			else
	//			{
	//				if (_aiOn)
	//					im->AdaptiveIntegration("UnitValue", vTemp);
	//				else if (_fgiOn)
	//					im->CalcConfigurationType3("UnitValue", vTemp);
	//				vol = vTemp[0];
	//			}
	//			//// print data
	//			//fout << idj	<< "\t" << node_j->GetX1() <<"\t" << node_j->GetX2() << "\t" << node_j->GetX3()
	//			//	<< "\t" << im->GetEndIntgPoint(1) << "\t" << im->GetEndIntgPoint(2) << "\t" << im->GetEndIntgPoint(3) 
	//			//	<< "\t" << vol << "\n";
	//			if (_aiOn)
	//				totalIntgPoint += im->GetEndIntgPoint(1)*im->GetEndIntgPoint(2)*im->GetEndIntgPoint(3);
	//			else // _fgiOn
	//				totalIntgPoint += im->GetEndIntgPoint(0);
	//			delete im;
	//			//// simplified code
	//			//vol = 0.0;
	//			//for (int i=0;i!=_intgPoint;++i)
	//			//	for (int j=0;j!=_intgPoint;++j)
	//			//		for (int k=0;k!=_intgPoint;++k)
	//			//		{
	//			//			// relative distance between Gaussian point and center node
	//			//			double xi1 = node_j->GetX1() + 0.5*_gridSpacing*_abscis[i] - node_i->GetX1();
	//			//			double xi2 = node_j->GetX2() + 0.5*_gridSpacing*_abscis[j] - node_i->GetX2();
	//			//			double xi3 = node_j->GetX3() + 0.5*_gridSpacing*_abscis[k] - node_i->GetX3();
	//			//			double r = sqrt(xi1*xi1 + xi2*xi2 + xi3*xi3);

	//			//			//// (1) discard Gaussian points outside the horizon
	//			//			//if (r <= _horizon)
	//			//			//{
	//			//			//	vol += (pow(_gridSpacing, 3.0)/8.0)*_weights[i]*_weights[j]*_weights[k];
	//			//			//}

	//			//			// (2) CCI beta inside the loop
	//			//			double beta;
	//			//			//double w = min(_weights[i], _weights[j]);
	//			//			//w = min(w, _weights[k]);
	//			//			//double w = max(_weights[i], _weights[j]);
	//			//			//w = max(w, _weights[k]);
	//			//			//double w = (_weights[i]+_weights[j]+_weights[k])/3.0;
	//			//			double w = pow(_weights[i]*_weights[j]*_weights[k], 1.0/3.0);
	//			//			double dx = w*_gridSpacing/2.0; // find the cell size belong to this Gaussian point
	//			//			if (r <= _horizon - 0.5*dx)
	//			//				beta = 1.0;
	//			//			else if (r > _horizon - 0.5*dx && r <= _horizon + 0.5*dx)
	//			//				beta = (_horizon+0.5*dx-r)/dx; // (original CCI method, scale the volume of boundary family nodes)
	//			//			else
	//			//				beta = 0.0;
	//			//			vol += (pow(_gridSpacing, 3.0)/8.0)*_weights[i]*_weights[j]*_weights[k]*beta;

	//			//			//// (3) CCI beta outside the loop
	//			//			//vol += (pow(_gridSpacing, 3.0)/8.0)*_weights[i]*_weights[j]*_weights[k];
	//			//		}
	//		}
	//		//// (3) CCI beta outside the loop
	//		//double beta;
	//		//double xi1 = node_j->GetX1() - node_i->GetX1();
	//		//double xi2 = node_j->GetX2() - node_i->GetX2();
	//		//double xi3 = node_j->GetX3() - node_i->GetX3();
	//		//double r = sqrt(xi1*xi1 + xi2*xi2 + xi3*xi3);
	//		//if (r <= _horizon - 0.5*_gridSpacing)
	//		//	beta = 1.0;
	//		//else if (r > _horizon - 0.5*_gridSpacing && r <= _horizon + 0.5*_gridSpacing)
	//		//	beta = (_horizon+0.5*_gridSpacing-r)/_gridSpacing; // (original CCI method, scale the volume of boundary family nodes)
	//		//else
	//		//	beta = 0.0;
	//		//vol *= beta;
	//		// sum to the total volume
	//		intgVol += vol;
	//		//// print nodal and volume data
	//		//fout << idj
	//		//	<< "\t" << node_j->GetX1() <<"\t" << node_j->GetX2() << "\t" << node_j->GetX3()
	//		//	<< "\t" << bond->GetMinDistance() << "\t" << vol << endl;	
	//	}
	//	intgVol += pow(_gridSpacing, 3.0); // add the volume of the source node
	//	double accVol = 4.0 * acos(-1.0) * pow(_horizon, 3.0) / 3.0; // get the accurate volume of the horizon sphere
	//	//double accVol = pow(_horizon*2, 3.0); // horizon as a cube
	//	if (_cciOn)
	//	{
	//		fout << "\n Results with CCI method on: \n"
	//			<< "total integration points = " << totalIntgPoint << "\n";
	//	}
	//	else if (_aiOn)
	//	{
	//		fout << "\n Results with trapzoidal quadrature on: \n"
	//			<< "error control accuracy = " << _eps << "\n"
	//			<< "starting trapezoidal index = " << _intgPoint << "\n"
	//			<< "total integration points = " << totalIntgPoint << "\n";
	//	}
	//	else // Gaussian integration
	//	{
	//		fout << "\n Results with Gaussian quadrature on: \n"
	//			<< "starting Gaussian points = " << _intgPoint << "\n"
	//			<< "total integration points = " << totalIntgPoint << "\n";
	//	}
	//	
	//	fout << "calculated volume = " << intgVol << "\n" 
	//		<< "accurate volume = " << accVol << "\n" 
	//		<< "error percentage = " << 100*abs(intgVol-accVol)/accVol << "\n";
	////}


	//// find the center source node
	//pdNode* node_s;
	//for (vector<pdNode*>::const_iterator itor_i=_nodes.begin();itor_i!=_nodes.end();++itor_i)
	//{
	//	pdNode* node = *itor_i;
	//	if (node->GetX1() == 0.5 * _gridSpacing && node->GetX2() == 0.5 * _gridSpacing
	//		&& node->GetX3() == 0.5 * _gridSpacing) // find the center node
	//	{
	//		node_s = node;
	//		break;
	//	}		
	//}	
	//UpdateNodeStatus(1.0, 1e-4); // to immediately apply the bc
	////if (!_mlsOn) // uncomment if debug routine does not include the InitialProblem(.) function
	////{
	////	SetGaussianPoints();
	////}
 //   double sigma[9];
	//int totalNumIntg;
 //   CalcArealForceDensity(node_s, sigma, totalNumIntg, fout);
	//fout << "grid spacing = " << _gridSpacing << "\n";
	//fout << "horizon = " << _horizon << "\n\n";
	//fout << "total number of integration points = " << totalNumIntg << "\n\n";
 //   fout << " sigma11 = " << sigma[0] << " sigma12 = " << sigma[1] << " sigma13 = " << sigma[2] << endl;
 //   fout << " sigma21 = " << sigma[3] << " sigma22 = " << sigma[4] << " sigma23 = " << sigma[5] << endl;
	//fout << " sigma31 = " << sigma[6] << " sigma32 = " << sigma[7] << " sigma33 = " << sigma[8] << endl;

	UpdateNodeStatus(1.0, 1e-4); // to immediately apply the bc
	//// output displacements on the trapezoidal points
	//for (vector<pdNode*>::const_iterator itor_i=_nodes.begin();itor_i!=_nodes.end();++itor_i)
	//{
	//	pdNode* node = *itor_i;
	//	fout << node->GetID() << "\t" << node->GetX1() << "\t" << node->GetX2() << "\t"	<< node->GetX3()
	//		<< "\t" << node->GetU1() << "\t" << node->GetU2() << "\t"	<< node->GetU3() << endl;
	//	int cbListSize = node->GetContribNodeListSize();
	//	int numTrap = node->GetNumTrapPoint();
	//	for (int p=0;p!=numTrap;++p)
	//	{
	//		double u1p = 0.0;
	//		double u2p = 0.0;
	//		double u3p = 0.0;			
	//		for (int i=0;i!=cbListSize;++i)
	//		{
	//			u1p += node->GetContribNode(i)->GetU1() * node->GetMLSCoeff(p)->GetCoeff(i);
	//			u2p += node->GetContribNode(i)->GetU2() * node->GetMLSCoeff(p)->GetCoeff(i);
	//			u3p += node->GetContribNode(i)->GetU3() * node->GetMLSCoeff(p)->GetCoeff(i);
	//		}
	//		fout << "  trap point " << p << "\t" << u1p << "\t" << u2p << "\t" << u3p << endl;
	//	}
	//	fout << endl;
	//}
	// find the center source node
	pdNode* node_s;
	for (vector<pdNode*>::const_iterator itor_i=_nodes.begin();itor_i!=_nodes.end();++itor_i)
	{
		pdNode* node = *itor_i;
		if (node->GetX1() == 0.5 * _gridSpacing && node->GetX2() == 0.5 * _gridSpacing
			&& node->GetX3() == 0.5 * _gridSpacing) // find the center node
		{
			node_s = node;
			break;
		}		
	}
	pdNode* node_j = node_s->GetFamilyNode(0);
	pdBond* bond_ij = node_s->GetFamilyBond(0);
	double pdForce[3];
	int numIntg;
	vector<pdNode*> cbList;
	node_j->GetContribNodeList(cbList);
	int ndNumber = cbList.size();
	// create a 2D array to store the displacements of contributing nodes
	double** cbDispl = new double*[ndNumber]; // allocate memory for array of displacement of column
	double* curPtr = new double[ndNumber*3]; // allocate total memory required
	for (int i=0;i!=ndNumber;++i)
	{
		*(cbDispl+i) = curPtr; // point the pointers to the right memory
		curPtr += 3;
	}
	for (int j=0;j!=ndNumber;++j)
	{
		cbDispl[j][0] = cbList[j]->GetU1();
		cbDispl[j][1] = cbList[j]->GetU2();
		cbDispl[j][2] = cbList[j]->GetU3();
	}
	double u[3] = {node_s->GetU1(), node_s->GetU2(), node_s->GetU3()};
	CalcBondForceIntg(bond_ij, u, node_s, cbDispl, node_j, 1.0, pdForce, numIntg);
	fout << pdForce[0] << "\t" << pdForce[1] << "\t" << pdForce[2] << "\t" << numIntg << endl; 

	
	fout.close();
}

void pdDataManager::CalcArealForceDensity(pdNode* node_s, double* sigma, int& totalNumIntg, ofstream& fout)
{
	/*
	
	This function calculates the stress tensor based on the areal force density at the source node "node_s".

	Input:
	  1. node_s: pointer to the source node.

	Output: 
	  1. sigma: array to sotre the nine components of the stress tensor.
      2. totalNumIntg: total number of integration points used in the calculation.

	*/

	double x1s = node_s->GetX1();
	double x2s = node_s->GetX2();
	double x3s = node_s->GetX3();

	// container used to store the ptr to the colinear points in three directions
	vector<pdNode*> node_x1;
	vector<pdNode*> node_x2;
	vector<pdNode*> node_x3;
	int ns = int(_horizon/_gridSpacing);
	// find the colinear points in three directions
	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{
		pdNode* node_i = *itor;
		double x1i = node_i->GetX1();
		double x2i = node_i->GetX2();
		double x3i = node_i->GetX3();
		for(int is=0;is!=ns;++is)
		{
			if ( x1i == x1s - _gridSpacing*(1+is) && x2i == x2s && x3i == x3s )
				node_x1.push_back(node_i);
			if ( x1i == x1s && x2i == x2s - _gridSpacing*(1+is) && x3i == x3s )
				node_x2.push_back(node_i);
			if ( x1i == x1s && x2i == x2s && x3i == x3s - _gridSpacing*(1+is) )
				node_x3.push_back(node_i);
		}
	}

	double ds = _horizon/ns;
	double tol = _tolerance;
	double dx = _gridSpacing;

	// initiate stress tensor
	double sig[9];
	memset(sig, 0, sizeof(sig));
	double f[3];
	totalNumIntg = 0;
	int numIntg;
	for (vector<pdNode*>::const_iterator itor=node_x1.begin();itor!=node_x1.end();++itor)
	{
		pdNode* node_i = *itor;
		CalcForceOfPlusSpace(node_i, node_s, 1, f, numIntg, fout);
		totalNumIntg += numIntg;
		// update stress
		sig[0] += f[0]*ds;
		sig[1] += f[1]*ds;
		sig[2] += f[2]*ds;
	}
	for (vector<pdNode*>::const_iterator itor=node_x2.begin();itor!=node_x2.end();++itor)
	{
		pdNode* node_i = *itor;
		CalcForceOfPlusSpace(node_i, node_s, 2, f, numIntg, fout);
		totalNumIntg += numIntg;
		// update stress
		sig[3] += f[0]*ds;
		sig[4] += f[1]*ds;
		sig[5] += f[2]*ds;
	}
	for (vector<pdNode*>::const_iterator itor=node_x3.begin();itor!=node_x3.end();++itor)
	{
		pdNode* node_i = *itor;
		CalcForceOfPlusSpace(node_i, node_s, 3, f, numIntg, fout);
		totalNumIntg += numIntg;
		// update stress
		sig[6] += f[0]*ds;
		sig[7] += f[1]*ds;
		sig[8] += f[2]*ds;
	}
	for (int i=0;i!=9;++i)
	{
		sigma[i] = sig[i];
	}
	node_x1.clear();
	node_x2.clear();
	node_x3.clear();
}

void pdDataManager::CalcForceOfPlusSpace(pdNode* node_i, pdNode* node_s, int norm, double* force, 
										 int& totalNumIntg, ofstream& fout)
{
	/*

	This function calculates the summation of bond force per unit volume for all family nodes of 
	colinear node "node_i" within the R+ space of the source node "node_s" in the direction of "norm".

	Input: 
	  1. node_s: pointer to the source node.
	  2. node_i: pointer to the colinear point in the given direction.
	  3. norm: given direction n.

	Output: 
	  1. force: array to store the three components of totoal bond force.
	  2. totalNumIntg: total number of integration points used in the calculation.

	*/

	double tol = _tolerance;
	double dx = _gridSpacing;
	// data for colinear node
	double x1i = node_i->GetX1();
	double x2i = node_i->GetX2();
	double x3i = node_i->GetX3();
	// loop over its family to find nodes which are within the R+ space of node_s
	double sum_forceij_x1 = 0.0; // sum of force
	double sum_forceij_x2 = 0.0;
	double sum_forceij_x3 = 0.0;
	totalNumIntg = 0;
	int numIntg;

	//int numfam = node_i->GetNumFamily();
	//for (int n=0;n!=numfam;++n)
	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{
		//pdNode* node_j = node_i->GetFamilyNode(n);
		pdNode* node_j = *itor;
		double x1j = node_j->GetX1();
		double x2j = node_j->GetX2();
		double x3j = node_j->GetX3();
		// check if this family node is in the R+ space of the source node node_s
		double j_plus;
		if (norm==1)
			j_plus = x1j - node_s->GetX1();
		else if (norm==2)
			j_plus = x2j - node_s->GetX2();
		else // norm==3
			j_plus = x3j - node_s->GetX3();
		//if (j_plus>=0.0)
		if (j_plus>=0.0 && node_j != node_i)
		{
			// calculate peridynamic force
			double pdForce[3];
			//pdBond* bond_ij = node_i->GetFamilyBond(n);			
			double dist_ij = sqrt(pow(x1j - x1i, 2.0) + pow(x2j - x2i, 2.0) + pow(x3j - x3i, 2.0)); 
			if (dist_ij > _horizon + 2.0 * _gridSpacing)
			{
				// a rough judge to skip calculation for two far away nodes
				continue;
			}
			else
			{
				// find the shortest and longest distance from node i to the cubic volume of node j
                double dist_min, dist_max;
                CalcDistanceToCube(node_i, node_j, dist_min, dist_max);
                if (dist_min < _horizon + _tolerance) // node j's cell has intersection volume with horizon sphere
                {
                    // only create needed bond to save memory					
                    pdBond* bond_ij = new pdBond(node_i, node_j); // here _numBonds serves as bond id
                    // save the min and max distance, set the initial bond breakage status to unbroken
                    bond_ij->SetMinDistance(dist_min);
                    bond_ij->SetMaxDistance(dist_max);
                    bond_ij->SetBreak(false);
					// add new bond to the container
					_bonds.push_back(bond_ij);
					if (_cciOn)
					{
						/*** CCI method (to compare with 1X1X1 Gaussian quadrature method) ***/
						CalcBondForceCCI(bond_ij, node_i, node_j, pdForce, numIntg);
						totalNumIntg += numIntg;
					}
					else
					{
						// Trapezoidal quadraqure or Gaussian quadraqure integration method
						CalcBondForceIntg(bond_ij, node_i, node_j, pdForce, numIntg);
						totalNumIntg += numIntg;
					}			
					//fout << node_i->GetID() << "\t" << node_j->GetID() 
					//	<< "\t" << forceij_x1 << "\t" << forceij_x2 << "\t" << forceij_x3 << endl;
					// sum
					sum_forceij_x1 += pdForce[0]; // force per unit volume
					sum_forceij_x2 += pdForce[1];
					sum_forceij_x3 += pdForce[2];
                }
			}			
		}
	}
	// output force
	force[0] = sum_forceij_x1;
	force[1] = sum_forceij_x2;
	force[2] = sum_forceij_x3;
}

void pdDataManager::CalcDistanceToCube(pdNode* node_i, pdNode* node_j, double& _min, double& _max)
{
	// this function finds the shortest distance from node i to the cubic volume posses by node j
	// to determine if node j is the family of node i
	// notice: node_i != node_j

	double delta = _horizon;
	double dx = _gridSpacing;
	double tol = _tolerance;
	// initial shortest and longest distance
	_min = MEGA;
	_max = TINY;
	// current position of two nodes. right now nodes don't have initial displacement (09/14/2009)
	double x1i = node_i->GetX1() + node_i->GetU1();
	double x2i = node_i->GetX2() + node_i->GetU2();
	double x3i = node_i->GetX3() + node_i->GetU3();
	double x1j = node_j->GetX1() + node_j->GetU1();
	double x2j = node_j->GetX2() + node_j->GetU2();
	double x3j = node_j->GetX3() + node_j->GetU3();
	// boundary limit of cubic volume of node j
	double x1lo = x1j - 0.5*dx;
	double x1hi = x1j + 0.5*dx;
	double x2lo = x2j - 0.5*dx;
	double x2hi = x2j + 0.5*dx;
	double x3lo = x3j - 0.5*dx;
	double x3hi = x3j + 0.5*dx;
	// calculate min and max distance to six faces
	double t_min, t_max; // temp store
	this->CalcDistanceToFace(1, x1lo, node_i, node_j, t_min, t_max);
	_min = min(_min, t_min);
	_max = max(_max, t_max);
	this->CalcDistanceToFace(1, x1hi, node_i, node_j, t_min, t_max);
	_min = min(_min, t_min);
	_max = max(_max, t_max);
	this->CalcDistanceToFace(2, x2lo, node_i, node_j, t_min, t_max);
	_min = min(_min, t_min);
	_max = max(_max, t_max);
	this->CalcDistanceToFace(2, x2hi, node_i, node_j, t_min, t_max);
	_min = min(_min, t_min);
	_max = max(_max, t_max);
	this->CalcDistanceToFace(3, x3lo, node_i, node_j, t_min, t_max);
	_min = min(_min, t_min);
	_max = max(_max, t_max);
	this->CalcDistanceToFace(3, x3hi, node_i, node_j, t_min, t_max);
	/**_min = min(*_min, t_min);
	*_max = max(*_max, t_max);*/
	_min = min(_min, t_min);
	_max = max(_max, t_max);
}

void pdDataManager::TestbdryEffFactor(pdNode* node_i, ofstream& fout)
{
	// This function tests force normalization factor for the source node node_i.
	
	double e11 = 1.0e-4;
	pdMaterial* mat_i = node_i->GetMaterial();
	double x1i = node_i->GetX1();
	double x2i = node_i->GetX2();
	double x3i = node_i->GetX3();
	double ecrit = mat_i->GetCriticalStretch();
	double spring = mat_i->GetSpringConstant();
	double elasticEnergyDensity = 0.0;
	// prescribed strain
	double eps11 = e11;
	double eps22 = e11;							
	double eps33 = e11;
	double eps12 = 0.0;
	double eps13 = 0.0;
	double eps21 = 0.0;
	double eps23 = 0.0;
	double eps31 = 0.0;
	double eps32 = 0.0;
	// create temp bc
	double val[9];
	val[0] = eps11;
	val[1] = eps12;
	val[2] = eps13;
	val[3] = eps21;
	val[4] = eps22;
	val[5] = eps23;
	val[6] = eps31;
	val[7] = eps32;
	val[8] = eps33;
	vector<double> value(val, val+sizeof(val)/sizeof(double));
	pdBdryCondition* bcTemp = new pdBdryConditionDisplGrad(-1, 1, 1, 1, value, 1.0, 2.0);
	// loop over its family nodes 
	int numfam = node_i->GetNumFamily();
	for (int n=0;n!=numfam;++n)
	{
		pdNode* node_j = node_i->GetFamilyNode(n);
		pdBond* bond = node_i->GetFamilyBond(n);
		double bdryEffij = 0.5*(node_i->GetBdryEffFactor() + node_j->GetBdryEffFactor());
		// copy the boundary condition of node j to a buffer
		pdBdryCondition* bcBuffer = node_j->GetBC();
		double elastic_j; // micropotential on node j
		if (!bond->IsBreak())
		{
			double x1j = node_j->GetX1();
			double x2j = node_j->GetX2();
			double x3j = node_j->GetX3();
			// xi
			double xi1 = x1j - x1i;
			double xi2 = x2j - x2i;
			double xi3 = x3j - x3i;						
			// calculate micropotential
			if (_cciOn)
			{
				// eta = strain*xi
				double eta1 = eps11*xi1 + eps12*xi2 + eps13*xi3;
				double eta2 = eps21*xi1 + eps22*xi2 + eps23*xi3;
				double eta3 = eps31*xi1 + eps32*xi2 + eps33*xi3;
				// r = abs(xi), p = abs(xi+eta)
				double r = sqrt(xi1*xi1 + xi2*xi2 + xi3*xi3);
				double p1 = xi1 + eta1;
				double p2 = xi2 + eta2;
				double p3 = xi3 + eta3;
				double p = sqrt(p1*p1 + p2*p2 + p3*p3);
				// calc micropotential based on PMB material
				double u = p - r;
				double stretch = u/r;
				double elastic;
				//if (r<=_horizon)
				//{
					if (u<=0 || (u>0 && stretch <= ecrit)) // bond in compression or not failed and in tension
						elastic = 0.5*(spring/r)*u*u;
					else // failed and in tension
						elastic = 0.0;
				//}
				//else
				//	elastic = 0.0;
				// calc volume reduction factor beta
				double beta = this->CalcBetaFactor(r, _horizon, _gridSpacing);
				double volj = beta * node_j->GetNodeVolume(); // the intersection volume between node i and j
				elastic_j = elastic * volj;
			}
			else
			{				
				node_j->SetBC(bcTemp);
				// use above given strain to calc displacement on the Gaussian points
				pdIntgManager* im = new pdIntgManager(_horizon, _gridSpacing, bond, node_i, node_j);
				im->SetIntgPoint(_intgPoint);
				im->SetTrapFlag(_aiOn);
				im->SetErrorControlFlag(_errorControlOn);
				im->SetErrorControlEPS(_eps);
				im->SetGaussianFlag(_fgiOn);
				im->SetGaussianAbscis(_abscis);
				im->SetGaussianWeight(_weights);
				im->SetMLSFlag(_mlsOn);
				
				double energy[1];
				if (_horizon <= bond->GetMinDistance()) // too far away (should never get here)
				{
					elastic_j = 0.0;
				}
				else if (_horizon >= bond->GetMaxDistance()) // cubic of node j is fully in the horizon of node i
				{
					im->CalcConfigurationType3("GetStrainEnergyDensity", energy);
					elastic_j = energy[0];
				}
				else
				{
					im->CalcConfigurationType3("GetStrainEnergyDensity", energy);
					elastic_j = energy[0];
				}
				delete im;
				node_j->SetBC(bcBuffer);
			}
			// sum up the elastic energy density 
			elasticEnergyDensity += 0.5 * elastic_j * bdryEffij;
		}
	}
	delete bcTemp;
	// print
	double p_grid = elasticEnergyDensity;
	double p_inf = 3.0*mat_i->GetYoungModulus()*e11*e11;
	fout << node_i->GetID() << "\t"
		<< node_i->GetX1() << "\t" << node_i->GetX2() << "\t" << node_i->GetX3() << "\t"
		<< node_i->GetNumFamily() << "\t" << node_i->GetBdryEffFactor() << "\t" << p_inf/p_grid << "\n";	
}

void pdDataManager::ReadExternalGridFile(ifstream &ext)
{
	// This function reads nodal position from an external grid file.

	ofstream outfile(_outFilePath.c_str(), ios::app);
	ext >> _numNodes;
	if (_myRank==0)
	{
		outfile << endl << "ReadExternalGridFile: total node number before deletion = " << _numNodes;
	}
	for (int i=0;i!=_numNodes;++i)
	{
		pdNode* node_i = new pdNode(i);
		double x1i, x2i, x3i;
		int mat_id;
		ext >> x1i >> x2i >> x3i >> mat_id;
		// check if this node is in a deletion space
		bool idel = false;
		for (vector<pdSpaceDeletion*>::const_iterator itor=_delSpaces.begin();itor!=_delSpaces.end();++itor)
		{
			// find space geometry type
			pdSpaceDeletion* dreg = *itor;
			int dreg_type = dreg->GetGeomType();
			// cuboid space, get boundary
			if (dreg_type==1)
			{
				double x1lo = dreg->GetX1Low();
				double x1hi = dreg->GetX1High();
				double x2lo = dreg->GetX2Low();
				double x2hi = dreg->GetX2High();
				double x3lo = dreg->GetX3Low();
				double x3hi = dreg->GetX3High();
				if (x1i>=x1lo && x1i<=x1hi && x2i>=x2lo && x2i<=x2hi
					&& x3i>=x3lo && x3i<=x3hi)
				{
					idel = true;
				}
			}
			// cylindrical space, get boundary
			else if (dreg_type==2)
			{
				double rad = dreg->GetRadius();
				double x1cen = dreg->GetX1Center();
				double x2cen = dreg->GetX2Center();
				double x3lo = dreg->GetX3Low();
				double x3hi = dreg->GetX3High();
				double dist = sqrt((x1i-x1cen)*(x1i-x1cen)+(x2i-x2cen)*(x2i-x2cen));
				if (dist<=rad && x3i>=x3lo && x3i<=x3hi)
				{
					idel = true;
				}
			}
		}
		if (!idel)
		{
			node_i->SetX1(x1i);
			node_i->SetX2(x2i);
			node_i->SetX3(x3i);
			node_i->SetMaterial(GetMaterial(mat_id));
		}
	}
	if (_myRank==0)
	{
		outfile << endl << "ReadExternalGridFile: total node number after deletion = " << _numNodes << endl;
	}
	outfile.close();
}

void pdDataManager::ReadEMUPlotFile(ifstream &finp, ofstream &fout) const
{
	// This function reads EMU plot file for nodal data.

	int node_id, nofail, nodbd, surfSegs, numNodes;
	double step, time, temp, x1tip, x2tip, x3tip, pe11, pe12, pe13, pe21, pe22, pe23;
	double x1, x2, x3, u1, u2, u3, v1, v2, v3, pe31, pe32, pe33;
	double x1clo, x1chi, x2clo, x2chi, x3clo, x3chi;
	double dmg, wt, yldfr, dil, anodty, dmgm, dmgi, ekn, vis, denst, volnod;
	double ecnode, stnode, bdryEff, scnode, edt, crit_exten, surf_dist;
	double crack_advance, crack_growth, dmg_fiber_tens, dmg_fiber_comp, dmg_matrix_tens;
	double dmg_matrix_comp, dmg_interlayer_tens, dmg_interlayer_comp, stiff, det_time, temp_node;
	string omitt, textline1, textline2, textline3;

	if (!finp)
	{
		if (_myRank==0)
		{
			cerr<<"ReadEMUPlotFile: EMU input file not exist!"<<endl;
		}
		exit(0);
	}

	// read the first line: total number of surface segment and node
	finp >> surfSegs >> numNodes;
	// read the run title and time
    getline(finp,textline1,'\n');
	getline(finp,textline2,'\n');
	getline(finp,textline3,'\n');
	// read penetrator data block
	finp >> step >> time;
	finp >> temp >> x1tip >> x2tip >> x3tip;
	finp >> pe11 >> pe12 >> pe13;
	finp >> pe21 >> pe22 >> pe23;
	finp >> pe31 >> pe32 >> pe33;
	finp >> x1clo >> x1chi >> x2clo >> x2chi;
	finp >> x3clo >> x3chi;
	fout << "\t" << step << "\t" << time << "\n";
	// read nodal data block
	for (int i=0;i!=numNodes;++i)
	{
		finp >> node_id >> temp >> temp >> x1 >> x2 >> x3 >> u1 >> u2 >> u3 >> v1 >> v2 >> v3 >> omitt;
		finp >> dmg >> wt >> yldfr >> dil >> anodty >> dmgm >> dmgi >> ekn >> nofail >> nodbd >> vis;
		finp >> ecnode >> stnode >> bdryEff >> scnode >> edt >> crit_exten >> surf_dist >> crack_advance;
		finp >> crack_growth >> denst >> volnod >> dmg_fiber_tens >> dmg_fiber_comp >> dmg_matrix_tens;
		finp >> dmg_matrix_comp >> dmg_interlayer_tens >> dmg_interlayer_comp >> stiff >> det_time >> temp_node;
		if (_myRank==0)
		{
			fout << node_id << "\t" << x1 << "\t" << x2 << "\t" << x3 
				<< "\t" << u1 << "\t" << u2 << "\t" << u3
				<< "\t" << v1 << "\t" << v2 << "\t" << v3 << "\n";
		}
	}
}

void pdDataManager::PostFormatData(ifstream& infile, ofstream& outfile) const
{
	/* This post-processing function reformats the output data.
	   Other datas are read by ReadInputFile() function.
	   All objects related to node need to be created here.
	 */

	int id, matID, bcID;
	double bdryEff, x1, x2, x3, u1, u2, u3, v1, v2, v3, dmg;
	string textline;
	// discard first four lines
	getline(infile,textline,'\n');
	getline(infile,textline,'\n');
	getline(infile,textline,'\n');
	getline(infile,textline,'\n');
	// read nodal data
	for (int i=0;i!=_numNodes;++i)
	{
		infile >> id >> x1 >> x2 >> x3 >> u1 >> u2 >> u3 >> v1 >> v2 >> v3
			>> matID >> bcID >> bdryEff >> dmg;
		pdNode* node_i = _nodes[id];
		// copy nodal data
		node_i->SetU1(u1);
		node_i->SetU2(u2);
		node_i->SetU3(u3);
		node_i->SetV1(v1);
		node_i->SetV2(v2);
		node_i->SetV3(v3);
		node_i->SetBdryEffFactor(bdryEff);
	}
	// write formatted data
	while(true)
	{
		cout << " Choose a format option: \n"
			<< " - Nodal data at every node on a cross-section plane (1) \n"			
			<< " - Nodal data at single node (2) \n"
			<< " - Nodal strain at every node on a cross-section plane (3) \n"
			<< " Your input: ";
		int opt;
		cin >> opt;
		if (opt==1)
		{
			int norm;
			double coord;
			cout << " Give the normal direction of the plane (1, 2, or 3): ";
			cin >> norm;
			cout << " Give the location of the plane: ";
			cin >> coord;
			outfile << " norm = " << norm << ", location = " << coord << endl;
			for (int i=0;i!=_numNodes;++i)
			{
				pdNode* node_i = _nodes[i];
				if (norm==1 && node_i->GetX1()==coord)
				{
					outfile << *(_nodes[i]);
					outfile << endl;
				}
				else if (norm==2 && node_i->GetX2()==coord)
				{
					outfile << *(_nodes[i]);
					outfile << endl;
				}
				else if (norm==3 && node_i->GetX3()==coord)
				{
					outfile << *(_nodes[i]);
					outfile << endl;
				}
			}
			outfile << endl;
		}
		else if (opt==2)
		{
			cout << " Give the node id: ";
			cin >> id;
			outfile << *(_nodes[id]);
			outfile << endl;
		}
		else if (opt==3)
		{
			int norm;
			double coord;
			cout << " Give the normal direction of the plane (1, 2, or 3): ";
			cin >> norm;
			cout << " Give the location of the plane: ";
			cin >> coord;
			outfile << " norm = " << norm << ", location = " << coord << endl;
			for (int i=0;i!=_numNodes;++i)
			{
				pdNode* node_i = _nodes[i];
				if (norm==1 && node_i->GetX1()==coord)
				{
					outfile << node_i->GetID()
						<< "\t" << node_i->GetX1() << "\t" << node_i->GetX2() << "\t" << node_i->GetX3();
					// calculate e11 strain
					int numFam = node_i->GetNumFamily();
					for (int j=0;j!=numFam;++j)
					{
						pdNode* node_j = node_i->GetFamilyNode(j);
						if (node_j->GetX1()==node_i->GetX1()+_gridSpacing
							&& node_j->GetX2()==node_i->GetX2()
							&& node_j->GetX3()==node_i->GetX3())
						{
							double e11 = (node_j->GetU1() - node_i->GetU1())/(node_j->GetX1() - node_i->GetX1());
							outfile << "\t" << e11;
							break;
						}
					}
					outfile << endl;
				}
				else if (norm==2 && node_i->GetX2()==coord)
				{
					outfile << node_i->GetID()
						<< "\t" << node_i->GetX1() << "\t" << node_i->GetX2() << "\t" << node_i->GetX3();
					// calculate e22 strain
					int numFam = node_i->GetNumFamily();
					for (int j=0;j!=numFam;++j)
					{
						pdNode* node_j = node_i->GetFamilyNode(j);
						if (node_j->GetX2()==node_i->GetX2()+_gridSpacing
							&& node_j->GetX1()==node_i->GetX1()
							&& node_j->GetX3()==node_i->GetX3())
						{
							double e22 = (node_j->GetU2() - node_i->GetU2())/(node_j->GetX2() - node_i->GetX2());
							outfile << "\t" << e22;
							break;
						}
					}
					outfile << endl;
				}
				else if (norm==3 && node_i->GetX3()==coord)
				{
					outfile << node_i->GetID()
						<< "\t" << node_i->GetX1() << "\t" << node_i->GetX2() << "\t" << node_i->GetX3();
					// calculate e33 strain
					int numFam = node_i->GetNumFamily();
					for (int j=0;j!=numFam;++j)
					{
						pdNode* node_j = node_i->GetFamilyNode(j);
						if (node_j->GetX3()==node_i->GetX3()+_gridSpacing
							&& node_j->GetX2()==node_i->GetX2()
							&& node_j->GetX1()==node_i->GetX1())
						{
							double e33 = (node_j->GetU3() - node_i->GetU3())/(node_j->GetX3() - node_i->GetX3());
							outfile << "\t" << e33;
							break;
						}
					}
					outfile << endl;
				}
			}
			outfile << endl;
		}
		// end of selection
		cout << " Another run? (y/n) ";
		string run;
		cin >> run;
		if (run=="n")
			break;
	}
}

void pdDataManager::PostAFD(ifstream& infile, ofstream& outfile)
{
	/* This post-processing function calculates areal force density at desired nodes.
	   Other datas are read by ReadInputFile() function.
	   All objects related to node need to be created here.
	 */

	int id, matID, bcID;
	double bdryEff, x1, x2, x3, u1, u2, u3, v1, v2, v3, dmg;
	string textline;
	// discard first four lines
	getline(infile,textline,'\n');
	getline(infile,textline,'\n');
	getline(infile,textline,'\n');
	getline(infile,textline,'\n');
	// read nodal data
	for (int i=0;i!=_numNodes;++i)
	{
		infile >> id >> x1 >> x2 >> x3 >> u1 >> u2 >> u3 >> v1 >> v2 >> v3
			>> matID >> bcID >> bdryEff >> dmg;
		pdNode* node_i = _nodes[id];
		// update nodal displacement and velocity
		node_i->SetU1(u1);
		node_i->SetU2(u2);
		node_i->SetU3(u3);
		node_i->SetV1(v1);
		node_i->SetV2(v2);
		node_i->SetV3(v3);
		node_i->SetBdryEffFactor(bdryEff);
	}
	// update the displacement of all the Gaussian points
	if (_mlsOn)
	{		
		for (vector<pdNode*>::iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
		{
			pdNode* node = *itor;
			CalcApproxDispl(node);
		}
	}
	else
	{
		if (_intgPoint==1) // one Gaussian point, so mls is turned off
		{
			// copy the nodal displacement to its Gaussian point since they are the same
			for (vector<pdNode*>::iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
			{
				pdNode* node = *itor;
				node->GetGaussianPoint(0)->SetU1(node->GetU1());
				node->GetGaussianPoint(0)->SetU2(node->GetU2());
				node->GetGaussianPoint(0)->SetU3(node->GetU3());
			}
		}
	}
	// calculate AFD
	while(true)
	{
		cout << " Choose an option to calculate the areal force density: \n"
			<< " - AFD at every node on a cross-section plane (1) \n"
			<< " - AFD at single node (2) \n"
			<< " Your input: ";
		int opt;
		cin >> opt;
		double sigma[9];
		if (opt==1)
		{
			int norm;
			double coord;
			cout << " Give the normal direction of the plane (1, 2, or 3): ";
			cin >> norm;
			cout << " Give the location of the plane: ";
			cin >> coord;
			outfile << " norm = " << norm << ", location = " << coord << endl;
			for (int i=0;i!=_numNodes;++i)
			{
				pdNode* node_i = _nodes[i];
				int numIntg;
				if (norm==1 && node_i->GetX1()==coord)
				{
					CalcArealForceDensity(node_i, sigma, numIntg, outfile);
					outfile << node_i->GetID() << "\t"
						<< node_i->GetX1() << "\t" << node_i->GetX2() << "\t" << node_i->GetX3() << "\t"
						<< sigma[0] << "\t" << sigma[1] << "\t" << sigma[2] << "\t"
						<< sigma[3] << "\t" << sigma[4] << "\t" << sigma[5] << "\t"
						<< sigma[6] << "\t" << sigma[7] << "\t" << sigma[8] << "\n";
				}
				else if (norm==2 && node_i->GetX2()==coord)
				{
					CalcArealForceDensity(node_i, sigma, numIntg, outfile);
					outfile << node_i->GetID() << "\t"
						<< node_i->GetX1() << "\t" << node_i->GetX2() << "\t" << node_i->GetX3() << "\t"
						<< sigma[0] << "\t" << sigma[1] << "\t" << sigma[2] << "\t"
						<< sigma[3] << "\t" << sigma[4] << "\t" << sigma[5] << "\t"
						<< sigma[6] << "\t" << sigma[7] << "\t" << sigma[8] << "\n";
				}
				else if (norm==3 && node_i->GetX3()==coord)
				{
					CalcArealForceDensity(node_i, sigma, numIntg, outfile);
					outfile << node_i->GetID() << "\t"
						<< node_i->GetX1() << "\t" << node_i->GetX2() << "\t" << node_i->GetX3() << "\t"
						<< sigma[0] << "\t" << sigma[1] << "\t" << sigma[2] << "\t"
						<< sigma[3] << "\t" << sigma[4] << "\t" << sigma[5] << "\t"
						<< sigma[6] << "\t" << sigma[7] << "\t" << sigma[8] << "\n";
				}
			}
		}
		else if (opt==2)
		{
			cout << " Give the node id: ";
			cin >> id;
			int numIntg;
			CalcArealForceDensity(_nodes[id], sigma, numIntg, outfile);
			outfile << " AFD at node " << _nodes[id]->GetX1() << "\t" << _nodes[id]->GetX2() << "\t" << _nodes[id]->GetX3() << endl;
			outfile << " sigma11 = " << sigma[0] << " sigma12 = " << sigma[1] << " sigma13 = " << sigma[2] << endl
				<< " sigma21 = " << sigma[3] << " sigma22 = " << sigma[4] << " sigma23 = " << sigma[5] << endl
				<< " sigma31 = " << sigma[6] << " sigma32 = " << sigma[7] << " sigma33 = " << sigma[8] << endl << endl;	
		}
		cout << " Another run? (y/n) ";
		string run;
		cin >> run;
		if (run=="n")
			break;
	}
}

string pdDataManager::GetProblemTitle() const
{
	return _title;
}

int pdDataManager::GetReadKeywordError() const
{
	return _rderr;
}

double pdDataManager::GetHorizon() const
{
	return _horizon;
}

int pdDataManager::GetNumNode() const
{
	return _numNodes;
}

int pdDataManager::GetNumBond() const
{
	return _numBonds;
}

int pdDataManager::GetNumMaterial() const
{
	return _numMats;
}

int pdDataManager::GetNumBC() const
{
	return _numBCs;
}

int pdDataManager::GetNumMatSpace() const
{
	return _numMatsSpaces;
}

int pdDataManager::GetNumDelSpace() const
{
	return _numDelSpaces;
}

int pdDataManager::GetNumBdSpace() const
{
	return _numBdSpaces;
}

double pdDataManager::GetGridCenterX1() const
{
	return _gridCenterX1;
}

double pdDataManager::GetGridCenterX2() const
{
	return _gridCenterX2;
}

double pdDataManager::GetGridCenterX3() const
{
	return _gridCenterX3;
}

int pdDataManager::GetGridDimensionX1() const
{
	return _gridDimensionX1;
}

int pdDataManager::GetGridDimensionX2() const
{
	return _gridDimensionX2;
}

int pdDataManager::GetGridDimensionX3() const
{
	return _gridDimensionX3;
}

double pdDataManager::GetGridSpacing() const
{
	return _gridSpacing;
}

double pdDataManager::GetBoundingX1Low() const
{
	return _gridBoundaryX1Low;
}

double pdDataManager::GetBoundingX1High() const
{
	return _gridBoundaryX1High;
}

double pdDataManager::GetBoundingX2Low() const
{
	return _gridBoundaryX2Low;
}

double pdDataManager::GetBoundingX2High() const
{
	return _gridBoundaryX2High;
}

double pdDataManager::GetBoundingX3Low() const
{
	return _gridBoundaryX3Low;
}

double pdDataManager::GetBoundingX3High() const
{
	return _gridBoundaryX3High;
}

int pdDataManager::GetTotalStep() const
{
	return _totalStep;
}

double pdDataManager::GetTotalTime() const
{
	return _totalTime;
}

int pdDataManager::GetStepDumpFrequency() const
{
	return _stepDumpFrequency;
}

bool pdDataManager::IsStateOn() const
{
	return _stateMatOn;
}

bool pdDataManager::IsExtGridFileOn() const
{
	return _extGridFileOn;
}

bool pdDataManager::IsFixedTimeStepOn() const
{
	return _fixedTimeStepOn;
}

pdNode* pdDataManager::GetNode(const int id) const
{
	return _nodes[id];
}

pdBond* pdDataManager::GetBond(const int id) const
{
	return _bonds[id];
}

pdMaterial* pdDataManager::GetMaterial(const int id) const
{
	return _mats[id].GetPtr();
}

pdBdryCondition* pdDataManager::GetBC(const int id) const
{
	return _BCs[id].GetPtr();
}

pdSpaceMaterial* pdDataManager::GetMaterialSpace(const int id) const
{
	return _matSpaces[id];
}

pdSpaceDeletion* pdDataManager::GetDelSpace(const int id) const
{
	return _delSpaces[id];
}

pdSpaceBdryCondition* pdDataManager::GetBdrySpace(const int id) const
{
	return _bdrySpaces[id];
}

pdSpacePlane* pdDataManager::GetPlane(const int id) const
{
	return _planes[id];
}

double pdDataManager::GetElastTotal() const
{
	return _elasticEnergyTotal;
}

double pdDataManager::GetKinetTotal() const
{
	return _kineticEnergyStep;
}

double pdDataManager::GetWorkBdryTotal() const
{
	return _externalWorkTotal;
}

double pdDataManager::GetElastStep() const
{
	return _elasticEnergyStep;
}

double pdDataManager::GetWorkBdryStep() const
{
	return _externalWorkStep;
}

double pdDataManager::GetX1force() const
{
	return _forcePlaneX1;
}

double pdDataManager::GetX2force() const
{
	return _forcePlaneX2;
}

double pdDataManager::GetX3force() const
{
	return _forcePlaneX3;
}

double pdDataManager::GetForce1() const
{
	return _forceX1;
}

double pdDataManager::GetForce2() const
{
	return _forceX2;
}

double pdDataManager::GetForce3() const
{
	return _forceX3;
}

bool pdDataManager::IsErrorControlOn() const
{
	return _errorControlOn;
}

bool pdDataManager::IsMLSOn() const
{
	return _mlsOn;
}

double pdDataManager::GetErroControlValue() const
{
	return _eps;
}

int pdDataManager::GetTrapIndex() const
{
	return _intgPoint;
}

void pdDataManager::PrintData(ofstream &out) const {}

void pdDataManager::PrintForceThroughPlane(ofstream &fout, int step) const
{
	for (vector<pdSpacePlane*>::const_iterator itor=_planes.begin();itor!=_planes.end();++itor)
	{
		pdSpacePlane* plane = *itor;
		if (_myRank==0)
		{
			fout << step << "\t"
				<< (plane->GetX1Low() + plane->GetX1High())/2 << "\t"
				<< (plane->GetX2Low() + plane->GetX2High())/2 << "\t"
				<< (plane->GetX3Low() + plane->GetX3High())/2 << "\t"
				<< plane->GetNorm() << "\t" <<plane->GetForce() << "\t" << plane->GetForce()/plane->GetArea() << endl;
		}
	}
}

void pdDataManager::SetWkDir()
{
	// this function is for user to define the working folder

	if (_myRank==0)
	{
		cout << " Please give the working folder: (no space between letters)" << endl;
	}
	cin >> _workDir;
	// erase the spaces before string head and any slash at the very end
	_workDir.erase(0, _workDir.find_first_not_of(' '));
	_workDir.erase(_workDir.find_last_not_of('\\')+1, string::npos);
	// start with the first slash
	string::iterator loc;
	loc = find(_workDir.begin(), _workDir.end(), '\\');
	while (loc!=_workDir.end())
	{
		_workDir.insert(loc, 1, '\\');
		loc = find(loc+2, _workDir.end(), '\\');
	}
	_workDir.insert(loc, 2, '\\');

	// set the default input and outprint file name
	_inFilePath = _workDir + "pds_input.txt";
	_outFilePath = _workDir + "pds_output.dat";
	// open the out data file
	ofstream outfile(_outFilePath.c_str(), ios::out);
	// check if the input file is in the working folder
	ifstream infile(_inFilePath.c_str(), ios::in);
	if(infile.fail())
	{
		if (_myRank==0)
		{
			cerr << "SetWkDir: input file not exits!" << endl;
		}
		exit(0);
	}
	infile.close();
}

string pdDataManager::GetWkDir() const
{
	return _workDir;
}

int pdDataManager::GetRank() const
{
	return _myRank;
}

string pdDataManager::KeywordAppend(const string keywd, const int i) const
{
	// This function appends string "i" to string "keywd" to create "keywdi"
	stringstream out;
	out << i;
	return (keywd + out.str());
}

int pdDataManager::CalcKroneckerDelta(const int i, const int j)
{
	// This function returns the Kronecker-delta value.
	return i==j;
}

double pdDataManager::CalcPermutationSymbol(const int i, const int j, const int k)
{
	// This function returns the permutation symbol value.
	return (0.5*(i-j)*(j-k)*(k-i));
}

void pdDataManager::CalcDistanceToFace(int norm, double x_face, pdNode* node_i, pdNode* node_j, double& dist_min, double& dist_max)
{
	// this function serves a sub-function for pdDataManager::CalcDistanceToCube.
	// It calculates the shortest and longest distance from node i to one face of the cubic of node j
	// input: norm is the direction of the norm for that face (either 1, 2 or 3)
	//           x_face is the RHS value of the surface function for that face (i.e. if norm=1, x_face is either x1j+0.5*dx or x1j-0.5*dx)
	// output: minimum and maximum distance from node i to that face

	double tol = 1e-4*_gridSpacing;
	int n1, n2, n3;
	if (norm==1)
	{
		n1 = 1;
		n2 = n3 = 0;
	}
	else if (norm==2)
	{
		n1 = n3 = 0;
		n2 = 1;
	}
	else
	{
		n1 = n2 = 0;
		n3 = 1;
	}
	double _x1i = node_i->GetX1() + node_i->GetU1();
    double _x2i = node_i->GetX2() + node_i->GetU2();
    double _x3i = node_i->GetX3() + node_i->GetU3();
    double _x1j = node_j->GetX1() + node_j->GetU1();
    double _x2j = node_j->GetX2() + node_j->GetU2();
    double _x3j = node_j->GetX3() + node_j->GetU3();
	// draw a line from node i to node j, check if it is intersecting with any of the six faces of the cubic
	if (abs((_x1i-_x1j)*n1+(_x2i-_x2j)*n2+(_x3i-_x3j)*n3)<tol)
	{
		// the line from node i to j is parallel to this plane, no intersecting point
		dist_min = MEGA;
		dist_max = TINY;
		return;
	}
	else
	{
		// boundary limit of cubic volume of node j
		double x1lo = _x1j - 0.5*_gridSpacing;
		double x1hi = _x1j + 0.5*_gridSpacing;
		double x2lo = _x2j - 0.5*_gridSpacing;
		double x2hi = _x2j + 0.5*_gridSpacing;
		double x3lo = _x3j - 0.5*_gridSpacing;
		double x3hi = _x3j + 0.5*_gridSpacing;
		// find the intersecting point and check if it is within the face limits
		double x_i = _x1i*n1 + _x2i*n2 + _x3i*n3;
		double x_j = _x1j*n1 + _x2j*n2 + _x3j*n3;
		double fac = (x_face - x_i)/(x_j - x_i);
		// coord for the intersecting point
		double x1int = _x1i + fac*(_x1j - _x1i);
		double x2int = _x2i + fac*(_x2j - _x2i);
		double x3int = _x3i + fac*(_x3j - _x3i);
		if (x1int>x1lo-tol && x1int<x1hi+tol && x2int>x2lo-tol && x2int<x2hi+tol && x3int>x3lo-tol && x3int<x3hi+tol)
		{
			// intersecting point is within the face limits, so the shortest path is from node i to one point on this face
			// draw a perpendicular line from node i to this face, find the foot
			double x1foot = (1-n1)*_x1i + n1*x_face;
			double x2foot = (1-n2)*_x2i + n2*x_face;
			double x3foot = (1-n3)*_x3i + n3*x_face;
			// find the relative location of the foot point to this face and calculate the relative distance
			double x1s = (1-n1)*GetShortestBetween(x1foot, x1lo, x1hi);
			double x2s = (1-n2)*GetShortestBetween(x2foot, x2lo, x2hi);
			double x3s = (1-n3)*GetShortestBetween(x3foot, x3lo, x3hi);
			double x1l = (1-n1)*GetLongestBetween(x1foot, x1lo, x1hi);
			double x2l = (1-n2)*GetLongestBetween(x2foot, x2lo, x2hi);
			double x3l = (1-n3)*GetLongestBetween(x3foot, x3lo, x3hi);
			// shortest and longest distance from foot to this face
			double ff_s = sqrt(x1s*x1s + x2s*x2s + x3s*x3s);
			double ff_l = sqrt(x1l*x1l + x2l*x2l + x3l*x3l);
			// distance from node i to foot
			double i_to_foot = sqrt((x1foot-_x1i)*(x1foot-_x1i) + (x2foot-_x2i)*(x2foot-_x2i) + (x3foot-_x3i)*(x3foot-_x3i));
			// shortest and longest distance from node i to this plane
			dist_min = sqrt(ff_s*ff_s + i_to_foot*i_to_foot);
			dist_max = sqrt(ff_l*ff_l + i_to_foot*i_to_foot);
			return;
		}
		else
		{
			// intersecting point is not within the plane limits
			dist_min = MEGA;
			dist_max = TINY;
			return;
		}
	}
}

double pdDataManager::GetShortestBetween(double x_in, double x_lo, double x_hi) const
{
	// Given two boundary limits x_lo and x_hi on any axis
	// This function returns the shortest distance from a point x_in to the boundary

	double tol = 1e-4*_gridSpacing;
	if (x_in<x_lo-tol)
	{
		return abs(x_in - x_lo);
	}
	else if (x_in>x_hi+tol)
	{
		return abs(x_in - x_hi);
	}
	else
	{
		return 0.0;
	}
}

double pdDataManager::GetLongestBetween(double x_in, double x_lo, double x_hi) const
{
	// Given two boundary limits x_lo and x_hi on any axis
	// This function returns the longest distance from a point x_in to the boundary

	double tol = 1e-4*_gridSpacing;
	if (x_in<x_lo-tol)
	{
		return abs(x_in - x_hi);
	}
	else if (x_in>x_hi+tol)
	{
		return abs(x_in - x_lo);
	}
	else // in the middle
	{
		if (x_in > 0.5*(x_hi + x_lo))
			return abs(x_in - x_lo);
		else
			return abs(x_in - x_hi);
	}
}

void pdDataManager::SetGaussianPointsValue()
{
	_abscis.clear(); // in case this function is called recursively
	_weights.clear();
	switch (_intgPoint)
	{
	case 2:
		{
			_abscis.push_back(1.0/sqrt(3.0));
			_abscis.push_back(-1.0/sqrt(3.0));
			_weights.push_back(1.0);
			_weights.push_back(1.0);
			break;
		}
	case 3:
		{
			_abscis.push_back(0.0);
			_abscis.push_back(sqrt(3.0/5.0));
			_abscis.push_back(-sqrt(3.0/5.0));
			_weights.push_back(8.0/9.0);
			_weights.push_back(5.0/9.0);
			_weights.push_back(5.0/9.0);
			break;
		}
	case 4:
		{
			_abscis.push_back(sqrt((3.0 - 2.0*(sqrt(6.0/5.0)))/7.0));
			_abscis.push_back(-sqrt((3.0 - 2.0*(sqrt(6.0/5.0)))/7.0));
			_abscis.push_back(sqrt((3.0 + 2.0*(sqrt(6.0/5.0)))/7.0));
			_abscis.push_back(-sqrt((3.0 + 2.0*(sqrt(6.0/5.0)))/7.0));
			_weights.push_back((18.0 + sqrt(30.0))/36.0);
			_weights.push_back((18.0 + sqrt(30.0))/36.0);
			_weights.push_back((18.0 - sqrt(30.0))/36.0);
			_weights.push_back((18.0 - sqrt(30.0))/36.0);
			break;
		}
	case 5:
		{
			_abscis.push_back(0.0);
			_abscis.push_back(sqrt(5.0 - 2.0*sqrt(10.0/7.0))/3.0);
			_abscis.push_back(-sqrt(5.0 - 2.0*sqrt(10.0/7.0))/3.0);
			_abscis.push_back(sqrt(5.0 + 2.0*sqrt(10.0/7.0))/3.0);
			_abscis.push_back(-sqrt(5.0 + 2.0*sqrt(10.0/7.0))/3.0);
			_weights.push_back(128.0/225.0);
			_weights.push_back((322.0 + 13.0*sqrt(70.0))/900.0);
			_weights.push_back((322.0 + 13.0*sqrt(70.0))/900.0);
			_weights.push_back((322.0 - 13.0*sqrt(70.0))/900.0);
			_weights.push_back((322.0 - 13.0*sqrt(70.0))/900.0);
			break;
		}
	case 6:
		{
			_abscis.push_back(0.2386191860831969086305017);
			_abscis.push_back(-0.2386191860831969086305017);
			_abscis.push_back(0.6612093864662645136613996);
			_abscis.push_back(-0.6612093864662645136613996);
			_abscis.push_back(0.9324695142031520278123016);
			_abscis.push_back(-0.9324695142031520278123016);
			_weights.push_back(0.4679139345726910473898703);
			_weights.push_back(0.4679139345726910473898703);
			_weights.push_back(0.3607615730481386075698335);
			_weights.push_back(0.3607615730481386075698335);
			_weights.push_back(0.1713244923791703450402961);
			_weights.push_back(0.1713244923791703450402961);
			break;
		}
	case 8:
		{
			_abscis.push_back(0.1834346424956498049394761);
			_abscis.push_back(0.5255324099163289858177390);
			_abscis.push_back(0.7966664774136267395915539);
			_abscis.push_back(0.9602898564975362316835609);
			_abscis.push_back(-0.1834346424956498049394761);
			_abscis.push_back(-0.5255324099163289858177390);
			_abscis.push_back(-0.7966664774136267395915539);
			_abscis.push_back(-0.9602898564975362316835609);
			_weights.push_back(0.3626837833783619829651504);
			_weights.push_back(0.3137066458778872873379622);
			_weights.push_back(0.2223810344533744705443560);
			_weights.push_back(0.1012285362903762591525314);
			_weights.push_back(0.3626837833783619829651504);
			_weights.push_back(0.3137066458778872873379622);
			_weights.push_back(0.2223810344533744705443560);
			_weights.push_back(0.1012285362903762591525314);
			break;
		}
	default: // n=1
		{
			_abscis.push_back(0.0);
			_weights.push_back(2.0);
			break;
		}
	}
}

double pdDataManager::CalcBetaFactor(double dist, double delta, double dx)
{
	/* This function calculated the boundary correct factor beta based on the inputs. 
	It is a one-point integration method originated from the EMU code.
	   dist: the distance between the source node and the family node 
	   delta: the horizon size
	   dx: the grid size (or the grid spacing)
	*/

	double innerBdry = delta - 0.5*dx;
	double outerBdry = delta + 0.5*dx; // modified CCI beta to inlucde all intersecting nodes as family
	if (_cciOn && !_cciMdf)
	{
		outerBdry = delta;  // original CCI beta
	}
	
	// beta calculation
	if (dist <= innerBdry)
	{
		return 1.0;
	}
	else if (dist > innerBdry && dist <= outerBdry)
	{
		return (delta + 0.5*dx - dist) / dx;
	}
	else
	{
		return 0.0;	
	}

}