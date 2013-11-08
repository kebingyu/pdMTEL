#include "pdDataManager.h"
//#include <stdlib.h>
//#include "random.h"
using namespace std;

double pdDataManager::CalcStableTimeStep()
{
	// This function calculates stable time step. If the "fixe_CalcStableTimeStep" keyword is shown in the input file, use it instead.
	return _fixedTimeStepOn ? _fixedTimeStep : _safeFactor*_minBondLength/_maxBulkSoundSpeed;
}

void pdDataManager::InitBond()
{
	// This function finds family nodes and set up all bonds.
	
	double tol = _tolerance;
	double dx = _gridSpacing;
	double delta = _horizon;
	// loop over all real nodes to find and initial all real bonds
	_minBondLength = _horizon;
	for (vector<pdNode*>::const_iterator itor_i=_nodes.begin();itor_i!=_nodes.end();++itor_i)
	{
		pdNode* node_i = *itor_i;
		double x1i = node_i->GetX1();
		double x2i = node_i->GetX2();
		double x3i = node_i->GetX3();
		// loop over nodes after node i (prevent duplicate searching)
		for (vector<pdNode*>::const_iterator itor_j=itor_i+1;itor_j!=_nodes.end();++itor_j)
		{
			pdNode* node_j = *itor_j;
			double x1j = node_j->GetX1();
			double x2j = node_j->GetX2();
			double x3j = node_j->GetX3();
			double dist_ij = sqrt(pow(x1j - x1i, 2.0) + pow(x2j - x2i, 2.0) + pow(x3j - x3i, 2.0)); 
			if (_horizon < dist_ij - 1.5*_gridSpacing)
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
                    // add new bond to the container
                    pdBond* bond_ij = new pdBond(node_i, node_j);
                    // save the min and max distance, set the initial bond breakage status to unbroken
                    bond_ij->SetMinDistance(dist_min);
                    bond_ij->SetMaxDistance(dist_max);
                    bond_ij->SetBreak(false);
                    _bonds.push_back(bond_ij);
                    // find the shortest bond length
                    _minBondLength = min(_minBondLength, dist_ij);
                    // store ptr to this bond to node i and j family container
                    node_i->AddFamilyBond(bond_ij);
                    node_j->AddFamilyBond(bond_ij);
                }
			}
		}
	}
	_numBonds = int(_bonds.size());
	ofstream outfile(_outFilePath.c_str(), ios::app);
	if (_myRank==0)
	{
		outfile << endl << "InitBond: total bonds number = " << _numBonds << endl;
		outfile << endl << "InitBond: bonds initialization done " << endl;
	}
}

void pdDataManager::SetGaussianPoint()
{
	// This function initials Gaussian points list for each node
	// (1) select Gaussian point abscissas and weights
	SetGaussianPointsValue();
	// (2) loop over node list
	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{
		pdNode* node_i = *itor;
		node_i->ClearGaussianPointList(); // in case this function is called recursively
		for (int i=0; i!=_intgPoint; ++i) // X1 direction
		{
			for (int j=0; j!=_intgPoint; ++j) // X2 direction
			{
				for (int k=0; k!=_intgPoint; ++k) // X3 direction
				{
					pdGauss* gPoint = new pdGauss();
					// Set global coordinates
					gPoint->SetX1(0.5*_gridSpacing*_abscis[i] + node_i->GetX1());
					gPoint->SetX2(0.5*_gridSpacing*_abscis[j] + node_i->GetX2());
					gPoint->SetX3(0.5*_gridSpacing*_abscis[k] + node_i->GetX3());
					gPoint->SetAbscissa(_abscis[i], _abscis[j], _abscis[k]);
					gPoint->SetWeight(_weights[i], _weights[j], _weights[k]);
					// currently set all Gaussian points to have zero initial displacement
					gPoint->SetU1(0.0);
					gPoint->SetU2(0.0);
					gPoint->SetU3(0.0);
					node_i->AddGaussianPoint(gPoint);
				}
			}
		}
	}
}

void pdDataManager::InitIntegrationPoint()
{
	/* 
	
	This function initials the integration points.
      1. If _mlsOn = true, set up the domain of influence for each node.
	  2. For fixed Gaussian integration, create Gaussian points and calculate the MLS coefficients for each Gaussian point if _mlsOn=ture.
	  3. For adaptive integration, calculate the MLS coefficients for each Gaussian point if _mlsOn=ture. 
	     The trapezoidal points are not actually created. Only the MLS coefficients for each set of trapezoidal point are stored.

	*/

	if (_mlsOn)
	{
		// Set up domain of influence for each node if MLS approximation is used.
		for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
		{
			pdNode* node_i = *itor;
			node_i->AddContribNode(node_i);
			double x1i = node_i->GetX1();
			double x2i = node_i->GetX2();
			double x3i = node_i->GetX3();
			/* 
			user need to set the horizon big enough so that it contains enough family nodes
			for a default horizon (equals 3 grid spacing), it is big enough
			nodes near the center have 250 family nodes, while nodes at the 3-edge corner have 50
			 */
			int numFam = node_i->GetNumFamily(); 
			for (int j = 0; j != numFam; ++j)
			{
				pdNode* node_j = node_i->GetFamilyNode(j);
				double x1j = node_j->GetX1();
				double x2j = node_j->GetX2();
				double x3j = node_j->GetX3();
				// get relative distance, right now no initial displacement boundary condition (07/29/2010)
				double dist_ij = sqrt(pow(x1j - x1i, 2.0) + pow(x2j - x2i, 2.0) + pow(x3j - x3i, 2.0)); 
				if (dist_ij <= _rw)
				{
					node_i->AddContribNode(node_j);
				}
			}
			// increase support domain size if not enough contributing nodes
			if (node_i->GetContribNodeListSize() < BASISDIM)
			{
				double _rwIncreased = _rw * 1.5; 
				node_i->ClearContribNodeList();
				node_i->AddContribNode(node_i);
				for (int j = 0; j != numFam; ++j)
				{
					pdNode* node_j = node_i->GetFamilyNode(j);
					double x1j = node_j->GetX1();
					double x2j = node_j->GetX2();
					double x3j = node_j->GetX3();
					double dist_ij = sqrt(pow(x1j - x1i, 2.0) + pow(x2j - x2i, 2.0) + pow(x3j - x3i, 2.0));
					if (dist_ij <= _rwIncreased)
					{
						node_i->AddContribNode(node_j);
					}
				}
			}
		}
	}

	if (_fgiOn)
	{
		// 1. Set up Gaussian points.
		SetGaussianPoint();
		if (_mlsOn)
		{
			// 2. calculate the coefficients of contributing nodes for each node using MLS approximation.
			// This is done only once before the dynamic solution starts.
			for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
			{
				pdNode* node = *itor;
				// get the contributing node list of this node 
				vector<pdNode*> cbList;
				node->GetContribNodeList(cbList);
				int gPointNum = node->GetGaussianPointListSize();
				for (int g=0;g!=gPointNum;++g)
				{
					pdGauss* gPoint = node->GetGaussianPoint(g);
					// get the global coordinate of Gaussian point
					double x = gPoint->GetX1();
					double y = gPoint->GetX2();
					double z = gPoint->GetX3();
					int dim = BASISDIM; // basis vector dimension
					int ndNumber = cbList.size(); // total number of contributing nodes
					// set up moment matrix (lhs)
					pdMatrix* lhs = new pdMatrix(dim, dim);
					lhs->Zero();
					for (int k = 0; k != ndNumber; ++k) // loop over contributing nodes list
					{
						pdVector* basisVec = this->GetBasis(cbList[k]->GetX1(), cbList[k]->GetX2(), cbList[k]->GetX3(), dim);		
						double weight = this->GetWeight(cbList[k], x, y, z, cbList[k]->GetSupportDomainSize());
						//cout << cbList[k]->GetID() << "\t" << weight << endl;
						for (int i = 0; i != dim; ++i)
						{
							for (int j = 0; j != dim; ++j)
							{				
								lhs->AddCoeff(i, j, weight * basisVec->GetCoeff(i) * basisVec->GetCoeff(j));
							}
						}
						delete basisVec;
					}

					// set up ds (rhs)
					pdVector* rhs = this->GetBasis(x, y, z, dim);
					pdVector* sol = new pdVector(dim);
					// solve system of equations, store the results in sol
					lhs->Gauss(rhs, sol);
					// check the condition numbers
					//double normMomentMatrix = lhs.RowSumNorm(); // infinite norm of the moment matrix
					//double normMomentMatrix = lhs->FBnorm(); // Frobenius norm
					//double normMomentMatrixInv = lhs.Norm(); // row-sum norm
					//double normMomentMatrixInv = lhs->FBnorm();
					//if (normMomentMatrix * normMomentMatrixInv > ILLCONDITION)
					//{
					//	//// for debug
					//	//fout << node->GetID() 
					//	//	<< "\t" << node->GetX1() << "\t" << node->GetX2() << "\t" << node->GetX3()
					//	//	<< "\t" << normMomentMatrix * normMomentMatrixInv << "\t" << 0 << endl;

					//	cout << "Moment matrix is ill-conditioned at node " << node->GetID() << endl;
					//	exit(0);
					//}
					//else
					//{
						// set up B matrix
						pdMatrix* b = new pdMatrix(0, 0);
						for (int i = 0; i != ndNumber; ++i)
						{
							pdVector* basisVec = this->GetBasis(cbList[i]->GetX1(), cbList[i]->GetX2(), cbList[i]->GetX3(), dim);
							basisVec->Mult(this->GetWeight(cbList[i], x, y, z, cbList[i]->GetSupportDomainSize()));
							b->AddColumn(basisVec);
						}
						pdMatrix* b2 = new pdMatrix(b->GetNumCols(), b->GetNumRows());
						b->Trans(b2);
						// store solved coefficient in Vector coeff
						pdVector* coeff = b2->Mult(sol);
						gPoint->SetCoeff(coeff);
						// delete temporary objects	
						delete lhs;
						delete rhs;
						delete sol;
						delete b;
						delete b2;
					//}
				}
			}
		}
	}
	else if (_aiOn)
	{
		if (_mlsOn)
		{
			// run the trapezoidal integration scheme once to calculate the MLS coefficient for each set of trapezoidal points
			double pdForce[3];
			for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
			{
				pdNode* node_i = *itor;
				int numFam = node_i->GetNumFamily();
				for (int j=0;j!=numFam;++j)
				{
					pdNode* node_j = node_i->GetFamilyNode(j);
					pdBond* bond_ij = node_i->GetFamilyBond(j);
					int numIntg = 0; // total number of integration points
					pdIntgManager* im = new pdIntgManager(_horizon, _gridSpacing, bond_ij, node_i, node_j);
					im->SetIntgPoint(_intgPoint);
					im->SetTrapFlag(true);
					im->SetErrorControlFlag(false);
					im->SetErrorControlEPS(0.0);
					im->SetGaussianFlag(false);
					im->SetMLSFlag(_mlsOn);
					im->SetInitMLSCoeffFlag(true);
					// use the unit integrand ("UnitValue") to calculate the MLS coefficient
					if (_horizon <= bond_ij->GetMinDistance()) // too far away (should never get here since we are doing calculation in node i's family)
					{
					}
					else if (_horizon >= bond_ij->GetMaxDistance()) // cell of node j is fully in the horizon sphere
					{			
						im->CalcConfigurationType3("UnitValue", pdForce);
						numIntg += im->GetEndIntgPoint(1)*im->GetEndIntgPoint(2)*im->GetEndIntgPoint(3);
					}
					else // node j's cell has intersection with the horizon sphere
					{
						im->AdaptiveIntegration("UnitValue", pdForce);
						numIntg += im->GetEndIntgPoint(1)*im->GetEndIntgPoint(2)*im->GetEndIntgPoint(3);
		
					}
					delete im;
				}
			}
		}
	}
	else
	{
		// implement other integration methods
	}
	ofstream outfile(_outFilePath.c_str(), ios::app);
	if (_myRank==0)
	{
		outfile << endl << "InitIntegrationPoint: integration points initialization done " << endl;
	}
}

void pdDataManager::InitProcBdry()
{
	// This function splits the whole space (including boundary margin) into cells owned by each processor (kebing 07/25/2009)
		
	// increments
	double dx1 = (_gridBoundaryX1High - _gridBoundaryX1Low)/_proc1;
	double dx2 = (_gridBoundaryX2High - _gridBoundaryX2Low)/_proc2;
	double dx3 = (_gridBoundaryX3High - _gridBoundaryX3Low)/_proc3;
	double high_limit, low_limit;

	// find processor boundaries in the x1-direction
    vector<double> bdry_x1_low, bdry_x1_high;	// containers for lower and upper limit 
	for (int ip=0;ip!=_proc1;++ip)	 // processor number starts from 0
	{
      if (ip==0)
	  {
		  low_limit = _gridBoundaryX1Low;
	  }
	  else
	  {
		  low_limit = bdry_x1_high[ip-1];
	  }
	  bdry_x1_low.push_back(low_limit);
	  high_limit = bdry_x1_low[ip] + dx1;
	  if (ip==_proc1-1 && high_limit<_gridBoundaryX1High)
	  {
		  high_limit = _gridBoundaryX1High;
	  }
	  bdry_x1_high.push_back(high_limit);
	}
		  
	// find processor boundaries in the x2-direction
    vector<double> bdry_x2_low, bdry_x2_high;	// containers for lower and upper limit 
	for (int ip=0;ip!=_proc2;++ip)	 // processor number starts from 0
	{
      if (ip==0)
	  {
		  low_limit = _gridBoundaryX2Low;
	  }
	  else
	  {
		  low_limit = bdry_x2_high[ip-1];
	  }
	  bdry_x2_low.push_back(low_limit);
	  high_limit = bdry_x2_low[ip] + dx2;
	  if (ip==_proc2-1 && high_limit<_gridBoundaryX2High)
	  {
		  high_limit = _gridBoundaryX2High;
	  }
	  bdry_x2_high.push_back(high_limit);
	}

	// find processor boundaries in the x3-direction
    vector<double> bdry_x3_low, bdry_x3_high;	// containers for lower and upper limit 
	for (int ip=0;ip!=_proc3;++ip)	 // processor number starts from 0
	{
      if (ip==0)
	  {
		  low_limit = _gridBoundaryX3Low;
	  }
	  else
	  {
		  low_limit = bdry_x3_high[ip-1];
	  }
	  bdry_x3_low.push_back(low_limit);
	  high_limit = bdry_x3_low[ip] + dx3;
	  if (ip==_proc3-1 && high_limit<_gridBoundaryX3High)
	  {
		  high_limit = _gridBoundaryX3High;
	  }
	  bdry_x3_high.push_back(high_limit);
	}
	
	// save the boundary limits for each processor
	_myBdryX1Low = bdry_x1_low[_myRank];
	_myBdryX1High = bdry_x1_high[_myRank];
	_myBdryX2Low = bdry_x2_low[_myRank];
	_myBdryX2High = bdry_x2_high[_myRank];
	_myBdryX3Low = bdry_x3_low[_myRank];
	_myBdryX3High = bdry_x3_high[_myRank];	

	// set needed boundary limits for each processor
	// (discuss:1. boundary limits at the edge may be same as need bdry 2. needbdry may exceed exceed other bdry limit)
	_myNeedBdryX1Low = _myBdryX1Low - _horizon;
	_myNeedBdryX1High = _myBdryX1High + _horizon;
	_myNeedBdryX2Low = _myBdryX2Low - _horizon;
	_myNeedBdryX2High = _myBdryX2High + _horizon;
	_myNeedBdryX3Low = _myBdryX3Low - _horizon;
	_myNeedBdryX3High = _myBdryX3High + _horizon;

	// set up table of what processors need to communicate
	bdry_x1_low.clear();
	bdry_x1_high.clear();
	bdry_x2_low.clear();
	bdry_x2_high.clear();
	bdry_x3_low.clear();
	bdry_x3_high.clear();
}

void pdDataManager::InitDynamicSolve()
{
	/* This function sets all energy and bond force storage to zero.
	   Also put the velocity of previous step to the velocity buff container.
	   This function is called at the beginning of each time step.
	*/

	_elasticEnergyStep = 0.0;
	_kineticEnergyStep = 0.0;
	_forceX1 = 0.0;
	_forceX2 = 0.0;
	_forceX3 = 0.0;
	// loop over all nodes
	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{
		pdNode* node_i = *itor;
		node_i->SetNumDmgBond(0); // recount damaged bond at each time step (05/30/09)
		node_i->SetNumYldBond(0);
		node_i->SetElastEnergyDensity(0.0);
		node_i->SetKinetEnergy(0.0);
		// store the velocity of the previous half time step into buffer
		node_i->SetV1Buffer(node_i->GetV1());
		node_i->SetV2Buffer(node_i->GetV2());
		node_i->SetV3Buffer(node_i->GetV3());
	}
	//for (vector<pdSpacePlane*>::const_iterator itor=_planes.begin();itor!=_planes.end();++itor)
	//{
	//	(*itor)->SetForce(0.0);
	//}
}

//void pdDataManager::CalcStateBasedForce(pdNode* node_i, const double dt)
//{
//	// This function calculates total force density per unit volume square due to peridynamic
//	// state interaction between node mi and its family node
//
////	int i, j, k, mj, n, numfam;
//	double xi[3], xj[3], ui[3], uj[3], vi[3], vj[3], vel_gra_eul[3][3], rate_of_def[3][3];
//	double dv[3][3], spin[3][3], z[3], w[3], temp[3][3], temp_inv[3][3], omega_vec[3];
//	double omega_tens[3][3], tempa[3][3], tempb[3][3], tempb_inv[3][3], rot_tens_old[3][3];
//	double scratch_matrix[3][3], rate_of_def_unrot[3][3], strain_inc[3][3], strain_inc_dev[3][3];
//	double stress_unrot_old[3][3], stress_unrot_new[3][3], strain_unrot[3][3], stress_piola[3][3], strain_cauchy[3][3];
////	double r, p, ecrit, det_shape, det_def_grad, trace_of_v, det_temp, det_tempb, detf;
////	double elastic, young_mod, poisson, bulk_mod, shear_mod, dilat_inc, dmg_zero_stress, fac;
////	double force_state_1, force_state_2, force_state_3, influence_state;
//	ofstream outfile(_outFilePath.c_str(), ios::app);
//
//	// find the current status of node mi
//	xi[0] = node_i->GetX1();
//	xi[1] = node_i->GetX2();
//	xi[2] = node_i->GetX3();
//	ui[0] = node_i->GetU1();
//	ui[1] = node_i->GetU2();
//	ui[2] = node_i->GetU3();
//	vi[0] = node_i->GetV1();
//	vi[1] = node_i->GetV2();
//	vi[2] = node_i->GetV3();
//	int numfam = node_i->GetNumFamily();
//	pdMaterial* mat_i = node_i->GetMaterial();
//
//	// ****** Step I ******
//	// Loop over family nodes to find:
//	// 1) deformation gradient tensor F 
//	// 2) displacement gradient tensor 
//	// 3) shape tensor inverse inv(K) 
//	// 4) velocity gradient tensor Fdot
//	// This step is same as "def_gradient" subroutine of state-based CCI.
//	memset(shape, 0, sizeof(shape));
//	memset(disp_grad, 0, sizeof(disp_grad));
//	memset(vel_grad, 0, sizeof(vel_grad));
//	memset(sum_disp_grad, 0, sizeof(sum_disp_grad));
//	memset(sum_vel_grad, 0, sizeof(sum_vel_grad));
//	memset(def_grad, 0, sizeof(def_grad));
//	memset(def_gra_inv, 0, sizeof(def_gra_inv));
//	double det_shape = 0.;
//	double det_def_grad = 0.;
//	double detf = 1.;
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			ident[i][j] = pdIntgManager::CalcKroneckerDelta(i, j);
//			def_grad[i][j] = ident[i][j];
//			def_gra_inv[i][j] = ident[i][j];
//			shape_inv[i][j] = ident[i][j];
//		}
//	}
//
//	// loop over family nodes
//	for (int n=0;n!=numfam;++n)
//	{
//		pdNode* node_j = node_i->GetFamilyNode(n);
//		pdBond* bond_ij = node_i->GetFamilyBond(n);
//		{
//			// find the current status of node mj
//			xj[0] = node_j->GetX1();
//			xj[1] = node_j->GetX2();
//			xj[2] = node_j->GetX3();
//			uj[0] = node_j->GetU1();
//			uj[1] = node_j->GetU2();
//			uj[2] = node_j->GetU3();
//			vj[0] = node_j->GetV1();
//			vj[1] = node_j->GetV2();
//			vj[2] = node_j->GetV3();
//			// find the state parameters between node mi and mj
//			for (int k=0;k!=3;++k)
//			{
//				refpos_state[k] = xj[k] - xi[k];
//				disp_state[k] = uj[k] - ui[k];
//				vel_state[k] = vj[k] - vi[k];
//				deform_state[k] = (xj[k] + uj[k]) - (xi[k] + ui[k]);
//			}
//			double r = sqrt(refpos_state[0]*refpos_state[0] + refpos_state[1]*refpos_state[1] + refpos_state[2]*refpos_state[2]);
//			double p = sqrt(deform_state[0]*deform_state[0] + deform_state[1]*deform_state[1] + deform_state[2]*deform_state[2]);
//			double exten_state = p - r;
//			double stretch_state = r>0. ? exten_state/r : 0.;
//			double vol_state = node_j->GetNodeVolume();
//			// break bonds based on bond stretch, this is same as "simple_dmg" subroutine of state-based CCI
//			broken_state_new = bond_ij->IsBreak_state();
//			double ecrit = mat_i->GetCriticalStretch(); // ??
//			if (stretch_state>ecrit && !broken_state_new)
//			{
//				bond_ij->SetBreak_state(true);
//			}
//			// find the influence function (either 0 or 1)
//			bond_ij->SetInfluence(0);
//			if (!bond_ij->IsBreak_state())
//			{
//				bond_ij->SetInfluence(1);
//			}
//			double influence_state = bond_ij->GetInfluence();
//			// find shape tensor K
//			if (!bond_ij->IsBreak_state())
//			{
//				for (int i=0;i!=3;++i)
//				{
//					for (int j=0;j!=3;++j)
//					{
//						t_shape[i][j] = refpos_state[i]*refpos_state[j]*vol_state*influence_state;
//						t_disp_grad[i][j] = disp_state[i]*refpos_state[j]*vol_state*influence_state;
//						t_vel_grad[i][j] = vel_state[i]*refpos_state[j]*vol_state*influence_state;
//						// sum up the tensors from each bond
//						shape[i][j] += t_shape[i][j];
//						sum_disp_grad[i][j] += t_disp_grad[i][j];
//						sum_vel_grad[i][j] += t_vel_grad[i][j];
//					}
//				}	
//			}
//	}
//	
//	// if there are too few unbroken bonds left to Get an accurate value of F, give up.
//	// otherwise, find the determinant of shape tensor
//	if(damage_state<=0.8 && nbond_not_broken>=10)
//	{
//		for (int i=0;i!=3;++i)
//		{
//			for (int j=0;j!=3;++j)
//			{
//				for (int k=0;k!=3;++k)
//				{
//					det_shape += pdIntgManager::CalcPermutationSymbol(i, j, k)*shape[0][i]*shape[1][j]*shape[2][k];
//				}
//			}
//		}
//		// if shape tensor is singular, give up. otherwise, find its inverse and gradient tensor
//		if (det_shape>0.)
//		{
//			// adjoint method
//			shape_inv[0][0] =  (shape[1][1]*shape[2][2] - shape[2][1]*shape[1][2])/det_shape;
//			shape_inv[0][1] = -(shape[0][1]*shape[2][2] - shape[2][1]*shape[0][2])/det_shape;
//			shape_inv[0][2] =  (shape[0][1]*shape[1][2] - shape[1][1]*shape[0][2])/det_shape;
//			shape_inv[1][0] = -(shape[1][0]*shape[2][2] - shape[2][0]*shape[1][2])/det_shape;
//			shape_inv[1][1] =  (shape[0][0]*shape[2][2] - shape[2][0]*shape[0][2])/det_shape;
//			shape_inv[1][2] = -(shape[0][0]*shape[1][2] - shape[1][0]*shape[0][2])/det_shape;
//			shape_inv[2][0] =  (shape[1][0]*shape[2][1] - shape[2][0]*shape[1][1])/det_shape;
//			shape_inv[2][1] = -(shape[0][0]*shape[2][1] - shape[2][0]*shape[0][1])/det_shape;
//			shape_inv[2][2] =  (shape[0][0]*shape[1][1] - shape[1][0]*shape[0][1])/det_shape;
//			// find the displacement gradient, velocity gradient and deformation gradient
//			for (int i=0;i!=3;++i)
//			{
//				for (int j=0;j!=3;++j)
//				{
//					for (int k=0;k!=3;++k)
//					{
//						disp_grad[i][j] += sum_disp_grad[i][k]*shape_inv[k][j];
//						vel_grad[i][j] += sum_vel_grad[i][k]*shape_inv[k][j];
//						def_grad[i][j] = disp_grad[i][j] + pdIntgManager::CalcKroneckerDelta(i, j);
//					}
//				}
//			}
//			// find the Jacobian (determinant of deformation gradient)
//			for (int i=0;i!=3;++i)
//			{
//				for (int j=0;j!=3;++j)
//				{
//					for (int k=0;k!=3;++k)
//					{
//						det_def_grad += pdIntgManager::CalcPermutationSymbol(i, j, k)*def_grad[0][i]*def_grad[1][j]*def_grad[2][k];
//					}
//				}
//			}
//			// if the deformation gradient tensor is singular, give up. otherwise, find its inverse
//			if (det_def_grad>0.)
//			{
//				def_gra_inv[0][0] =  (def_grad[1][1]*def_grad[2][2] - def_grad[2][1]*def_grad[1][2])/det_def_grad;
//				def_gra_inv[0][1] = -(def_grad[0][1]*def_grad[2][2] - def_grad[2][1]*def_grad[0][2])/det_def_grad;
//				def_gra_inv[0][2] =  (def_grad[0][1]*def_grad[1][2] - def_grad[1][1]*def_grad[0][2])/det_def_grad;
//				def_gra_inv[1][0] = -(def_grad[1][0]*def_grad[2][2] - def_grad[2][0]*def_grad[1][2])/det_def_grad;
//				def_gra_inv[1][1] =  (def_grad[0][0]*def_grad[2][2] - def_grad[2][0]*def_grad[0][2])/det_def_grad;
//				def_gra_inv[1][2] = -(def_grad[0][0]*def_grad[1][2] - def_grad[1][0]*def_grad[0][2])/det_def_grad;
//				def_gra_inv[2][0] =  (def_grad[1][0]*def_grad[2][1] - def_grad[2][0]*def_grad[1][1])/det_def_grad;
//				def_gra_inv[2][1] = -(def_grad[0][0]*def_grad[2][1] - def_grad[2][0]*def_grad[0][1])/det_def_grad;
//				def_gra_inv[2][2] =  (def_grad[0][0]*def_grad[1][1] - def_grad[1][0]*def_grad[0][1])/det_def_grad;
//				detf = det_def_grad;
//			}
//			else
//			{
//				for (int i=0;i!=3;++i)
//				{
//					for (int j=0;j!=3;++j)
//					{
//						disp_grad[i][j] = 0.;
//						vel_grad[i][j] = 0.;
//						def_grad[i][j] = pdIntgManager::CalcKroneckerDelta(i, j);
//					}
//				}
//			}
//		}
//	}
//
//	// ****** Step II ******
//	// Find the left polar decomposition of the deformation gradient tensor and the unrotated rate of deformation tensor. 
//	// The following procedure basically follows Flanagan & Taylor, Stress integration with finite rotations, 
//	// Computational Methods in Applied Mechanics and Engineering vol 62 (1987) pp 305-320
//	// but with one main differences: The left stretch tensor is found from V = F R' instead of integrating V over time.
//	// Find:
//	// Rate of deformation tensor = rate_of_def = D 
//	// Eulerian velocity gradient = vel_gra_eul = L
//	// Spin tensor = spin = W
//	// Rotation tensor = rot_tens = R
//	// Unrotated rate of deformation tensor = rate_of_def_unrot = d
//	// This step is same as "stretch_rates" subroutine of state-based CCI.
//
//    // ****** Step 2.1: find D and W ******
//	memset(vel_gra_eul, 0, sizeof(vel_gra_eul));
//	memset(rate_of_def, 0, sizeof(rate_of_def));
//	memset(spin, 0, sizeof(spin));
//	// find the Eulerian velocity gradient, L, from L = Fdot Finv
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			for (int k=0;k!=3;++k)
//			{
//				vel_gra_eul[i][j] += vel_grad[i][k]*def_gra_inv[k][j];
//			}
//		}
//	}	
//	// find the rate of deformation tensor, D, from D = sym L
//	// find the spin tensor, W, from W = skw L
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			rate_of_def[i][j] = 0.5*(vel_gra_eul[i][j] + vel_gra_eul[j][i]);
//			spin[i][j] = 0.5*(vel_gra_eul[i][j] - vel_gra_eul[j][i]);
//		}
//	}
//	
//	// ****** Step 2.2: find z, omega, and Omega ******
//	memset(dv, 0, sizeof(dv));
//	double trace_of_v = 0.;
//	double det_temp = 0.;
//	// find the product DV
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			for (int k=0;k!=3;++k)
//			{
//				dv[i][j] += rate_of_def[i][k]*(node_i->GetStress()->Get_left_stretch(k, j));
//			}
//		}
//	}
//	// find z, the vector dual of DV-VD = 2*dual(DV)
//	z[0] = dv[2][1] - dv[1][2];
//	z[1] = -dv[2][0] + dv[0][2];
//	z[2] = dv[1][0] - dv[0][1];
//	// find w, the vector dual of W
//	w[0] = 0.5*(spin[2][1] - spin[1][2]);
//	w[1] = 0.5*(-spin[2][0] + spin[0][2]);
//	w[2] = 0.5*(spin[1][0] - spin[0][1]);
//	// find omega vector (Flanagan & Taylor eqn 12)
//	for (int k=0;k!=3;++k)
//	{
//		trace_of_v += node_i->GetStress()->Get_left_stretch(k, k);
//	}
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			// I*tr(V) - V
//			temp[i][j] = trace_of_v*ident[i][j] - node_i->GetStress()->Get_left_stretch(i, j);
//		}
//	}
//	// find the determinant of temp tensor
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			for (int k=0;k!=3;++k)
//			{
//				det_temp += pdIntgManager::CalcPermutationSymbol(i, j, k)*temp[0][i]*temp[1][j]*temp[2][k];
//			}
//		}
//	}
//	// if it is not singular, find its inverse
//	if (det_temp!=0.)
//	{
//		// temp_inv = inv( I*tr(V) - V )
//		temp_inv[0][0] =  (temp[1][1]*temp[2][2] - temp[2][1]*temp[1][2])/det_temp;
//		temp_inv[0][1] = -(temp[0][1]*temp[2][2] - temp[2][1]*temp[0][2])/det_temp;
//		temp_inv[0][2] =  (temp[0][1]*temp[1][2] - temp[1][1]*temp[0][2])/det_temp;
//		temp_inv[1][0] = -(temp[1][0]*temp[2][2] - temp[2][0]*temp[1][2])/det_temp;
//		temp_inv[1][1] =  (temp[0][0]*temp[2][2] - temp[2][0]*temp[0][2])/det_temp;
//		temp_inv[1][2] = -(temp[0][0]*temp[1][2] - temp[1][0]*temp[0][2])/det_temp;
//		temp_inv[2][0] =  (temp[1][0]*temp[2][1] - temp[2][0]*temp[1][1])/det_temp;
//		temp_inv[2][1] = -(temp[0][0]*temp[2][1] - temp[2][0]*temp[0][1])/det_temp;
//		temp_inv[2][2] =  (temp[0][0]*temp[1][1] - temp[1][0]*temp[0][1])/det_temp;
//	}
//	else
//	{
//		if (_myRank==0)
//		{
//			outfile << "CalcStateBasedForce: singular temp_inv at node = " << node_i->GetID()+1 << endl;
//			cout << "CalcStateBasedForce: singular temp_inv at node = " << node_i->GetID()+1 << endl;
//		}
//		exit(0);
//	}
//	// find omega_vec = w + (inv(I*tr(V) - V)) z
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			omega_vec[i] = w[i] + temp_inv[i][j]*z[j];
//		}
//	}
//	// find Omega tensor from its vector dual, omega_vec
//	omega_tens[0][0] = 0.;
//	omega_tens[0][1] = -omega_vec[2];
//	omega_tens[0][2] = omega_vec[1];
//	omega_tens[1][0] = omega_vec[2];
//	omega_tens[1][1] = 0.;
//	omega_tens[1][2] = -omega_vec[0];
//	omega_tens[2][0] = -omega_vec[1];
//	omega_tens[2][1] = omega_vec[0];
//	omega_tens[2][2] = 0.;
//	
//	// ****** Step 2.3: update R ******
//	memset(scratch_matrix, 0, sizeof(scratch_matrix));
//	memset(temp, 0, sizeof(temp));
//	double det_tempb = 0.;
//	// increment the rotation tensor, R
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			// TempA = I + 0.5*dt*Omega, TempB = I - 0.5*dt*Omega
//			tempa[i][j] = ident[i][j] + 0.5*dt*omega_tens[i][j];
//			tempb[i][j] = ident[i][j] - 0.5*dt*omega_tens[i][j];
//		}
//	}
//	// find the determinant of tempb tensor
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			for (int k=0;k!=3;++k)
//			{
//				det_tempb += pdIntgManager::CalcPermutationSymbol(i, j, k)*tempb[0][i]*tempb[1][j]*tempb[2][k];
//			}
//		}
//	}
//	// if it is not singular, find its inverse
//	if (det_tempb!=0.)
//	{
//		tempb_inv[0][0] =  (tempb[1][1]*tempb[2][2] - tempb[2][1]*tempb[1][2])/det_tempb;
//		tempb_inv[0][1] = -(tempb[0][1]*tempb[2][2] - tempb[2][1]*tempb[0][2])/det_tempb;
//		tempb_inv[0][2] =  (tempb[0][1]*tempb[1][2] - tempb[1][1]*tempb[0][2])/det_tempb;
//		tempb_inv[1][0] = -(tempb[1][0]*tempb[2][2] - tempb[2][0]*tempb[1][2])/det_tempb;
//		tempb_inv[1][1] =  (tempb[0][0]*tempb[2][2] - tempb[2][0]*tempb[0][2])/det_tempb;
//		tempb_inv[1][2] = -(tempb[0][0]*tempb[1][2] - tempb[1][0]*tempb[0][2])/det_tempb;
//		tempb_inv[2][0] =  (tempb[1][0]*tempb[2][1] - tempb[2][0]*tempb[1][1])/det_tempb;
//		tempb_inv[2][1] = -(tempb[0][0]*tempb[2][1] - tempb[2][0]*tempb[0][1])/det_tempb;
//		tempb_inv[2][2] =  (tempb[0][0]*tempb[1][1] - tempb[1][0]*tempb[0][1])/det_tempb;
//	}
//	else
//	{
//		if (_myRank==0)
//		{
//			outfile << "CalcStateBasedForce: singular tempb_inv at node = " << node_i->GetID()+1 << endl;
//			cout << "CalcStateBasedForce: singular tempb_inv at node = " << node_i->GetID()+1 << endl;
//		}
//		exit(0);
//	}
//	// save the rotation tensor of previous step
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			rot_tens_old[i][j] = node_i->GetStress()->Get_rot_tens(i, j);
//		}
//	}
//	// find inv(tempB)*tempA
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			for (int k=0;k!=3;++k)
//			{
//				scratch_matrix[i][j] += tempb_inv[i][k]*tempa[k][j];
//			}
//		}
//	}
//	// update the new rotation tensor
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			for (int k=0;k!=3;++k)
//			{
//				
//				temp[i][j] += scratch_matrix[i][k]*rot_tens_old[k][j];
//			}
//			node_i->GetStress()->Set_rot_tens(i, j, temp[i][j]);
//		}
//	}
//	
//	// ****** Step 2.4: find d ******
//	memset(scratch_matrix, 0, sizeof(scratch_matrix));
//	memset(rate_of_def_unrot, 0, sizeof(rate_of_def_unrot));
//	// find the unrotated rate of deformation tensor, d, from d = R' D R
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			for (int k=0;k!=3;++k)
//			{
//				// D*R
//				scratch_matrix[i][j] += rate_of_def[i][k]*node_i->GetStress()->Get_rot_tens(k, j);
//			}
//		}
//	}
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			for (int k=0;k!=3;++k)
//			{
//				// R'*D*R
//				rate_of_def_unrot[i][j] += node_i->GetStress()->Get_rot_tens(k, i)*scratch_matrix[k][j];
//			}
//		}
//	}
//	
//	// ****** Step 2.5: find V ******
//	memset(tempa, 0, sizeof(tempa));
//	// find the left stretch tensor from V = F R'
//	// use tempa[3][3] as left_stretch buffer
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			for (int k=0;k!=3;++k)
//			{
//				tempa[i][j] += def_grad[i][k]*node_i->GetStress()->Get_rot_tens(j, k);
//			}
//			node_i->GetStress()->Set_left_stretch(i, j, tempa[i][j]);
//		}
//	}
//
//	// ****** Step III ******
//	// Find the new force state at node mi. The material is linear elastic ordinary state-based
//	// The constitutive model is applied in an unrotated configuration, then the resulting 
//	// stress tensor is rotated to the current configuration according to R, the rotation 
//	// tensor from the left polar decomposition (see Flanagan & Taylor, 1987)	
//	// This is same as "model100" subroutine of state-based CCI.
//	double young_mod = mat_i->GetYoungModulus();
//	double poisson = mat_i->GetNu();
//	double bulk_mod = mat_i->GetBulkModulus();
//	double shear_mod = mat_i->GetShearModulus();
//	memset(strain_inc_dev, 0, sizeof(strain_inc_dev));
//	memset(temp, 0, sizeof(temp));
//	memset(stress_unrot_old, 0, sizeof(stress_unrot_old));
//	memset(tempa, 0, sizeof(tempa));
//	double elastic = 0.;
//	// find the unrotated strain increment tensor from d = R'*D*R, the unrotated rate of deformation tensor
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			strain_inc[i][j] = dt*rate_of_def_unrot[i][j];
//		}
//	}
//	// find the dilatation increment
//	double dilat_inc = strain_inc[0][0] + strain_inc[1][1] + strain_inc[2][2];
//	// find the deviatoric strain increment
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			strain_inc_dev[i][j] = strain_inc[i][j] - ident[i][j]*dilat_inc/3.;
//		}
//	}
//    // find the old unrotated Cauchy stress tensor
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			for (int k=0;k!=3;++k)
//			{
//				// T R
//				temp[i][j] += (node_i->GetStress()->Get_stress_cauchy(i, k))*rot_tens_old[k][j];
//			}
//		}
//	}	
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			for (int k=0;k!=3;++k)
//			{
//				// sigma_old = R' T R
//				stress_unrot_old[i][j] += rot_tens_old[k][i]*temp[k][j];
//			}
//		}
//	}
//	// find the new unrotated Cauchy stress tensor
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			// sigma_new = sigma_old + 2G*eps_dev + K*tr(eps)*I
//			stress_unrot_new[i][j] = stress_unrot_old[i][j] + 2.*shear_mod*strain_inc_dev[i][j] + ident[i][j]*bulk_mod*dilat_inc;
//		}
//	}
//	// find the new rotated Cauchy stress tensor
//	memset(temp, 0, sizeof(temp));
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			for (int k=0;k!=3;++k)
//			{
//				// sigma_new R'
//				temp[i][j] += stress_unrot_new[i][k]*(node_i->GetStress()->Get_rot_tens(j, k));
//			}
//		}
//	}
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			for (int k=0;k!=3;++k)
//			{
//				// T_new = R sigma_new R'
//				tempa[i][j] += (node_i->GetStress()->Get_rot_tens(i, k))*temp[k][j];
//			}
//			node_i->GetStress()->Set_stress_cauchy(i, j, tempa[i][j]);
//		}
//	}
//	// find the elastic energy density
//	strain_unrot[0][0] = (stress_unrot_new[0][0] - poisson*(stress_unrot_new[1][1] + stress_unrot_new[2][2]))/young_mod;
//	strain_unrot[1][1] = (stress_unrot_new[1][1] - poisson*(stress_unrot_new[2][2] + stress_unrot_new[0][0]))/young_mod;
//	strain_unrot[2][2] = (stress_unrot_new[2][2] - poisson*(stress_unrot_new[0][0] + stress_unrot_new[1][1]))/young_mod;
//	strain_unrot[0][1] = stress_unrot_new[0][1]/(2.*shear_mod);
//	strain_unrot[1][2] = stress_unrot_new[1][2]/(2.*shear_mod);
//	strain_unrot[2][0] = stress_unrot_new[2][0]/(2.*shear_mod);
//	strain_unrot[1][0] = strain_unrot[0][1];
//	strain_unrot[2][1] = strain_unrot[1][2];
//	strain_unrot[0][2] = strain_unrot[2][0];
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			elastic += 0.5*strain_unrot[i][j]*stress_unrot_new[i][j];
//		}
//	}
//	// (kebing) starts
//	// store the elastic energy density
//	node_i->AddElastEnergyDensity(elastic);
//	// store the rotated strain, tempa is the buff of new rotated cauchy stress
//	memset(strain_cauchy, 0, sizeof(strain_cauchy));
//	strain_cauchy[0][0] = (tempa[0][0] - poisson*(tempa[1][1] + tempa[2][2]))/young_mod;
//	strain_cauchy[1][1] = (tempa[1][1] - poisson*(tempa[2][2] + tempa[0][0]))/young_mod;
//	strain_cauchy[2][2] = (tempa[2][2] - poisson*(tempa[0][0] + tempa[1][1]))/young_mod;
//	strain_cauchy[0][1] = tempa[0][1]/(2.*shear_mod);
//	strain_cauchy[1][2] = tempa[1][2]/(2.*shear_mod);
//	strain_cauchy[2][0] = tempa[2][0]/(2.*shear_mod);
//	strain_cauchy[1][0] = strain_cauchy[0][1];
//	strain_cauchy[2][1] = strain_cauchy[1][2];
//	strain_cauchy[0][2] = strain_cauchy[2][0];
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			node_i->GetStress()->Set_strain_cauchy(i, j, strain_cauchy[i][j]);
//		}
//	}
//	// (kebing) ends
//
//	// ****** Step IV ******
//	// Find the Piola stress tensor from the Cauchy stress tensor
//	// if the deformation gradient tensor is singular, give up. otherwise, find its inverse
//	// This is same as "cauchy_to_piola" subroutine of state-based CCI.
//	memset(stress_piola, 0, sizeof(stress_piola));
//	if (detf>0.)
//	{
//		for (int i=0;i!=3;++i)
//		{
//			for (int j=0;j!=3;++j)
//			{
//				for (int k=0;k!=3;++k)
//				{
//					stress_piola[i][j] += detf*node_i->GetStress()->Get_stress_cauchy(i, k)*def_gra_inv[j][k];
//				}
//				node_i->GetStress()->Set_stress_piola(i, j, stress_piola[i][j]);
//			}
//		}
//	}
//
//	// ****** Step V ******
//	// Find the force state at node mi from the expansion of the Piola stress tensor (see equation 142 of Silling et al, 2007)
//	// This is same as "piola_to_force_state" subroutine of state-based CCI.
//	memset(temp, 0, sizeof(temp));
//	for (int i=0;i!=3;++i)
//	{
//		for (int j=0;j!=3;++j)
//		{
//			for (int k=0;k!=3;++k)
//			{
//				temp[i][j] += stress_piola[i][k]*shape_inv[k][j];
//			}
//		}
//	}
//	// stress is reduced with increasing damage.
//	// dmg_zero_stress is the damage at which the stress becomes zero.
//	double dmg_zero_stress = 0.8;
//	double fac = 1. - pow(min(1., damage_state/dmg_zero_stress), 2.0);
//	// find the force state from eqn (142) of Silling et al, "Peridynamic states and constitutive modeling", J. Elast vol 88 pp 151-184.
//	for (int n=0;n!=numfam;++n)
//	{
//		pdNode* node_j = node_i->GetFamilyNode(n);
//		pdBond* bond_ij = node_i->GetFamilyBond(n);
//		double force_state_1 = 0.;
//		double force_state_2 = 0.;
//		double force_state_3 = 0.;
//			xj[0] = node_j->GetX1();
//			xj[1] = node_j->GetX2();
//			xj[2] = node_j->GetX3();
//			refpos_state[0] = xj[0] - xi[0];
//			refpos_state[1] = xj[1] - xi[1];
//			refpos_state[2] = xj[2] - xi[2];
//			double influence_state = bond_ij->GetInfluence();
//			force_state_1 = fac*influence_state*(temp[0][0]*refpos_state[0] + temp[0][1]*refpos_state[1] + temp[0][2]*refpos_state[2]);
//			force_state_2 = fac*influence_state*(temp[1][0]*refpos_state[0] + temp[1][1]*refpos_state[1] + temp[1][2]*refpos_state[2]);
//			force_state_3 = fac*influence_state*(temp[2][0]*refpos_state[0] + temp[2][1]*refpos_state[1] + temp[2][2]*refpos_state[2]);
//		// store the force value into family node array
//		bond_ij->SetForce_state_1(force_state_1);
//		bond_ij->SetForce_state_2(force_state_2);
//		bond_ij->SetForce_state_3(force_state_3);
//	}
//
//	outfile << endl;
//	outfile.close();
//}

void pdDataManager::CalcBondForceCCI(pdBond* bond_ij, pdNode* node_i, pdNode* node_j, double* pdForce, int& numIntg)
{
	// calculate bond-based peridynamic force using CCI method
	numIntg = 0; // total number of integration points
	double xi1 = node_j->GetX1() - node_i->GetX1();
	double xi2 = node_j->GetX2() - node_i->GetX2();
	double xi3 = node_j->GetX3() - node_i->GetX3();
	double eta1 = node_j->GetU1() - node_i->GetU1();
	double eta2 = node_j->GetU2() - node_i->GetU2();
	double eta3 = node_j->GetU3() - node_i->GetU3();
	double r = sqrt(xi1*xi1 + xi2*xi2 + xi3*xi3);
	double p = sqrt(pow(xi1+eta1, 2.0) + pow(xi2+eta2, 2.0) + pow(xi3+eta3, 2.0));
	pdMaterial* mat_i = node_i->GetMaterial();
	pdMaterial* mat_j = node_j->GetMaterial();
	double ecrit = min(mat_i->GetCriticalStretch(), mat_j->GetCriticalStretch());
	double spring = min(mat_i->GetSpringConstant(), mat_j->GetSpringConstant());
	double bdryEffij = 0.5*(node_i->GetBdryEffFactor() + node_j->GetBdryEffFactor());
	//1. exclude nodes outside the horizon
	//if (r <= _horizon)
	//2. no if test
	{
		// Get bond extension and bond stretch
		double u = p - r;
		double stretch = u / r; // current bond stretch
		//bond_ij->SetStretch(stretch);
		if (u <= 0.0 || (u > 0.0 && stretch <= ecrit)) // bond in compression or not failed and in tension
		{
			// calc volume reduction factor beta 
			double beta = this->CalcBetaFactor(r, _horizon, _gridSpacing);
			if (beta>0)
			{
				++numIntg;
			}
			double vol = node_j->GetNodeVolume();			
			pdForce[0] = bdryEffij * vol * beta * spring * stretch * (xi1 + eta1) / p;
			pdForce[1] = bdryEffij * vol * beta * spring * stretch * (xi2 + eta2) / p;
			pdForce[2] = bdryEffij * vol * beta * spring * stretch * (xi3 + eta3) / p;
			//// apply linear damping coefficient (relative velocity)
			//double dv1 = node_j->GetV1() - node_i->GetV1();
			//double dv2 = node_j->GetV2() - node_i->GetV2();
			//double dv3 = node_j->GetV3() - node_i->GetV3();
			//pdForce[0] = bdryEffij * vol * beta * spring * (stretch) * (xi1 + eta1) / p  + _dampCoeff * (dv1)/r;
			//pdForce[1] = bdryEffij * vol * beta * spring * (stretch) * (xi2 + eta2) / p  + _dampCoeff * (dv2)/r;
			//pdForce[2] = bdryEffij * vol * beta * spring * (stretch) * (xi3 + eta3) / p  + _dampCoeff * (dv3)/r;			
			// store the elastic energy density (each node shares half of the total energy of one bond)
			double elast = bdryEffij * 0.5 * (0.5 * beta * (spring/r) * u * u);
			node_i->AddElastEnergyDensity(elast);	
		}
		else // failed and in tension
		{
			pdForce[0] = 0.0;
			pdForce[1] = 0.0;
			pdForce[2] = 0.0;
			// update bond breakage status
			bond_ij->SetBreak(true);
			node_i->AddNumDmgBond(1);
		}
	}
	//else
	//{
	//	pdForce[0] = 0.0;
	//	pdForce[1] = 0.0;
	//	pdForce[2] = 0.0;
	//	// update bond breakage status
	//	bond_ij->SetBreak(true);
	//	node_i->AddNumDmgBond(1);
	//}	
}



void pdDataManager::CalcBondForceIntg(pdBond* bond_ij, pdNode* node_i, pdNode* node_j, double* pdForce, int& numIntg)
{
	/*
	
	This function uses fixed Gaussian integration or adaptive trapezoidal integration to carry out 
	the force integration over volume.

	Input:
	  1. bond_ij: pointer to the bond of interest.
	  2. node_i, node_j: pointer to two nodes of this bond.

	Output:
	  1. forceij_xi (i=1, 2, 3): three force components.
	  2. numIntg: total number of integration points used. 

	*/

	numIntg = 0; // total number of integration points
	pdIntgManager* im = new pdIntgManager(_horizon, _gridSpacing, bond_ij, node_i, node_j);
	im->SetIntgPoint(_intgPoint);
	im->SetTrapFlag(_aiOn);
	im->SetErrorControlFlag(_errorControlOn);
	im->SetErrorControlEPS(_eps);
	im->SetGaussianFlag(_fgiOn);
	if (_fgiOn)
	{
		im->SetGaussianAbscis(_abscis);
		im->SetGaussianWeight(_weights);
	}
	im->SetMLSFlag(_mlsOn);

	if (_horizon <= bond_ij->GetMinDistance()) // too far away (should never get here since we are doing calculation in node i's family)
	{
		pdForce[0] = 0.0;
		pdForce[1] = 0.0;
		pdForce[2] = 0.0;
	}
	else if (_horizon >= bond_ij->GetMaxDistance()) // cell of node j is fully in the horizon sphere
	{			
		im->CalcConfigurationType3("CalcBondBasedForce", pdForce);
		// integration points calculation
		if (_aiOn)
			numIntg += im->GetEndIntgPoint(1)*im->GetEndIntgPoint(2)*im->GetEndIntgPoint(3);
		else // fixed Gaussian integration
			numIntg += im->GetEndIntgPoint(0);
	}
	else // node j's cell has intersection with the horizon sphere
	{
		if (_aiOn)
		{
			// Adaptive trapezoidal integration
			im->AdaptiveIntegration("CalcBondBasedForce", pdForce);
			numIntg += im->GetEndIntgPoint(1)*im->GetEndIntgPoint(2)*im->GetEndIntgPoint(3);
		}
		else if (_fgiOn)
		{
			// Fixed Gaussian integration
			im->CalcConfigurationType3("CalcBondBasedForce", pdForce);
			numIntg += im->GetEndIntgPoint(0);
		}
	}
	delete im;
}

void pdDataManager::CalcBondForceIntg(pdBond* bond_ij, const double* u, pdNode* node_i, double** cbDispl, pdNode* node_j, const double time, double* pdForce, int& numIntg)
{
	/*
	
	This function uses fixed Gaussian integration or adaptive trapezoidal integration to carry out 
	the force integration over volume.

	Input:
	  1. bond_ij: pointer to the bond of interest.
	  2. node_i, node_j: pointer to two nodes of this bond.

	Output:
	  1. forceij_xi (i=1, 2, 3): three force components.
	  2. numIntg: total number of integration points used. 

	*/

	numIntg = 0; // total number of integration points
	pdIntgManager* im = new pdIntgManager(_horizon, _gridSpacing, bond_ij, u, node_i, node_j);
	im->SetIntgPoint(_intgPoint);
	im->SetTrapFlag(_aiOn);
	im->SetErrorControlFlag(_errorControlOn);
	im->SetErrorControlEPS(_eps);
	// enable AI method with MLS approximation
	im->SetMLSFlag(_mlsOn);
	im->SetCBNodeDispl(cbDispl);
	im->SetCurrentRunTime(time);
	im->SetGaussianFlag(_fgiOn);
	if (_fgiOn)
	{
		im->SetGaussianAbscis(_abscis);
		im->SetGaussianWeight(_weights);
	}

	if (_horizon <= bond_ij->GetMinDistance()) // too far away (should never get here since we are doing calculation in node i's family)
	{
		pdForce[0] = 0.0;
		pdForce[1] = 0.0;
		pdForce[2] = 0.0;
	}
	else if (_horizon >= bond_ij->GetMaxDistance()) // cell of node j is fully in the horizon sphere
	{			
		im->CalcConfigurationType3("CalcBondBasedForce", pdForce);
		// integration points calculation
		if (_aiOn)
		{
			numIntg += im->GetEndIntgPoint(1)*im->GetEndIntgPoint(2)*im->GetEndIntgPoint(3);
		}
		else // fixed Gaussian integration
		{
			numIntg += im->GetEndIntgPoint(0);
		}
	}
	else // node j's cell has intersection with the horizon sphere
	{
		if (_aiOn)
		{
			//// Adaptive trapezoidal integration
			im->AdaptiveIntegration("CalcBondBasedForce", pdForce);
			numIntg += im->GetEndIntgPoint(1)*im->GetEndIntgPoint(2)*im->GetEndIntgPoint(3);
		}
		else if (_fgiOn)
		{
			// Fixed Gaussian integration
			im->CalcConfigurationType3("CalcBondBasedForce", pdForce);
			numIntg += im->GetEndIntgPoint(0);
		}
	}
	delete im;	
}

void pdDataManager::CalcBondForceCCI(const double* xi, const double* eta, const double ecrit, const double spring, 
									   const double bdryEffij, const double effWeight, double* pdForce, double& elast)
{	
    // This function uses CCI method to calculate peridynamic force.
	// the input 'effWeight' is actually a concept from Gaussian integration.
	// since CCI is an one-point integration, it is essentially same as 1*1*1 Gaussian.
	// therefore, for CCI method, the effWeight is always equal to 2. 
	// this function is also used by fixed Gaussian integration method when a direct summation scheme is used.

	double xi1 = xi[0];
	double xi2 = xi[1];
	double xi3 = xi[2];
	double eta1 = eta[0];
	double eta2 = eta[1];
	double eta3 = eta[2];
    // get |xi| and |xi+eta|
    double r = sqrt(xi1*xi1 + xi2*xi2 + xi3*xi3);
    double p = sqrt(pow(xi1+eta1, 2.0) + pow(xi2+eta2, 2.0) + pow(xi3+eta3, 2.0));
	// get scale factor beta	
	double dx = effWeight*_gridSpacing / 2.0; // get the effective subcell size of the Gaussian point
	double beta = CalcBetaFactor(r, _horizon, dx);
	/*if (beta > 0.0)
	{
		++_intgPointEnd[0];
	}*/

	// calculate the peridynamic force
    // Get bond extension and bond stretch
    double u = p - r;
    double stretch = u / r; // current bond stretch
	//_bond->SetStretch(stretch);
    if (u <= 0.0 || (u > 0.0 && stretch <= ecrit)) // bond in compression or not failed and in tension
    {
		if (_cciOn)
		{
			double vol = _gridSpacing*_gridSpacing*_gridSpacing; // currently use global grid spacing		
			pdForce[0] = bdryEffij * vol * beta * spring * stretch * (xi1 + eta1) / p;
			pdForce[1] = bdryEffij * vol * beta * spring * stretch * (xi2 + eta2) / p;
			pdForce[2] = bdryEffij * vol * beta * spring * stretch * (xi3 + eta3) / p;
		}
		else
		{
			// volume is multiplied outside the function
			pdForce[0] = bdryEffij * beta * spring * stretch * (xi1 + eta1) / p;
			pdForce[1] = bdryEffij * beta * spring * stretch * (xi2 + eta2) / p;
			pdForce[2] = bdryEffij * beta * spring * stretch * (xi3 + eta3) / p;
		}
		// save the elastic energy density (each node shares half of the total energy of one bond)
		// the energy need to be scaled based on the effective weight
		elast = bdryEffij * 0.5 * (0.5 * beta * (spring/r) * u * u) * pow(0.5*effWeight, 3);
    }
    else // failed and in tension
    {
        pdForce[0] = 0.0;
        pdForce[1] = 0.0;
        pdForce[2] = 0.0;
        elast = 0.0;
    }
}

void pdDataManager::CalcBdryEffFactor(pdNode* node_i)
{
	// This function calculates boundary effect compensation factor for the source node node_i. 
	// It must be called after function InitMaterial().
	
	double e11 = 1.0e-4;
	pdMaterial* mat_i = node_i->GetMaterial();
	double x1i = node_i->GetX1();
	double x2i = node_i->GetX2();
	double x3i = node_i->GetX3();
	double ecrit = mat_i->GetCriticalStretch();
	double spring = mat_i->GetSpringConstant();
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
	// create a temporary boundary condition pointer
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
	double elasticEnergyDensity = 0.0;
	int numfam = node_i->GetNumFamily();
	for (int n=0;n!=numfam;++n)
	{
		pdNode* node_j = node_i->GetFamilyNode(n);
		pdBond* bond = node_i->GetFamilyBond(n);		
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
			//if (_cciOn)
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
				if (u<=0 || (u>0 && stretch <= ecrit)) // bond in compression or not failed and in tension
					elastic = 0.5*(spring/r)*u*u;
				else // failed and in tension
					elastic = 0.0;
				// calc volume reduction factor beta
				double beta = this->CalcBetaFactor(r, _horizon, _gridSpacing);		
				elastic_j = elastic * beta * node_j->GetNodeVolume();
				//cout << r << "\t" << bond->GetMinDistance() << "\t" << beta  
				//	<< "\t" << elastic << "\t" << elastic_j << "\n";
			}
			//else
			//{
			//	// backup the initial boundary condition of node j
			//	pdBdryCondition* bcBuffer = node_j->GetBC();
			//	// use the temporary boundary condition
			//	node_j->SetBC(bcTemp);
			//	// use above given strain and 3*3*3 Gaussian to calc the displacement on the integration points
			//	// The 3*3*3 Gaussian points have been set up in function pdDataManager::InitBdryEffectFactor
			//	pdIntgManager* im = new pdIntgManager(_horizon, _gridSpacing, bond, node_i, node_j);
			//	im->SetIntgPoint(3);
			//	im->SetTrapFlag(false);
			//  im->SetErrorControlFlag(false);
			//  im->SetErrorControlEPS(0.0);
			//	im->SetGaussianFlag(true);
			//	im->SetGaussianAbscis(_abscis);
			//	im->SetGaussianWeight(_weights);
			//	im->SetMLSFlag(false);

			//	if (_horizon <= bond->GetMinDistance()) // too far away (should never get here)
			//	{
			//		elastic_j = 0.0;
			//	}
			//	else if (_horizon >= bond->GetMaxDistance()) // cubic of node j is fully in the horizon of node i
			//	{
			//		elastic_j = im->CalcConfigurationType3("GetStrainEnergyDensity");
			//	}
			//	else
			//	{
			//		//if (_fgiOn)
			//			elastic_j = im->CalcConfigurationType3("GetStrainEnergyDensity");
			//		//else if (_aiOn)
			//		//	elastic_j = im->AdaptiveIntegration("GetStrainEnergyDensity");
			//	}
			//	delete im;
			//	// restore the initial boundary condition of node j
			//	node_j->SetBC(bcBuffer);
			//}
			// sum up the elastic energy density 
			elasticEnergyDensity += 0.5 * elastic_j;
		}
	}
	delete bcTemp;
	double p_grid = elasticEnergyDensity;
	double p_inf = 3.0*mat_i->GetYoungModulus()*e11*e11;
	node_i->SetBdryEffFactor(p_inf/p_grid);	
}

double pdDataManager::CalcBdryEffFactorTypeOne(double h, pdMaterial* mat)
{
	double fullEnergy = 3.0*mat->GetYoungModulus();
	double partialEnergy = fullEnergy -
	0.5*mat->GetSpringConstant()*PI*(pow(_horizon, 4)/4.0 - h*pow(_horizon, 3)/3.0 + pow(h, 4)/12.0);
	return (fullEnergy/partialEnergy);
}

//void pdDataManager::CalcForceThroughPlane(pdBond* bond_ij, const double f[], const double vol_sq)
//{
//	// This function loops through all planes this bond passes through and add average force to each plane
//	// input: f[0], f[1], f[2] is the bond force on this bond in three directions
//	// input: vol_sq is the input of volume_node_i*boolean_volume_node_j
//
//	pdNode* node_i = bond_ij->GetNodeI();
//	pdNode* node_j = bond_ij->GetNodeJ();
//	int num_plane = bond_ij->GetNumberPlane();
//
//	for (int i=0;i!=num_plane;++i)
//	{
//		pdSpacePlane* plane = bond_ij->GetPlane(i);
//		double x1lo = plane->GetX1Low();
//		double x1hi = plane->GetX1High();
//		double x2lo = plane->GetX2Low();
//		double x2hi = plane->GetX2High();
//		double x3lo = plane->GetX3Low();
//		double x3hi = plane->GetX3High();
//		int dir = abs(plane->GetNorm());
//		int sign_norm = plane->GetNorm()/dir;
//		double xi, xj, xplane;
//		if (dir==1)
//		{
//			xi = node_i->GetX1();
//			xj = node_j->GetX1();
//			xplane = x1lo;
//		}
//		else if (dir==2)
//		{
//			xi = node_i->GetX2();
//			xj = node_j->GetX2();
//			xplane = x2lo;
//		}
//		else
//		{
//			xi = node_i->GetX3();
//			xj = node_j->GetX3();
//			xplane = x3lo;
//		}
//		double sign_force = sign_norm*(xj-xplane)/abs(xj-xplane);
//		// check if node i and j are not at the same side of the plane
//		if (dir==1)
//		{
//			plane->AddForce(sign_force*f[0]*vol_sq/num_plane);
//		//	plane->AddForce(sign_force*f[0]*vol_sq);
//		}
//		else if (dir==2)
//		{
//			plane->AddForce(sign_force*f[1]*vol_sq/num_plane);
//		//	plane->AddForce(sign_force*f[1]*vol_sq);
//		}
//		else
//		{
//			plane->AddForce(sign_force*f[2]*vol_sq/num_plane);
//		//	plane->AddForce(sign_force*f[2]*vol_sq);
//		}
//	}
//}

pdVector* pdDataManager::GetBasis(double x, double y, double z, int dim)
{
    // calculate the polynomial basis of order 2 and store in basisArray array (dimension=10) 
	// x, y, z are the coordinates of one point
	pdVector* vec = new pdVector(dim);
    vec->SetCoeff(0, 1.0);
    vec->SetCoeff(1, x);
    vec->SetCoeff(2, y);
    vec->SetCoeff(3, z);
    vec->SetCoeff(4, x * y);
    vec->SetCoeff(5, x * z);
    vec->SetCoeff(6, y * z);
    vec->SetCoeff(7, x * x);
    vec->SetCoeff(8, y * y);
    vec->SetCoeff(9, z * z);
    return vec;
}

double pdDataManager::GetWeight(pdNode* n, double x, double y, double z, double rw)
{
	// n is the one of the contributing node
	// x, y, z are the coordinates of integration point
	double xi1 = x - n->GetX1();
	double xi2 = y - n->GetX2();
	double xi3 = z - n->GetX3();
	double r = sqrt(xi1*xi1 + xi2*xi2 + xi3*xi3);
    double d = r / rw;
    if (d > 1.0)
        return 0.0;
    else
        return (1 - 6.0 * d * d + 8.0 * pow(d, 3.0) - 3.0 * pow(d, 4.0)); // quartic spline weight function
}

bool pdDataManager::ConsistencyCheck(double x, double y, double z, vector<pdNode*> cbList, pdVector* coeff)
{
	// unity consistency check
	double uniSum = 0.0;
	for (int i = 0; i != coeff->GetNumRows(); ++i)
	    uniSum += coeff->GetCoeff(i);
	// approx. order consistency check
	double appOrder[BASISDIM-1];
	memset(appOrder, 0, sizeof(appOrder));
	for (int i = 0; i != cbList.size(); ++i)
	{
	    // linear order consistency
	    appOrder[0] += (cbList[i]->GetX1() - x) * coeff->GetCoeff(i);
	    appOrder[1] += (cbList[i]->GetX2() - y) * coeff->GetCoeff(i);
	    appOrder[2] += (cbList[i]->GetX3() - z) * coeff->GetCoeff(i);
	    // second order consistency
	    appOrder[3] += (cbList[i]->GetX1() * cbList[i]->GetX1() - x * x) * coeff->GetCoeff(i);
	    appOrder[4] += (cbList[i]->GetX2() * cbList[i]->GetX2() - y * y) * coeff->GetCoeff(i);
	    appOrder[5] += (cbList[i]->GetX3() * cbList[i]->GetX3() - z * z) * coeff->GetCoeff(i);
	    appOrder[6] += (cbList[i]->GetX1() * cbList[i]->GetX2() - x * y) * coeff->GetCoeff(i);
	    appOrder[7] += (cbList[i]->GetX2() * cbList[i]->GetX3() - y * z) * coeff->GetCoeff(i);
	    appOrder[8] += (cbList[i]->GetX3() * cbList[i]->GetX1() - z * x) * coeff->GetCoeff(i);
	}
	double firstNorm = 0.0;
	double secondNorm = 0.0;
	for (int i = 0; i != 3; ++i)
	    firstNorm += appOrder[i] * appOrder[i];
	for (int i = 0; i != POLYBASISDIM-1; ++i)
	    secondNorm += appOrder[i] * appOrder[i];
	double tol = 1.0e-8;
	if (abs(uniSum-1.0) > tol || firstNorm > tol || secondNorm > tol)
		return false;
	else
		return true;
}

double pdDataManager::CalcDampingCoeff(const double currentTimeStep, pdVector* absKVector)
{
	/*

	 Notice:
	   1. This function calculates the damping coefficient at the current time step for dynamic relaxation.
	   2. The force function is linearized so that a linear relationship of f with relative displacement 
	      can be acquired.
	   3. The stiffness matrix and fictitious diagonal mass matrices are not actually created since 
	      large number of nodes will require a huge amount of memory to construct K and M matrix, 
	      which will violate the size limit of array in C++ (2G bytes for 32-bit system).

	 Input: 
	   1. currentTimeStep: current time step length.

	 Output:
	   1. dampCoeff: the damping coefficient at current time step.
	   2. absKVector: the diagonal mass matrix. To save memory, only the diagonal terms are stored. 
	*/

	/*
    //The following code shows how the global stiffness matrix and mass matrix are constructed.

	int nrow = _numNodes*3; // dimension of stiffness and mass matrix (nrow by nrow square)
	pdVector* uVector = new pdVector(nrow);
	pdMatrix* kMatrix = new pdMatrix(nrow, nrow);
	pdMatrix* mMatrix = new pdMatrix(nrow, nrow);
	// loop over all nodes
	for (int i=0;i!=_numNodes;++i)
	{
		pdNode* node_i = _nodes[i];
		// Get the row number of node i in three dofs
		int irowx = _rowOfNode[i];
		int irowy = irowx + 1;
		int irowz = irowx + 2;
		// Get the node coordinate
		double x1i = node_i->GetX1();
		double x2i = node_i->GetX2();
		double x3i = node_i->GetX3();
		pdMaterial* mat_i = node_i->GetMaterial();
		// loop over the family nodes of node i
		int numfam = node_i->GetNumFamily();
		for (int n=0;n!=numfam;++n)
		{
			pdNode* node_j = node_i->GetFamilyNode(n);
			int j = node_j->GetID();
			// Get the row number of node j in three dofs
			int jrowx = _rowOfNode[j];
			int jrowy = jrowx + 1;
			int jrowz = jrowx + 2;
			// Get node coordinate
			double x1j = node_j->GetX1();
			double x2j = node_j->GetX2();
			double x3j = node_j->GetX3();
			// Get xi and its magnitude
			double xi1 = x1j - x1i;
			double xi2 = x2j - x2i;
			double xi3 = x3j - x3i;
			double r = sqrt(xi1*xi1 + xi2*xi2 + xi3*xi3);
			// Get the component of c
			pdMaterial* mat_j = node_j->GetMaterial();
			double spring = min(mat_i->GetSpringConstant(), mat_j->GetSpringConstant());
			double c11 = spring*xi1*xi1/pow(r, 3.0);
			double c12 = spring*xi1*xi2/pow(r, 3.0);
			double c13 = spring*xi1*xi3/pow(r, 3.0);
			double c21 = spring*xi2*xi1/pow(r, 3.0);
			double c22 = spring*xi2*xi2/pow(r, 3.0);
			double c23 = spring*xi2*xi3/pow(r, 3.0);
			double c31 = spring*xi3*xi1/pow(r, 3.0);
			double c32 = spring*xi3*xi2/pow(r, 3.0);
			double c33 = spring*xi3*xi3/pow(r, 3.0);
			double volj = node_j->GetNodeVolume();
			// fill the stiffness matrix
			// *** degree of freedom in x1 direction ***
			kMatrix->AddCoeff(irowx, irowx,  c11*volj);
			kMatrix->AddCoeff(irowx, jrowx, -c11*volj);
			kMatrix->AddCoeff(irowx, irowy, c12*volj);
			kMatrix->AddCoeff(irowx, jrowy, -c12*volj);
			kMatrix->AddCoeff(irowx, irowz, c13*volj);
			kMatrix->AddCoeff(irowx, jrowz, -c13*volj);
			// *** degree of freedom in x2 direction ***
			kMatrix->AddCoeff(irowy, irowx, c21*volj);
			kMatrix->AddCoeff(irowy, jrowx, -c21*volj);
			kMatrix->AddCoeff(irowy, irowy, c22*volj);
			kMatrix->AddCoeff(irowy, jrowy, -c22*volj);
			kMatrix->AddCoeff(irowy, irowz, c23*volj);
			kMatrix->AddCoeff(irowy, jrowz, -c23*volj);
			// *** degree of freedom in x3 direction ***
			kMatrix->AddCoeff(irowz, irowx, c31*volj);
			kMatrix->AddCoeff(irowz, jrowx, -c31*volj);
			kMatrix->AddCoeff(irowz, irowy, c32*volj);
			kMatrix->AddCoeff(irowz, jrowy, -c32*volj);
			kMatrix->AddCoeff(irowz, irowz, c33*volj);
			kMatrix->AddCoeff(irowz, jrowz, -c33*volj);
		}
		// fill the global displacement vector
		double u1i = node_i->GetU1();
		double u2i = node_i->GetU2();
		double u3i = node_i->GetU3();
		uVector->SetCoeff(irowx, u1i);
		uVector->SetCoeff(irowy, u2i);
		uVector->SetCoeff(irowz, u3i);
	}
	double nom = (kMatrix->Mult(uVector))->Mult(uVector); // {u}^T*[K]*{u}

	// assemble mass matrix
	for (int i=0;i!=nrow;++i)
	{
		double kRowSum =0.0;
		for (int j=0;j!=nrow;++j)
		{
			kRowSum += abs(kMatrix->GetCoeff(i, j));
		}
		mMatrix->SetCoeff(i, i, 0.25*pow(currentTimeStep, 2)*kRowSum);
	}
	double dem = (mMatrix->Mult(uVector))->Mult(uVector); // {u}^T*[M]*{u}

	*/
	
	// However, large grid will encounter the problem of array size exceeds limit of 2G bytes 
	// with above method. Therefore, the K and M matrices are not actually created.

	double dampCoeff; // the damping coefficient at current time step
	int nrow = _numNodes*3; // dimension of stiffness and mass matrix (nrow by nrow square)
	double dtSquare = pow(currentTimeStep, 2);
	pdVector* uVector = new pdVector(nrow);
	pdVector* kuVector = new pdVector(nrow); // [K]*{u}
	pdVector* muVector = new pdVector(nrow); // [M]*{u}
	pdMatrix* diagKSumMatrix = new pdMatrix(9, _numNodes); // the diagonal terms of |K_{ij}|
	// first loop over nodes to set up the global displacement vector
	for (int i=0;i!=_numNodes;++i)
	{
		pdNode* node_i = _nodes[i];
		// Get the row number of node i in three dof
		int irowx = _rowOfNode[i];
		int irowy = irowx + 1;
		int irowz = irowx + 2;
		uVector->SetCoeff(irowx, node_i->GetU1());
		uVector->SetCoeff(irowy, node_i->GetU2());
		uVector->SetCoeff(irowz, node_i->GetU3());
	}
	// zero other vectors and diagKSumMatrix
	kuVector->Zero();
	absKVector->Zero();// sum of |K_{ij}|
	muVector->Zero();
	diagKSumMatrix->Zero();

	// loop over nodes
	for (int i=0;i!=_numNodes;++i)
	{
		pdNode* node_i = _nodes[i];
		// Get the row number of node i in three dofs
		int irowx = _rowOfNode[i];
		int irowy = irowx + 1;
		int irowz = irowx + 2;
		// Get the node coordinate
		double x1i = node_i->GetX1();
		double x2i = node_i->GetX2();
		double x3i = node_i->GetX3();
		pdMaterial* mat_i = node_i->GetMaterial();
		// loop over the family nodes of node i
		int numfam = node_i->GetNumFamily();
		for (int n=0;n!=numfam;++n)
		{
			pdNode* node_j = node_i->GetFamilyNode(n);
			// Get the row number of node j in three dof
			int jrowx = _rowOfNode[node_j->GetID()];
			int jrowy = jrowx + 1;
			int jrowz = jrowx + 2;
			// Get node coordinate
			double x1j = node_j->GetX1();
			double x2j = node_j->GetX2();
			double x3j = node_j->GetX3();
			// Get xi and its magnitude
			double xi1 = x1j - x1i;
			double xi2 = x2j - x2i;
			double xi3 = x3j - x3i;
			double r = sqrt(xi1*xi1 + xi2*xi2 + xi3*xi3);
			// Get the component of c
			pdMaterial* mat_j = node_j->GetMaterial();
			double spring = min(mat_i->GetSpringConstant(), mat_j->GetSpringConstant());
			double c11 = spring*xi1*xi1/pow(r, 3.0);
			double c12 = spring*xi1*xi2/pow(r, 3.0);
			double c13 = spring*xi1*xi3/pow(r, 3.0);
			double c21 = spring*xi2*xi1/pow(r, 3.0);
			double c22 = spring*xi2*xi2/pow(r, 3.0);
			double c23 = spring*xi2*xi3/pow(r, 3.0);
			double c31 = spring*xi3*xi1/pow(r, 3.0);
			double c32 = spring*xi3*xi2/pow(r, 3.0);
			double c33 = spring*xi3*xi3/pow(r, 3.0);
			double volj = node_j->GetNodeVolume();
			// fill kuVector ([K]*{u})
			// *** degree of freedom in X1 direction ***
			kuVector->AddCoeff(irowx, c11*volj*uVector->GetCoeff(irowx));
			kuVector->AddCoeff(irowx, -c11*volj*uVector->GetCoeff(jrowx));
			kuVector->AddCoeff(irowx, c12*volj*uVector->GetCoeff(irowy));
			kuVector->AddCoeff(irowx, -c12*volj*uVector->GetCoeff(jrowy));
			kuVector->AddCoeff(irowx, c13*volj*uVector->GetCoeff(irowz));
			kuVector->AddCoeff(irowx, -c13*volj*uVector->GetCoeff(jrowz));
			// *** degree of freedom in X2 direction ***
			kuVector->AddCoeff(irowy, c21*volj*uVector->GetCoeff(irowx));
			kuVector->AddCoeff(irowy, -c21*volj*uVector->GetCoeff(jrowx));
			kuVector->AddCoeff(irowy, c22*volj*uVector->GetCoeff(irowy));
			kuVector->AddCoeff(irowy, -c22*volj*uVector->GetCoeff(jrowy));
			kuVector->AddCoeff(irowy, c23*volj*uVector->GetCoeff(irowz));
			kuVector->AddCoeff(irowy, -c23*volj*uVector->GetCoeff(jrowz));
			// *** degree of freedom in X3 direction ***
			kuVector->AddCoeff(irowz, c31*volj*uVector->GetCoeff(irowx));
			kuVector->AddCoeff(irowz, -c31*volj*uVector->GetCoeff(jrowx));
			kuVector->AddCoeff(irowz, c32*volj*uVector->GetCoeff(irowy));
			kuVector->AddCoeff(irowz, -c32*volj*uVector->GetCoeff(jrowy));
			kuVector->AddCoeff(irowz, c33*volj*uVector->GetCoeff(irowz));
			kuVector->AddCoeff(irowz, -c33*volj*uVector->GetCoeff(jrowz));
			// fill absKVector and diagKSumMatrix
			// *** degree of freedom in X1 direction ***
			diagKSumMatrix->AddCoeff(0, i, c11*volj);    //(irowx, irowx) in the stiffness matrix
			absKVector->AddCoeff(irowx, abs(-c11*volj)); // (irowx, jrowx)
			diagKSumMatrix->AddCoeff(1, i, c12*volj);    // (irowx, irowy)
			absKVector->AddCoeff(irowx, abs(-c12*volj)); // (irowx, jrowy)
			diagKSumMatrix->AddCoeff(2, i, c13*volj);    // (irowx, irowz)
			absKVector->AddCoeff(irowx, abs(-c13*volj)); // (irowx, jrowz)
			// *** degree of freedom in X2 direction ***
			diagKSumMatrix->AddCoeff(3, i, c21*volj);    // (irowy, irowx)
			absKVector->AddCoeff(irowy, abs(-c21*volj)); // (irowy, jrowx)
			diagKSumMatrix->AddCoeff(4, i, c22*volj);    // (irowy, irowy)
			absKVector->AddCoeff(irowy, abs(-c22*volj)); // (irowy, jrowy)
			diagKSumMatrix->AddCoeff(5, i, c23*volj);    // (irowy, irowz)
			absKVector->AddCoeff(irowy, abs(-c23*volj)); // (irowy, jrowz)
			// *** degree of freedom in X3 direction ***
			diagKSumMatrix->AddCoeff(6, i, c31*volj);    // (irowz, irowx)
			absKVector->AddCoeff(irowz, abs(-c31*volj)); // (irowz, jrowx)
			diagKSumMatrix->AddCoeff(7, i, c32*volj);    // (irowz, irowy)
			absKVector->AddCoeff(irowz, abs(-c32*volj)); // (irowz, jrowy)
			diagKSumMatrix->AddCoeff(8, i, c33*volj);    // (irowz, irowz)
			absKVector->AddCoeff(irowz, abs(-c33*volj)); // (irowz, jrowz)	
		}		
	}
	// fill the rest of absKVector
	for (int i=0;i!=_numNodes;++i)
	{
		pdNode* node_i = _nodes[i];
		// Get the row number of node i in three dof
		int irowx = _rowOfNode[i];
		int irowy = irowx + 1;
		int irowz = irowx + 2;
		// X1
		absKVector->AddCoeff(irowx, abs(diagKSumMatrix->GetCoeff(0, i))
			+abs(diagKSumMatrix->GetCoeff(1, i))+abs(diagKSumMatrix->GetCoeff(2, i)));
		// X2
		absKVector->AddCoeff(irowy, abs(diagKSumMatrix->GetCoeff(3, i))
			+abs(diagKSumMatrix->GetCoeff(4, i))+abs(diagKSumMatrix->GetCoeff(5, i)));
		// X3
		absKVector->AddCoeff(irowz, abs(diagKSumMatrix->GetCoeff(6, i))
			+abs(diagKSumMatrix->GetCoeff(7, i))+abs(diagKSumMatrix->GetCoeff(8, i)));
	}
	// right now, absKVecotor[i] = sum of |K_{ij}|

	// assemble muVector [M]*{u}
	for (int i=0;i!=nrow;++i)
	{
		muVector->SetCoeff(i, 0.25*dtSquare*absKVector->GetCoeff(i)*uVector->GetCoeff(i));
	}
	absKVector->Mult(0.25*dtSquare); // times dt^2*(1/4) to get the diagonal mass matrix
	double nom = uVector->Mult(kuVector); // {u}^T*[K]*{u}
	double dem = uVector->Mult(muVector); // {u}^T*[M]*{u}
	if (dem == 0.0)
	{
		dampCoeff = 0.0;
	}
	else if (nom/dem < 0.0)
	{
		dampCoeff = 0.0;
	}
	else
	{
		double omega = sqrt(nom/dem); // lowest frequency of the system
		dampCoeff = 2*omega;
	}
	
	//// debug
	//if (_myRank==0)
	//{
	//	string ofpath = _workDir + "pds_debug.dat";
	//	ofstream fout(ofpath.c_str(), ios::app);
	//	fout << nom << "\t" << dem << "\t" << omega << "\t" << dampCoeff << endl;
	//}

	delete kuVector;
	delete muVector;
	delete uVector;
	delete diagKSumMatrix;
	return dampCoeff;
}


void pdDataManager::CalcVolumeToPointForceDensity(const pdNode* node_i, double * pdForce)
{
	/* 
	
	This function calculates the volume-to-point force density for node 'node_i'.

	Input: 
	  1. 

	*/
}