#include "pdIntgManager.h"

void pdIntgManager::SetTrapFlag(const bool flag)
{
	_aiOn = flag;
}

void pdIntgManager::SetErrorControlFlag(const bool flag)
{
	_errorControlOn = flag;
}

void pdIntgManager::SetErrorControlEPS(const double eps)
{
	_eps = eps;
}

void pdIntgManager::SetGaussianFlag(const bool flag)
{
	_fgiOn = flag;
}

void pdIntgManager::SetGaussianAbscis(const vector<double>& abscis)
{
	_abscis.resize(_intgPoint);
	copy(abscis.begin(), abscis.end(), _abscis.begin());
}

void pdIntgManager::SetGaussianWeight(const vector<double>& weights)
{
	_weights.resize(_intgPoint);
	copy(weights.begin(), weights.end(), _weights.begin());
}

void pdIntgManager::SetMLSFlag(const bool flag)
{
	_mlsOn = flag;
}

void pdIntgManager::SetIntgPoint(const int p)
{
	_intgPoint = p;
}



void pdIntgManager::SetInitMLSCoeffFlag(const bool flag)
{
	_initMLSCoeffOn = flag;
}

void pdIntgManager::UnitValue(double x, double y, double z, double* unitValue)
{
	if (_fgiOn)
	{
		double dist = sqrt(x*x+y*y+z*z); //family coordinates
		// scale factor beta based on the location of this Gaussian point
		pdGauss* gPoint = _node_j->GetGaussianPoint(_intgPointCounter);
		double w1 = gPoint->GetWeightX1();
		double w2 = gPoint->GetWeightX2();
		double w3 = gPoint->GetWeightX3();
		double w_eff = pow(w1*w2*w3, 1.0/3.0); // geometric average
		double dx = w_eff*_dx/2.0; // effective cell size of this Gaussian point
		double beta = CalcBetaFactor(dist, _delta, dx);
		if (beta>0.0)
		{
			++_intgPointEnd[0];
		}
		unitValue[0] = beta;
		unitValue[1] = beta;
		unitValue[2] = beta;
	}
	else if (_aiOn)
	{
		unitValue[0] = 1.0;
		unitValue[1] = 1.0;
		unitValue[2] = 1.0;
		if (_initMLSCoeffOn)
		{
			// set up MLS coefficient list for this trapezoidal points
			this->InitMLSCoeff(x, y, z, _node_j);
		}
	}
	else
	{
		unitValue[0] = 1.0;
		unitValue[1] = 1.0;
		unitValue[2] = 1.0;
	}
}

void pdIntgManager::CalcBondBasedForce(const double x1j, const double x2j, const double x3j, double* pdForce)
{	
    /* 
	
	This function calculate bond-based force and return all three force components.
	  1. Note the given coordinate (x1j, x2j, x3j) of integration point is in the local coordinate.
      2. The displacements of the trapezoidal points are updated by:
	     1) If _mlsOn=true and the node where this trapezoidal point locates is applied with boundary condition,
			use the boundary condition value to update.
		 2) If _mlsOn=true and the node where this trapezoidal point locates is not applied with boundary condition,
		    use MLS approximation with the trial nodal displacements of RK4 steps to update.
	     3) If _mlsOn=false, this node must be applied with boundary condition and use the BC value to update.
	  3. The displacements of the Gaussian points are update in pdDataManager::CalcIntgPointDisplRK4(.) function. 

	*/

    // get relative position xi
    double xi1 = x1j - _x1i;
	double xi2 = x2j - _x2i;
	double xi3 = x3j - _x3i;
	// get relative displacement eta
	double eta1, eta2, eta3;
	if (_aiOn)
	{
		if (_mlsOn)
		{			
			double u1j, u2j, u3j; // displacements of the trapezoidal point
			if (_node_j->GetBC())
			{
				// update displacements of trapezoidal points using the boundary condition associated with node_j
				pdBdryCondition* bc = _node_j->GetBC();
				// find boundary condition type, direction flag and time-dependant factor
				pdBdryCondition::BCType bc_type = bc->GetType();
				int dir_x1 = bc->GetDirFlagX1();
				int dir_x2 = bc->GetDirFlagX2();
				int dir_x3 = bc->GetDirFlagX3();
				double fac = bc->GetTimeFactor(_curRunTime);
				// type 1: prescribed displacement boundary conditon
				if (bc_type==pdBdryCondition::DISPL)
				{
					if (dir_x1==1)
					{					
						u1j = bc->GetValue(1)*fac;
					}
					if (dir_x2==1)
					{
						u2j = bc->GetValue(2)*fac;
					}
					if (dir_x3==1)
					{
						u3j = bc->GetValue(3)*fac;
					}
				}
				// type 3: prescribed displacement gradient (12/06/09)
				else if (bc_type == pdBdryCondition::DISPLGRAD)
				{
					// get the global coordinate of this trapezoidal point
					double x1g = x1j + _node_i->GetX1();
					double x2g = x2j + _node_i->GetX2();
					double x3g = x3j + _node_i->GetX3();
					if (dir_x1==1)
					{
						u1j = (bc->GetValue(1, 1)*x1g + bc->GetValue(1, 2)*x2g + bc->GetValue(1, 3)*x3g)*fac;
					}
					if (dir_x2==1)
					{
						u2j = (bc->GetValue(2, 1)*x1g + bc->GetValue(2, 2)*x2g + bc->GetValue(2, 3)*x3g)*fac;
					}
					if (dir_x3==1)
					{
						u3j = (bc->GetValue(3, 1)*x1g + bc->GetValue(3, 2)*x2g + bc->GetValue(3, 3)*x3g)*fac;
					}
				}
				else
				{
					// implement other boundary condition types
				}
			}
			else
			{
				// get the approximated displacement at this trapezoidal (integration) point using MLS approximation
				// the displacements of the contributing nodes of node_j are the trial displacements from RK4 steps
				pdVector* mlsCoeff = _node_j->GetMLSCoeff(_intgPointCounter);
				int ndNumber = _node_j->GetContribNodeListSize();
				u1j = u2j = u3j = 0.0;
				for (int i = 0; i != ndNumber; ++i)
				{
					u1j += _cbNodeDispl[i][0] * mlsCoeff->GetCoeff(i);
					u2j += _cbNodeDispl[i][1] * mlsCoeff->GetCoeff(i);
					u3j += _cbNodeDispl[i][2] * mlsCoeff->GetCoeff(i);
				}
			}
			eta1 = u1j - _u1i;
			eta2 = u2j - _u2i;
			eta3 = u3j - _u3i;
		}
		else
		{
			// the displacements of trapezoidal points are prescribed by boundary conditions
			pdBdryCondition* bc = _node_j->GetBC(); 
			if (bc)
			{
				if (bc->GetType()==pdBdryCondition::DISPLGRAD) // prescribed displacement gradient
				{
					// {eta} = [strain] * {xi}
					eta1 = bc->GetValue(1, 1)*xi1 + bc->GetValue(1, 2)*xi2 + bc->GetValue(1, 3)*xi3;
					eta2 = bc->GetValue(2, 1)*xi1 + bc->GetValue(2, 2)*xi2 + bc->GetValue(2, 3)*xi3;
					eta3 = bc->GetValue(3, 1)*xi1 + bc->GetValue(3, 2)*xi2 + bc->GetValue(3, 3)*xi3;
				}
				else
				{
					// implement other BC type
				}
			}
			else
			{
				// fatal error, no boundary condition is applied on node j (the family node)
				cerr << "pdIntgManager::GetStrainEnergyDensity: no boundary condition is applied on node " << _node_j->GetID() << endl;
				exit(0);
			}
		}
	}
	else if (_fgiOn)
	{
		// get the approximated displacements of the Gaussian points
		// the approximated displacements has been updated in pdDataManager::CalcIntgPointDisplRK4 function

		// get the index for this Gaussian point (it is located in node j's cell)
		//int gIndex = _idx3 + _intgPoint*_idx2 + _intgPoint*_intgPoint*_idx1;
		pdGauss* gPoint = _node_j->GetGaussianPoint(_intgPointCounter);
		double u1j = gPoint->GetU1();
		double u2j = gPoint->GetU2();
		double u3j = gPoint->GetU3();
		eta1 = u1j - _u1i;
		eta2 = u2j - _u2i;
		eta3 = u3j - _u3i;
	}
	else
	{
		// implement other integration scheme
	}
	
	// get |xi| and |xi+eta|
    double r = sqrt(xi1*xi1 + xi2*xi2 + xi3*xi3);
    double p = sqrt(pow(xi1+eta1, 2.0) + pow(xi2+eta2, 2.0) + pow(xi3+eta3, 2.0));
    // find peridynamic force between node i and integration point (per unit volume square)
	pdMaterial* mat_i = _node_i->GetMaterial();
    pdMaterial* mat_j = _node_j->GetMaterial();
    double ecrit = min(mat_i->GetCriticalStretch(), mat_j->GetCriticalStretch());
	double spring = min(mat_i->GetSpringConstant(), mat_j->GetSpringConstant());
	double bdryEffij = 0.5*(_node_i->GetBdryEffFactor() + _node_j->GetBdryEffFactor());
	double beta; //scale factor beta for boundary nodes
	if (_aiOn)
	{
		beta = 1.0;
	}
	else if (_fgiOn)
	{
		// scale factor (beta) based on the location of this Gaussian point
		pdGauss* gPoint = _node_j->GetGaussianPoint(_intgPointCounter);
		double w1 = gPoint->GetWeightX1(); // weight of Gaussian points
		double w2 = gPoint->GetWeightX2();
		double w3 = gPoint->GetWeightX3();
		double w_eff = pow(w1 * w2 * w3, 1.0 / 3.0); // take the geometric average
		double dx = w_eff*_dx / 2.0; // get the effective cell size of this Gaussian point
		beta = this->CalcBetaFactor(r, _delta, dx); // this part can be improved to increase accuracy (07/14/2011)
		if (beta > 0.0)
		{
			++_intgPointEnd[0];
		}
	}
	else
	{
		// implement other integration scheme
	}
	// calculate the peridynamic force
    // Get bond extension and bond stretch
    double u = p - r;
    double stretch = u / r; // current bond stretch
	//_bond->SetStretch(stretch);
    if (u <= 0.0 || (u > 0.0 && stretch <= ecrit)) // bond in compression or not failed and in tension
    {
        pdForce[0] = bdryEffij * beta * spring * stretch * (xi1 + eta1) / p;
        pdForce[1] = bdryEffij * beta * spring * stretch * (xi2 + eta2) / p;
        pdForce[2] = bdryEffij * beta * spring * stretch * (xi3 + eta3) / p;
		// save the elastic energy density (each node shares half of the total energy of one bond)
		double elast = bdryEffij * 0.5 * (0.5 * beta * (spring/r) * u * u);
		if (_fgiOn)
		{
			pdGauss* gPoint = _node_j->GetGaussianPoint(_intgPointCounter);
			double w1 = gPoint->GetWeightX1(); // weight of Gaussian points
			double w2 = gPoint->GetWeightX2();
			double w3 = gPoint->GetWeightX3();
			double w_eff = pow(w1 * w2 * w3, 1.0 / 3.0); // take the geometric average
			elast *= pow(0.5*w_eff, 3);
		}
		else if (_aiOn)
		{
			// the elastic energy needs to be scaled for each trapezoidal point
			// like divide by the total number of trapezoidal points
			// however, the total number of trapezoidal points is not stored and 
			// it can not be seen in this function scope
		}
		else
		{
			// implement other integration methods
		}
		_node_i->AddElastEnergyDensity(elast);
    }
    else // failed and in tension
    {
        pdForce[0] = 0.0;
        pdForce[1] = 0.0;
        pdForce[2] = 0.0;
        // update bond breakage status
        _bond->SetBreak(true);
        _node_i->AddNumDmgBond(1);
    }
}

void pdIntgManager::GetStrainEnergyDensity(double x1j, double x2j, double x3j, double* energy)
{	
    /* 
	
	This function calculate strain energy density at given integration point.
	Currently this function is called once before the dynamic solution starts to calculate the boundary effect compensation factor (07/04/2011).

	*/

    // get xi
    double xi1 = x1j - _x1i;
	double xi2 = x2j - _x2i;
	double xi3 = x3j - _x3i;
	// get eta
	double eta1, eta2, eta3;	
	if (_aiOn)
	{
		if (_mlsOn)
		{
			// update displacements of trapezoidal points using MLS
		}
		else
		{
			// the displacements of trapezoidal points are prescribed by boundary conditions
			pdBdryCondition* bc = _node_j->GetBC(); 
			if (bc)
			{
				if (bc->GetType()==pdBdryCondition::DISPLGRAD) // prescribed displacement gradient
				{
					// {eta} = [strain] * {xi}
					eta1 = bc->GetValue(1, 1)*xi1 + bc->GetValue(1, 2)*xi2 + bc->GetValue(1, 3)*xi3;
					eta2 = bc->GetValue(2, 1)*xi1 + bc->GetValue(2, 2)*xi2 + bc->GetValue(2, 3)*xi3;
					eta3 = bc->GetValue(3, 1)*xi1 + bc->GetValue(3, 2)*xi2 + bc->GetValue(3, 3)*xi3;
				}
				else
				{
					// implement other BC type
				}
			}
			else
			{
				// fatal error, no boundary condition is applied on node j (the family node)
				cerr << "pdIntgManager::GetStrainEnergyDensity: no boundary condition is applied on node " << _node_j->GetID() << endl;
				exit(0);
			}
		}
	}
	else if (_fgiOn)
	{
		if (_mlsOn)
		{
			// get the approximated displacement at integration point
			// get the index for this Gaussian point (it is located in node j's cell)
			//int gIndex = _idx3 + _intgPoint*_idx2 + _intgPoint*_intgPoint*_idx1;
			pdGauss* gPoint = _node_j->GetGaussianPoint(_intgPointCounter);
			double u1j = gPoint->GetU1();
			double u2j = gPoint->GetU2();
			double u3j = gPoint->GetU3();
			eta1 = u1j - _u1i;
			eta2 = u2j - _u2i;
			eta3 = u3j - _u3i;
		}
		else
		{
			// the displacements of trapezoidal points are prescribed by boundary conditions
			pdBdryCondition* bc = _node_j->GetBC(); 
			if (bc)
			{
				if (bc->GetType()==pdBdryCondition::DISPLGRAD) // prescribed displacement gradient
				{
					// {eta} = [strain] * {xi}
					eta1 = bc->GetValue(1, 1)*xi1 + bc->GetValue(1, 2)*xi2 + bc->GetValue(1, 3)*xi3;
					eta2 = bc->GetValue(2, 1)*xi1 + bc->GetValue(2, 2)*xi2 + bc->GetValue(2, 3)*xi3;
					eta3 = bc->GetValue(3, 1)*xi1 + bc->GetValue(3, 2)*xi2 + bc->GetValue(3, 3)*xi3;
				}
				else
				{
					// implement other BC type
				}
			}
			else
			{
				// fatal error, no boundary condition is applied on node j (the family node)
				cerr << "pdIntgManager::GetStrainEnergyDensity: no boundary condition is applied on node " << _node_j->GetID() << endl;
				exit(0);
			}
		}		
	}
	// get |xi| and |xi+eta|
    double r = sqrt(xi1*xi1 + xi2*xi2 + xi3*xi3);
    double p = sqrt(pow(xi1+eta1, 2.0) + pow(xi2+eta2, 2.0) + pow(xi3+eta3, 2.0));
    // find peridynamic force between node i and integration point (per unit volume square)
	pdMaterial* mat_i = _node_i->GetMaterial();
    pdMaterial* mat_j = _node_j->GetMaterial();
    double ecrit = min(mat_i->GetCriticalStretch(), mat_j->GetCriticalStretch());
	double spring = min(mat_i->GetSpringConstant(), mat_j->GetSpringConstant());
	double elastic;
    // Get bond extension and bond stretch
    double u = p - r;
    double stretch = u / r; // current bond stretch
    if (u <= 0.0 || (u > 0.0 && stretch <= ecrit)) // bond in compression or not failed and in tension
    {
        elastic = 0.5*(spring/r)*u*u;		
    }
    else // failed and in tension
    {
        elastic = 0.0;
    }
	if (_fgiOn)
	{
		// scale factor beta based on the location of this Gaussian point
		pdGauss* gPoint = _node_j->GetGaussianPoint(_intgPointCounter);
		double w1 = gPoint->GetWeightX1();
		double w2 = gPoint->GetWeightX2();
		double w3 = gPoint->GetWeightX3();
		double w_eff = pow(w1*w2*w3, 1.0/3.0); // geometric average
		double dx = w_eff*_dx/2.0; // effective cell size of this Gaussian point
		double beta = this->CalcBetaFactor(r, _delta, dx);
		// return value
		energy[0] = beta*elastic;
	}
	else // trapezoidal integration method, beta==1.0
	{
		energy[0] = elastic;
	}
}

// Y limit methods (single argument)
double pdIntgManager::YLimitXLow(double x)
{
    return (_x1j - 0.5 * _dx);
}
double pdIntgManager::YLimitXHigh(double x)
{
    return (_x1j + 0.5 * _dx);
}
double pdIntgManager::YLimitYLow(double x)
{
    return (_x2j - 0.5 * _dx);
}
double pdIntgManager::YLimitYHigh(double x)
{
    return (_x2j + 0.5 * _dx);
}
double pdIntgManager::YLimitZLow(double x)
{
    return (_x3j - 0.5 * _dx);
}
double pdIntgManager::YLimitZHigh(double x)
{
    return (_x3j + 0.5 * _dx);
}        
double pdIntgManager::YLimitHighXYPos(double y)
{
    // return positive solution: x^2 + y^2 = _delta^2 - (x3j+0.5*_dx)^2 based on input y
    double x = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(_x3j + 0.5 * _dx, 2.0)), 2, y);
    if (x < 0.0)
        return 0.0;
    else
        return x;
}
double pdIntgManager::YLimitHighXYNeg(double y)
{
    return -this->YLimitHighXYPos(y);
}
double pdIntgManager::YLimitHighXZPos(double z)
{
    // return positive solution: x^2 + z^2 = _delta^2 - (x2j+0.5*_dx)^2 based on input z
    double x = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(_x2j + 0.5 * _dx, 2.0)), 2, z);
    if (x < 0.0)
        return 0.0;
    else
        return x;
}
double pdIntgManager::YLimitHighXZNeg(double z)
{
    return -this->YLimitHighXZPos(z);
}
double pdIntgManager::YLimitHighYXPos(double x)
{
    // return positive solution: x^2 + y^2 = _delta^2 - (x3j+0.5*_dx)^2 based on input x
    double y = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(_x3j + 0.5 * _dx, 2.0)), 1, x);
    if (y < 0.0)
        return 0.0;
    else
        return y;
}
double pdIntgManager::YLimitHighYXNeg(double x)
{
    return -this->YLimitHighYXPos(x);
}
double pdIntgManager::YLimitHighYZPos(double z)
{
    // return positive solution: z^2 + y^2 = _delta^2 - (x1j+0.5*_dx)^2 based on input z
    double y = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(_x1j + 0.5 * _dx, 2.0)), 1, z);
    if (y < 0.0)
        return 0.0;
    else
        return y;
}
double pdIntgManager::YLimitHighYZNeg(double z)
{
    return -this->YLimitHighYZPos(z);
}
double pdIntgManager::YLimitHighZXPos(double x)
{
    // return positive solution: x^2 + z^2 = _delta^2 - (x2j+0.5*_dx)^2 based on input x
    double z = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(_x2j + 0.5 * _dx, 2.0)), 2, x);
    if (z < 0.0)
        return 0.0;
    else
        return z;
}
double pdIntgManager::YLimitHighZXNeg(double x)
{
    return -this->YLimitHighZXPos(x);
}
double pdIntgManager::YLimitHighZYPos(double y)
{
    // return positive solution: z^2 + y^2 = _delta^2 - (x1j+0.5*_dx)^2 based on input y
    double z = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(_x1j + 0.5 * _dx, 2.0)), 1, y);
    if (z < 0.0)
        return 0.0;
    else
        return z;
}
double pdIntgManager::YLimitHighZYNeg(double y)
{
    return -this->YLimitHighZYPos(y);
}
double pdIntgManager::YLimitLowXYPos(double y)
{
    // return postive x solution: x^2 + y^2 = delta^2 - (x3j - 0.5*dx)^2 based on input y
    double x = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(_x3j - 0.5 * _dx, 2.0)), 2, y);
    if (x < 0.0)
        return 0.0;
    else
        return x;
}
double pdIntgManager::YLimitLowXYNeg(double y)
{
    return -this->YLimitLowXYPos(y);
}
double pdIntgManager::YLimitLowXZPos(double z)
{
    // return postive x solution of : x^2 + z^2 = delta^2 - (x2j - 0.5*dx)^2 based on input z
    double x = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(_x2j - 0.5 * _dx, 2.0)), 2, z);
    if (x < 0.0)
        return 0.0;
    else
        return x;
}
double pdIntgManager::YLimitLowXZNeg(double z)
{
    return -this->YLimitLowXZPos(z);
}
double pdIntgManager::YLimitLowYXPos(double x)
{
    // return postive y solution: x^2 + y^2 = delta^2 - (x3j - 0.5*dx)^2 based on input x
    double y = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(_x3j - 0.5 * _dx, 2.0)), 1, x);
    if (y < 0.0)
        return 0.0;
    else
        return y;
}
double pdIntgManager::YLimitLowYXNeg(double x)
{
    return -this->YLimitLowYXPos(x);
}
double pdIntgManager::YLimitLowYZPos(double z)
{
    // return postive y solution: y^2 + z^2 = delta^2 - (x1j - 0.5*dx)^2 based on input z
    double y = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(_x1j - 0.5 * _dx, 2.0)), 1, z);
    if (y < 0.0)
        return 0.0;
    else
        return y;
}
double pdIntgManager::YLimitLowYZNeg(double z)
{
    return -this->YLimitLowYZPos(z);
}
double pdIntgManager::YLimitLowZXPos(double x)
{
    // return postive x solution of : x^2 + z^2 = delta^2 - (x2j - 0.5*dx)^2 based on input x
    double z = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(_x2j - 0.5 * _dx, 2.0)), 2, x);
    if (z < 0.0)
        return 0.0;
    else
        return z;
}
double pdIntgManager::YLimitLowZXNeg(double x)
{
    return -this->YLimitLowZXPos(x);
}
double pdIntgManager::YLimitLowZYPos(double y)
{
    // return postive y solution: y^2 + z^2 = delta^2 - (x1j - 0.5*dx)^2 based on input y
    double z = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(_x1j - 0.5 * _dx, 2.0)), 1, y);
    if (z < 0.0)
        return 0.0;
    else
        return z;
}
double pdIntgManager::YLimitLowZYNeg(double y)
{
    return -this->YLimitLowZYPos(y);
}

// Z limit methods (double argument)
double pdIntgManager::ZLimitZSpherePos(double x, double y)
{
    // return positive  z solution of: x^2 + y^2 + z^2 = delta^2
    double z_sq = _delta * _delta - pow(x, 2.0) - pow(y, 2.0);
    if (z_sq < 0.0)
        return 0.0;
    else
        return sqrt(z_sq);
}
double pdIntgManager::ZLimitZSphereNeg(double x, double y)
{
    return -this->ZLimitZSpherePos(x, y);
}
double pdIntgManager::ZLimitZLow(double x, double y)
{
    return (_x3j - 0.5 * _dx);
}
double pdIntgManager::ZLimitZHigh(double x, double y)
{
    return (_x3j + 0.5 * _dx);
}
double pdIntgManager::ZLimitYSpherePos(double x, double y)
{
    // return positive y solution of: x^2 + y^2 + z^2 = delta^2
    double y_sq = _delta * _delta - pow(x, 2.0) - pow(y, 2.0);
    if (y_sq < 0.0)
        return 0.0;
    else
        return sqrt(y_sq);
}
double pdIntgManager::ZLimitYSphereNeg(double x, double y)
{
    return -this->ZLimitYSpherePos(x, y);
}
double pdIntgManager::ZLimitYLow(double x, double y)
{
    return (_x2j - 0.5 * _dx);
}
double pdIntgManager::ZLimitYHigh(double x, double y)
{
    return (_x2j + 0.5 * _dx);
}
double pdIntgManager::ZLimitXSpherePos(double x, double y)
{
    // return positive x solution of: x^2 + y^2 + z^2 = delta^2
    double x_sq = _delta * _delta - pow(x, 2.0) - pow(y, 2.0);
    if (x_sq < 0.0)
        return 0.0;
    else
        return sqrt(x_sq);
}
double pdIntgManager::ZLimitXSphereNeg(double x, double y)
{
    return -this->ZLimitXSpherePos(x, y);
}
double pdIntgManager::ZLimitXLow(double x, double y)
{
    return (_x1j - 0.5 * _dx);
}
double pdIntgManager::ZLimitXHigh(double x, double y)
{
    return (_x1j + 0.5 * _dx);
}

FUNC_T pdIntgManager::SelectIntegrand(string funcName) // select integrand function
{
	if (funcName=="UnitValue")
		return &pdIntgManager::UnitValue;
	else if (funcName=="CalcBondBasedForce")
		return &pdIntgManager::CalcBondBasedForce;
	else if (funcName=="GetStrainEnergyDensity")
		return &pdIntgManager::GetStrainEnergyDensity;
	else
		return NULL;
}

FUNC_S pdIntgManager::SelectYLimitFunc(string funcName) // select integration limits for Y
{
	if (funcName=="YLimitXLow")
		return &pdIntgManager::YLimitXLow;
	else if (funcName=="YLimitXHigh")
		return &pdIntgManager::YLimitXHigh;
	else if (funcName=="YLimitYLow")
		return &pdIntgManager::YLimitYLow;
	else if (funcName=="YLimitYHigh")
		return &pdIntgManager::YLimitYHigh;
	else if (funcName=="YLimitZLow")
		return &pdIntgManager::YLimitZLow;
	else if (funcName=="YLimitZHigh")
		return &pdIntgManager::YLimitZHigh;
	else if (funcName=="YLimitHighXYPos")
		return &pdIntgManager::YLimitHighXYPos;
	else if (funcName=="YLimitHighXYNeg")
		return &pdIntgManager::YLimitHighXYNeg;
	else if (funcName=="YLimitHighXZPos")
		return &pdIntgManager::YLimitHighXZPos;
	else if (funcName=="YLimitHighXZNeg")
		return &pdIntgManager::YLimitHighXZNeg;
	else if (funcName=="YLimitHighYXPos")
		return &pdIntgManager::YLimitHighYXPos;
	else if (funcName=="YLimitHighYXNeg")
		return &pdIntgManager::YLimitHighYXNeg;
	else if (funcName=="YLimitHighYZPos")
		return &pdIntgManager::YLimitHighYZPos;
	else if (funcName=="YLimitHighYZNeg")
		return &pdIntgManager::YLimitHighYZNeg;
	else if (funcName=="YLimitHighZXPos")
		return &pdIntgManager::YLimitHighZXPos;
	else if (funcName=="YLimitHighZXNeg")
		return &pdIntgManager::YLimitHighZXNeg;
	else if (funcName=="YLimitHighZYPos")
		return &pdIntgManager::YLimitHighZYPos;
	else if (funcName=="YLimitHighZYNeg")
		return &pdIntgManager::YLimitHighZYNeg;
	else if (funcName=="YLimitLowXYPos")
		return &pdIntgManager::YLimitLowXYPos;
	else if (funcName=="YLimitLowXYNeg")
		return &pdIntgManager::YLimitLowXYNeg;
	else if (funcName=="YLimitLowXZPos")
		return &pdIntgManager::YLimitLowXZPos;
	else if (funcName=="YLimitLowXZNeg")
		return &pdIntgManager::YLimitLowXZNeg;
	else if (funcName=="YLimitLowYXPos")
		return &pdIntgManager::YLimitLowYXPos;
	else if (funcName=="YLimitLowYXNeg")
		return &pdIntgManager::YLimitLowYXNeg;
	else if (funcName=="YLimitLowYZPos")
		return &pdIntgManager::YLimitLowYZPos;
	else if (funcName=="YLimitLowYZNeg")
		return &pdIntgManager::YLimitLowYZNeg;
	else if (funcName=="YLimitLowZXPos")
		return &pdIntgManager::YLimitLowZXPos;
	else if (funcName=="YLimitLowZXNeg")
		return &pdIntgManager::YLimitLowZXNeg;
	else if (funcName=="YLimitLowZYPos")
		return &pdIntgManager::YLimitLowZYPos;
	else if (funcName=="YLimitLowZYNeg")
		return &pdIntgManager::YLimitLowZYNeg;
	else
		return NULL;
}

FUNC_D pdIntgManager::SelectZLimitFunc(string funcName) // select integration limits for Z
{
	if (funcName=="ZLimitXLow")
		return &pdIntgManager::ZLimitXLow;
	else if (funcName=="ZLimitXHigh")
		return &pdIntgManager::ZLimitXHigh;
	else if (funcName=="ZLimitYLow")
		return &pdIntgManager::ZLimitYLow;
	else if (funcName=="ZLimitYHigh")
		return &pdIntgManager::ZLimitYHigh;
	else if (funcName=="ZLimitZLow")
		return &pdIntgManager::ZLimitZLow;
	else if (funcName=="ZLimitZHigh")
		return &pdIntgManager::ZLimitZHigh;
	else if (funcName=="ZLimitXSpherePos")
		return &pdIntgManager::ZLimitXSpherePos;
	else if (funcName=="ZLimitXSphereNeg")
		return &pdIntgManager::ZLimitXSphereNeg;
	else if (funcName=="ZLimitYSpherePos")
		return &pdIntgManager::ZLimitYSpherePos;
	else if (funcName=="ZLimitYSphereNeg")
		return &pdIntgManager::ZLimitYSphereNeg;
	else if (funcName=="ZLimitZSpherePos")
		return &pdIntgManager::ZLimitZSpherePos;
	else if (funcName=="ZLimitZSphereNeg")
		return &pdIntgManager::ZLimitZSphereNeg;
	else
		return NULL;
}

FUNC_S pdIntgManager::SelectSingleArguFunc(string funcName)
{
	if (funcName=="YLimitXLow")
		return &pdIntgManager::YLimitXLow;
	else if (funcName=="YLimitXHigh")
		return &pdIntgManager::YLimitXHigh;
	else if (funcName=="YLimitYLow")
		return &pdIntgManager::YLimitYLow;
	else if (funcName=="YLimitYHigh")
		return &pdIntgManager::YLimitYHigh;
	else if (funcName=="YLimitZLow")
		return &pdIntgManager::YLimitZLow;
	else if (funcName=="YLimitZHigh")
		return &pdIntgManager::YLimitZHigh;
	else
		return NULL;
}

double pdIntgManager::SolveBinaryQuadric(double c_x, double c_y, double r_sq, int p, double xy)
{
	// this method solves equation (X-c_x)^2 + (Y-c_y)^2 = r_sq with input xy :
    // if p=1, X=xy is know, return abs(Y-c_y) if solveable, return -1 if unsolveable
    // if p=2, Y=xy is known, return abs(X-c_x) if solveable, return -1 if unsolveable
	
	if (p==1)
	{
		double y_sq = r_sq - pow(xy-c_x, 2.0);
		if (y_sq<0.0)
			return -1.0;
		else
			return sqrt(y_sq);
	}
	else
	{
		double x_sq = r_sq - pow(xy-c_y, 2.0);
		if (x_sq<0.0)
			return -1.0;
		else
			return sqrt(x_sq);
	}
}

int pdIntgManager::FindIntersectionPointOnEightEdges(vector<string> funcName, double* iPoints)
{
    // initial binary operation
    int b = 0; // store the intersecting data, 1 for intersected, 0 for not intersected
    int bm = 1;

    // The followings identify intersection points of horizon sphere and the edges on two surfaces of the family node cell
    // function1 for intersection of the horizon sphere with the high surface of node j's cell: 
    //     (x-x1i)^2 + (y-x2i)^2 + (z-x3i)^2 = _delta^2, z is determined by the surface function of funcName[5] (zhigh)
    // function2 for intersection of the horizon sphere with the low surface of node j's cell: 
    //     (x-x1i)^2 + (y-x2i)^2 + (z-x3i)^2 = _delta^2, z is determined by the surface function of funcName[4] (zlow)
    // because we are using local coordinate, x1i=x2i=x3i=0
    // value of xlow, xhigh, ylow, yhigh determine the location of the edges

    double xlow = (this->*SelectSingleArguFunc(funcName[0]))(0);
    double xhigh = (this->*SelectSingleArguFunc(funcName[1]))(0);
    double ylow = (this->*SelectSingleArguFunc(funcName[2]))(0);
    double yhigh = (this->*SelectSingleArguFunc(funcName[3]))(0);
    double zlow = (this->*SelectSingleArguFunc(funcName[4]))(0); // low arc
    double zhigh = (this->*SelectSingleArguFunc(funcName[5]))(0); // high arc
    // relative location of node j (norm)
    double n1, n2, n3;
    if (_x1j == 0.0)
        n1 = 1.0;
    else
        n1 = _x1j / abs(_x1j);
    if (_x2j == 0.0)
        n2 = 1.0;
    else
        n2 = _x2j / abs(_x2j);
    if (_x3j == 0.0)
        n3 = 1.0;
    else
        n3 = _x3j / abs(_x3j);
    // T1 (top left edge): z=zhigh, x=xlow, find the Y coordinate of intersection point
    double temp = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(zhigh, 2.0)), 1, xlow);
    if (temp >= 0.0)
    {
        temp *= n2; //adjust solution based the location of node j in the family cooridnate
        if (temp >= ylow && temp <= yhigh) // function1 intersects with T1
        {
            b = b | bm;
            *iPoints = temp;
        }
        else  // intersection point is not within the edge
        {
        }
    }
    else // no intersection point
    {
    }
    // T3 (top right edge): z=zhigh, x=xhigh, find the Y coordinate of intersection point
    bm <<= 1;
    temp = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(zhigh, 2.0)), 1, xhigh);
    if (temp >= 0.0)
    {
        temp *= n2;
        if (temp >= ylow && temp <= yhigh)
        {
            b = b | bm;
            *(iPoints+1) = temp;
        }
    }
    // T2 (top front edge): z=zhigh, y=ylow, find the X coordinate of intersection point
    bm <<= 1;
    temp = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(zhigh, 2.0)), 2, ylow);
    if (temp >= 0.0) 
    {
        temp *= n1;
        if (temp >= xlow && temp <= xhigh)
        {
            b = b | bm;
            *(iPoints+2) = temp;
        }
    }
    // T4 (top back edge): z=zhigh, y=yhigh, find the X coordinate of intersection point
    bm <<= 1;
    temp = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(zhigh, 2.0)), 2, yhigh);
    if (temp >= 0.0)
    {
        temp *= n1;
        if (temp >= xlow && temp <= xhigh)
        {
            b = b | bm;
            *(iPoints+3) = temp;
        }
    }
    // B1 (bottom left edge): z=zlow, x=xlow, find the Y coordinate of intersection point
    bm <<= 1;
    temp = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(zlow, 2.0)), 1, xlow);
    if (temp >= 0.0)
    {
        temp *= n2;
        if (temp >= ylow && temp <= yhigh)
        {
            b = b | bm;
            *(iPoints+4) = temp;
        }
    }
    // B2 (bottom right edge): z=zlow, x=xhigh, find the Y coordinate of intersection point
    bm <<= 1;
    temp = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(zlow, 2.0)), 1, xhigh);
    if (temp >= 0.0)
    {
        temp *= n2;
        if (temp >= ylow && temp <= yhigh)
        {
            b = b | bm;
            *(iPoints+5) = temp;
        }
    }
    // B3 (bottom front edge): z=zlow, y=ylow, find the X coordinate of intersection point
    bm <<= 1;
    temp = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(zlow, 2.0)), 2, ylow);
    if (temp >= 0.0)
    {
        temp *= n1;
        if (temp >= xlow && temp <= xhigh)
        {
            b = b | bm;
            *(iPoints+6) = temp;
        }
    }
    // B4 (bottom back edge): z=zlow, y=yhigh, find the X coordinate of intersection point
    bm <<= 1;
    temp = SolveBinaryQuadric(0.0, 0.0, (_delta * _delta - pow(zlow, 2.0)), 2, yhigh);
    if (temp >= 0.0)
    {
        temp *= n1;
        if (temp>= xlow && temp<= xhigh)
        {
            b = b | bm;
            *(iPoints+7) = temp;
        }
    }
    // return the binary results
    return b;
}

int pdIntgManager::FindIntersectionPointOnFourEdges(vector<string> funcName, double* iPoints)
{
    /* 
      This function is used to find intersection points on four other edges if horizon sphere does not 
      intersects with one of the surfaces when method "FindIntersectionPointOnEightEdges" is called
    */
    // initial binary operation
    int b = 0; // store the intersecting data, 1 for intersected, 0 for not intersected
    int bm = 1;

    // The followings identify intersection points of horizon sphere and the four edges defined in funcName
    // x^2 + y^2 + z^2 = _delta^2, x and y are given by the funcName
    // value of xlow, xhigh, ylow, yhigh determine the location of the edges

    double xlow = (this->*SelectSingleArguFunc(funcName[0]))(0);
    double xhigh = (this->*SelectSingleArguFunc(funcName[1]))(0);
    double ylow = (this->*SelectSingleArguFunc(funcName[2]))(0);
    double yhigh = (this->*SelectSingleArguFunc(funcName[3]))(0);
    double zlow = (this->*SelectSingleArguFunc(funcName[4]))(0);
    double zhigh = (this->*SelectSingleArguFunc(funcName[5]))(0);
    // relative location of node j (norm)
    double n1, n2, n3;
    if (_x1j == 0.0)
        n1 = 1.0;
    else
        n1 = _x1j / abs(_x1j);
    if (_x2j == 0.0)
        n2 = 1.0;
    else
        n2 = _x2j / abs(_x2j);
    if (_x3j == 0.0)
        n3 = 1.0;
    else
        n3 = _x3j / abs(_x3j);
    // edge V1: x=xlow, y=ylow, find the Z coordinate of intersection point
    double temp = _delta * _delta - xlow * xlow - ylow * ylow;
    if (temp >= 0.0)
    {
        double z_sq = sqrt(temp) * n3; //adjust solution based the location of node j in the family cooridnate
        if (z_sq >= zlow && z_sq <= zhigh)
        {
            b = b | bm;
            iPoints[0] = z_sq;
        }
    }
    bm <<= 1;
    // edge V2: x=xhigh, y=ylow, find the Z coordinate of intersection point
    temp = _delta * _delta - xhigh * xhigh - ylow * ylow;
    if (temp >= 0.0)
    {
        double z_sq = sqrt(temp) * n3; //adjust solution based the location of node j in the family cooridnate
        if (z_sq >= zlow && z_sq <= zhigh)
        {
            b = b | bm;
            iPoints[1] = z_sq;
        }
    }
    bm <<= 1;
    // edge V3: x=xhigh, y=yhigh, find the Z coordinate of intersection point
    temp = _delta * _delta - xhigh * xhigh - yhigh * yhigh;
    if (temp >= 0.0)
    {
        double z_sq = sqrt(temp) * n3; //adjust solution based the location of node j in the family cooridnate
        if (z_sq >= zlow && z_sq <= zhigh)
        {
            b = b | bm;
            iPoints[2] = z_sq;
        }
    }
    bm <<= 1;
    // edge V4: x=xlow, y=yhigh, find the Z coordinate of intersection point
    temp = _delta * _delta - xlow * xlow - yhigh * yhigh;
    if (temp >= 0.0)
    {
        double z_sq = sqrt(temp) * n3; //adjust solution based the location of node j in the family cooridnate
        if (z_sq >= zlow && z_sq <= zhigh)
        {
            b = b | bm;
            iPoints[3] = z_sq;
        }
    }            
    // return the binary results
    return b;
}

void pdIntgManager::CalcConfigurationType1Sub1(string dir1, string dir2, string dir3, string itg, double dist, double* pdForce)
{
    // string dir1, dir2, dir3 are the integration order (outer to inner)
    double low = -sqrt(_delta * _delta - pow((dist - 0.5 * _dx), 2.0));
    double high = sqrt(_delta * _delta - pow((dist - 0.5 * _dx), 2.0));
    vector<string> flag; // table index for limit function
    double innerCoord = this->GetCoordAtDirection(dir3);
    if (innerCoord > 0.0)
    {
        flag.push_back("YLimitLow" + dir2 + dir1 + "Neg");
        flag.push_back("YLimitLow" + dir2 + dir1 + "Pos");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "SpherePos");
        flag.push_back(dir1 + dir2 + dir3);
    }
    else
    {
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Neg");
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Pos");
        flag.push_back("ZLimit" + dir3 + "SphereNeg");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
    }
    CalcIntg3D(itg, low, high, flag, pdForce);
} 

void pdIntgManager::CalcConfigurationType1Sub2(string dir1, string dir2, string dir3, string itg, double dist, double* pdForce)
{
	double low = -sqrt(_delta * _delta - pow((dist - 0.5 * _dx), 2.0) - 0.25 * _dx * _dx);
    double high = sqrt(_delta * _delta - pow((dist - 0.5 * _dx), 2.0) - 0.25 * _dx * _dx);
    double outerCoord = this->GetCoordAtDirection(dir1);
    double innerCoord = this->GetCoordAtDirection(dir3);
    vector<string> flag; // table index for limit function
	double pdForceStep[3];
	pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
    if (innerCoord > 0.0)
    {
        flag.push_back("YLimitLow" + dir2 + dir1 + "Neg");
        flag.push_back("YLimitLow" + dir2 + dir1 + "Pos");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "SpherePos");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, high, outerCoord + 0.5 * _dx, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
		flag.clear();
        flag.push_back("YLimitLow" + dir2 + dir1 + "Neg");
        flag.push_back("YLimitLow" + dir2 + dir1 + "Pos");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "SpherePos");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, outerCoord - 0.5 * _dx, low, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "SpherePos");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, low, high, flag, pdForceStep);
        ArrayAddition(pdForceStep, pdForce, 3);
    }
    else
    {
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Neg");
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Pos");
        flag.push_back("ZLimit" + dir3 + "SphereNeg");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, high, outerCoord + 0.5 * _dx, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Neg");
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Pos");
        flag.push_back("ZLimit" + dir3 + "SphereNeg");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, outerCoord - 0.5 * _dx, low, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "SphereNeg");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, low, high, flag, pdForceStep);
        ArrayAddition(pdForceStep, pdForce, 3);
    }
}

void pdIntgManager::CalcConfigurationType1Sub3(string dir1, string dir2, string dir3, string itg, double* pdForce)
{
    double outerCoord = this->GetCoordAtDirection(dir1);
    double innerCoord = this->GetCoordAtDirection(dir3);
    vector<string> flag; // table index for limit function
    if (innerCoord > 0.0)
    {
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "SpherePos");
        flag.push_back(dir1 + dir2 + dir3);
    }
    else
    {
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "SphereNeg");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
    }
    CalcIntg3D(itg, outerCoord - 0.5 * _dx, outerCoord + 0.5 * _dx, flag, pdForce);
}

void pdIntgManager::CalcConfigurationType1Sub4(string dir1, string dir2, string dir3, string itg, double dist, double* pdForce)
{
    double low = -sqrt(_delta * _delta - pow(dist + 0.5 * _dx, 2.0));
    double high = sqrt(_delta * _delta - pow(dist + 0.5 * _dx, 2.0));
    double outerCoord = this->GetCoordAtDirection(dir1);
    double innerCoord = this->GetCoordAtDirection(dir3);
    vector<string> flag; // table index for limit function
	double pdForceStep[3];
	pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
    if (innerCoord > 0.0)
    {
        // center trunk
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Pos");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "SpherePos");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, low, high, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Neg");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "SpherePos");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, low, high, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "SpherePos");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, high, outerCoord + 0.5 * _dx, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "SpherePos");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, outerCoord - 0.5 * _dx, low, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        // tip
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Neg");
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Pos");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "SpherePos");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, low, high, flag, pdForceStep);
        ArrayAddition(pdForceStep, pdForce, 3);
    }
    else
    {
        flag.push_back("YLimitLow" + dir2 + dir1 + "Pos");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "SphereNeg");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, low, high, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimitLow" + dir2 + dir1 + "Neg");
        flag.push_back("ZLimit" + dir3 + "SphereNeg");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, low, high, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "SphereNeg");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, high, outerCoord + 0.5 * _dx, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "SphereNeg");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, outerCoord - 0.5 * _dx, low, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimitLow" + dir2 + dir1 + "Neg");
        flag.push_back("YLimitLow" + dir2 + dir1 + "Pos");
        flag.push_back("ZLimit" + dir3 + "SphereNeg");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, low, high, flag, pdForceStep);
        ArrayAddition(pdForceStep, pdForce, 3);
    }
}

void pdIntgManager::CalcConfigurationType1Sub5(string dir1, string dir2, string dir3, string itg, double dist, double* pdForce)
{
    double low = -sqrt(_delta * _delta - pow(dist + 0.5 * _dx, 2.0) - 0.25 * _dx * _dx);
    double high = sqrt(_delta * _delta - pow(dist + 0.5 * _dx, 2.0) - 0.25 * _dx * _dx);
    double outerCoord = this->GetCoordAtDirection(dir1);
    double innerCoord = this->GetCoordAtDirection(dir3);
    vector<string> flag; // table index for limit function
	double pdForceStep[3];
	pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
    if (innerCoord > 0.0)
    {
        // four corner
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Pos");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "SpherePos");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, outerCoord - 0.5 * _dx, low, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Pos");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "SpherePos");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, high, outerCoord + 0.5 * _dx, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Neg");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "SpherePos");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, outerCoord - 0.5 * _dx, low, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Neg");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "SpherePos");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, high, outerCoord + 0.5 * _dx, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        // center trunk
        flag.clear();
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Neg");
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Pos");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, outerCoord - 0.5 * _dx, low, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Neg");
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Pos");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, high, outerCoord + 0.5 * _dx, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, low, high, flag, pdForceStep);
        ArrayAddition(pdForceStep, pdForce, 3);
    }
    else
    {
        // four corner
        flag.push_back("YLimitLow" + dir2 + dir1 + "Pos");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "SphereNeg");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, outerCoord - 0.5 * _dx, low, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimitLow" + dir2 + dir1 + "Pos");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "SphereNeg");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, high, outerCoord + 0.5 * _dx, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimitLow" + dir2 + dir1 + "Neg");
        flag.push_back("ZLimit" + dir3 + "SphereNeg");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, outerCoord - 0.5 * _dx, low, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimitLow" + dir2 + dir1 + "Neg");
        flag.push_back("ZLimit" + dir3 + "SphereNeg");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, high, outerCoord + 0.5 * _dx, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        // center trunk
        flag.clear();
        flag.push_back("YLimitLow" + dir2 + dir1 + "Neg");
        flag.push_back("YLimitLow" + dir2 + dir1 + "Pos");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, outerCoord - 0.5 * _dx, low, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimitLow" + dir2 + dir1 + "Neg");
        flag.push_back("YLimitLow" + dir2 + dir1 + "Pos");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, high, outerCoord + 0.5 * _dx, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, low, high, flag, pdForceStep);
        ArrayAddition(pdForceStep, pdForce, 3);
    }
}

double pdIntgManager::GetCoordAtDirection(string dir)
{
    if (dir=="X")
        return _x1j;
	else if (dir=="Y")
        return _x2j;
	else if (dir=="Z")
        return _x3j;
	else
		return 0.0; // should never get here
}

void pdIntgManager::CalcConfigurationType2Sub1(string dir1, string dir2, string dir3, string itg, double* pdForce)
{
    double outerCoord = this->GetCoordAtDirection(dir1);
    double midCoord = this->GetCoordAtDirection(dir2);
    double innerCoord = this->GetCoordAtDirection(dir3);
    vector<string> flag; // table index for limit function
    if (midCoord > 0.0 && innerCoord < 0.0)
    {                
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Pos");
        flag.push_back("ZLimit" + dir3 + "SphereNeg");
        flag.push_back("ZLimit" + dir3 + "High");
    }
    else if (midCoord > 0.0 && innerCoord > 0.0)
    {              
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimitLow" + dir2 + dir1 + "Pos");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "SpherePos");
    }
    else if (midCoord < 0.0 && innerCoord > 0.0)
    {              
        flag.push_back("YLimitLow" + dir2 + dir1 + "Neg");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "SpherePos");
    }
    else // (midCoord < 0.0 && innerCoord < 0.0)
    {              
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Neg");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "SphereNeg");
        flag.push_back("ZLimit" + dir3 + "High");
    }
    flag.push_back(dir1 + dir2 + dir3);
    CalcIntg3D(itg, outerCoord - 0.5 * _dx, outerCoord + 0.5 * _dx, flag, pdForce);
}

void pdIntgManager::CalcConfigurationType2Sub2(string dir1, string dir2, string dir3, string itg, double* pdForce)
{
    double outerCoord = this->GetCoordAtDirection(dir1);
    double midCoord = this->GetCoordAtDirection(dir2);
    double innerCoord = this->GetCoordAtDirection(dir3);
    vector<string> flag; // table index for limit function
	double pdForceStep[3];
	pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
    if (midCoord > 0.0 && innerCoord < 0.0)
    {
        flag.push_back("YLimitLow" + dir2 + dir1 + "Pos");
        flag.push_back("YLimit" + dir2 + "High");      
        flag.push_back("ZLimit" + dir3 + "SphereNeg");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, outerCoord - 0.5 * _dx, outerCoord + 0.5 * _dx, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();     
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimitLow" + dir2 + dir1 + "Pos");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, outerCoord - 0.5 * _dx, outerCoord + 0.5 * _dx, flag, pdForceStep);
        ArrayAddition(pdForceStep, pdForce, 3);
    }
    else if (midCoord > 0.0 && innerCoord > 0.0)
    {
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Pos");
        flag.push_back("YLimit" + dir2 + "High");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "SpherePos");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, outerCoord - 0.5 * _dx, outerCoord + 0.5 * _dx, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();      
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Pos");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, outerCoord - 0.5 * _dx, outerCoord + 0.5 * _dx, flag, pdForceStep);
        ArrayAddition(pdForceStep, pdForce, 3);
    }
    else if (midCoord < 0.0 && innerCoord > 0.0)
    {         
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Neg");
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "SpherePos");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, outerCoord - 0.5 * _dx, outerCoord + 0.5 * _dx, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimitHigh" + dir2 + dir1 + "Neg");
        flag.push_back("YLimit" + dir2 + "High");  
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, outerCoord - 0.5 * _dx, outerCoord + 0.5 * _dx, flag, pdForceStep);
        ArrayAddition(pdForceStep, pdForce, 3);
    }
    else // (midCoord < 0.0 && innerCoord < 0.0)
    {               
        flag.push_back("YLimit" + dir2 + "Low");
        flag.push_back("YLimitLow" + dir2 + dir1 + "Neg");
        flag.push_back("ZLimit" + dir3 + "SphereNeg");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, outerCoord - 0.5 * _dx, outerCoord + 0.5 * _dx, flag, pdForceStep);
		ArrayAddition(pdForceStep, pdForce, 3);
        flag.clear();
        flag.push_back("YLimitLow" + dir2 + dir1 + "Neg");
        flag.push_back("YLimit" + dir2 + "High");         
        flag.push_back("ZLimit" + dir3 + "Low");
        flag.push_back("ZLimit" + dir3 + "High");
        flag.push_back(dir1 + dir2 + dir3);
        CalcIntg3D(itg, outerCoord - 0.5 * _dx, outerCoord + 0.5 * _dx, flag, pdForceStep);
        ArrayAddition(pdForceStep, pdForce, 3);
    }           
}

void pdIntgManager::GetAppDispl(double x, double y, double z, vector<pdNode*> cbList, double* appDispl)
{
    /* This methods returns approximated displacement at integration point using Moving Least Square method.
     * Input:  x, y, z are the coordinates of integration point. Because the coordinate of node j is now 
	 *         in the family coordinate, but the coordinate from the nodes in the contributing list is in the
	 *         global coordinate, the location of integration points need to be shifted back to global coordinate
     *         rw is the support domain size
     *         cbList is the contributing node list
     * Output: appDispl (array) is the approximated displacement at integration point
     */

	// shift the integration point back to global coordinate
	x += _node_i->GetX1();
	y += _node_i->GetX2();
	z += _node_i->GetX3();
    int dim = POLYBASISDIM; // basis vector dimension
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
    //double normMomentMatrix = lhs.RowSumNorm(); // infinite norm of the moment matrix
    double normMomentMatrix = lhs->FBnorm(); // Frobenius norm
    // set up ds (rhs)
    pdVector* rhs = this->GetBasis(x, y, z, dim);
    pdVector* sol = new pdVector(dim);
    // solve system of equations, store the results in sol
    lhs->Gauss(rhs, sol);
	// check if the moment matrix's condition number
    //double normMomentMatrixInv = lhs.Norm(); // row-sum norm
    double normMomentMatrixInv = lhs->FBnorm();
    if (normMomentMatrix * normMomentMatrixInv > ILLCONDITION)
    {
        cout << "Moment matrix is ill-conditioned";
		exit(0);
    }
	else
	{
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
		//bool ifConsistent = this->ConsistencyCheck(x, y, z, cbList, coeff);
		// get the approximated displacement at this integration point
		appDispl[0] = 0.0;
		appDispl[1] = 0.0;
		appDispl[2] = 0.0;
		for (int i = 0; i != ndNumber; ++i)
		{
			appDispl[0] += cbList[i]->GetU1() * coeff->GetCoeff(i);
			appDispl[1] += cbList[i]->GetU2() * coeff->GetCoeff(i);
			appDispl[2] += cbList[i]->GetU3() * coeff->GetCoeff(i);
		}
		// delete temporary objects	
		delete lhs;
		delete rhs;
		delete sol;
		delete b;
		delete b2;
		delete coeff;
	}
}

pdVector* pdIntgManager::GetBasis(double x, double y, double z, int dim)
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

double pdIntgManager::GetWeight(pdNode* n, double x, double y, double z, double rw)
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

bool pdIntgManager::ConsistencyCheck(double x, double y, double z, vector<pdNode*> cbList, pdVector* coeff)
{
	// unity consistency check
	double uniSum = 0.0;
	for (int i = 0; i != coeff->GetNumRows(); ++i)
	    uniSum += coeff->GetCoeff(i);
	// approx. order consistency check
	double appOrder[POLYBASISDIM-1];
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

int pdIntgManager::GetEndIntgPoint(int idx)
{
	// return the endding integration point number in the direction given by the input index idx
	return _intgPointEnd[idx];
}

double pdIntgManager::CalcBetaFactor(double dist, double delta, double dx)
{
	/* 
	
	This function calculated the boundary correct factor beta based on the distance between integration point and souce node 
	It is originated from the EMU code.

	Input:
	  1. dist: the distance between the source node and integration point.
	  2. delta: the horizon size.
	  3. dx: the effective subcell size.

    Output:
	  1. beta.

	*/

	if (dist <= delta - 0.5*dx)
	{
		return 1.0;
	}
	else if (dist > delta - 0.5*dx && dist <= delta + 0.5*dx)
	{
		return (delta + 0.5*dx - dist) / dx;
	}
	else
	{
		return 0.0;	
	}
}

void pdIntgManager::sobseq(int *n, double x[])
{
	// quasi-random point generator (from Numerical Recipes in C: 2nd edition P312-313)
	const int MAXBIT = 30;
	const int MAXDIM = 6;
	// when n is negative, internally initializes a set of MAXBIT direction numbers for each of MAXDIM different Sobol's sequences.
	// when n is positive (but <=MAXDIM), returns as the vector x[1,...n]
	// the next value from n of these sequences. (n must not be changed between initializations)
	int j, k, l;
	unsigned long i, im, ipp;
	static double fac;
	static unsigned long in, ix[MAXDIM+1], *iu[MAXBIT+1];
	static unsigned long mdeg[MAXDIM+1] = {0,1,2,3,3,4,4};
	static unsigned long ip[MAXDIM+1] = {0,0,1,1,2,1,4};
	static unsigned long iv[MAXDIM*MAXBIT+1] = {0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};

	if (*n<0) // initialize, don't return a vector
	{
		for (k=1;k<=MAXDIM;k++)
			ix[k] = 0;
		in = 0;
		if (iv[1]!=1)
			return;
		fac = 1.0/(1L<< MAXBIT);
		for (j=1, k=0;j<=MAXBIT;j++, k+=MAXDIM)
			iu[j] = &iv[k];
		// to allow both 1D and 2D addressing
		for (k=1;k<=MAXDIM;k++)
		{
			for (j=1;j<=mdeg[k];j++)
				iu[j][k] <<= (MAXBIT-j);
			// stored values only require normalization
			for (j=mdeg[k]+1;j<=MAXBIT;j++) // use the recurrence to get other values
			{
				ipp = ip[k];
				i = iu[j-mdeg[k]][k];
				i ^= (i >> mdeg[k]);
				for (l=mdeg[k]-1;l>=1;l--)
				{
					if (ipp & 1)
						i ^= iu[j-1][k];
					ipp >>= 1;
				}
				iu[j][k] = i;
			}
		}
	}
	else // calculate the next vector in the sequence
	{
		im = in++;
		for (j=1;j<=MAXBIT;j++) // find the rightmost zero bit
		{
			if ( !(im & 1) )
				break;
			im >>= 1;
		}
		if (j>MAXBIT)
			cout << "MAXBIT too small in sobseq" << endl;
		im = (j-1)*MAXDIM;
		for (k=1;k<=min(*n, MAXDIM);k++) // XOR the appropriate direction number into each component of the vector and convert to a floating number
		{
			ix[k] ^= iv[im+k];
			x[k] = ix[k]*fac;
		}
	}
}

double pdIntgManager::qgaus(FUNC_S func, double a, double b)
{
	// returns the integral of the function func between a and b, by ten-point Gauss-Legendra integration
	// the function is evaluated exactly ten times at interior points in the range of integration

	// the _abscissas and _weights, first value of each array not used
	static double x[] = {0.0, 0.1488743389, 0.4333953941, 0.6794095682, 0.8650633666, 0.9739065285};
	static double w[] = {0.0, 0.2955242247, 0.2692667193, 0.2190863625, 0.1494513491, 0.0666713443};

	double xm = 0.5*(b+a);
	double xr = 0.5*(b-a);
	double s = 0.;
	// will be twice the average value of the function, since the ten _weights (five numbers above each used twice) sum to 2
	for (int j=1;j<=5;j++)
	{
		double dx = xr*x[j];
		s += w[j]*((this->*func)(xm+dx) + (this->*func)(xm-dx));
	}
	return s *= xr; // scale the answer to the range of integration
}

void pdIntgManager::ArrayAddition(double *source, double *target, const int length)
{
	// This function adds the elements in array 'target' to array 'source'.
	// element index starts from 0 to 'length'-1
	// user has to make sure the size of array 'source' and 'target' are both equal to 'length'.
	try
	{
		/*int tarSize = sizeof(target)/sizeof(double);
		int sorSize = sizeof(source)/sizeof(double);
		if (tarSize >= length && sorSize >= length)
		{*/
			for (int i=0;i!=length;++i)
				target[i] += source[i];
		/*}
		else
		{
			throw "pdIntgManager::ArrayAddition: array too short";
		}*/
	}
	catch (char* str)
	{
		cerr << str;
		exit(0);
	}
}

bool pdIntgManager::IsConverge(const double pRun, const double cRun, const int endTrapIndex, int& endIntgPoint)
{
	/*

	This function judges if the error from two consecutive runs is less than the desired error.

	Input:
	  1. pRun: previous run result.
	  2. cRun: current run result.
	  3. endTrapIndex: the ending trapezoidal index.

    Output:
	  1. If the results converges, return true.
	  2. endIntgPoint: save the ending trapezoidal (integration) points into the '_intgPointEnd' array.

	*/

	double dss = cRun - pRun;
	double ss = pRun;
	// If the difference from two consecutive runs are greater than the tolerance "TOL"
	// compare the relative error with the desired accuracy "_eps"
	if (abs(ss) > TOL)
	{
		if (abs(dss) <= _eps * abs(ss))
		{
			endIntgPoint = int(pow(2, double(endTrapIndex-1))) + 1;
			return true;
		}
	}
	// If the difference from two consecutive runs are smaller than the tolerance "TOL"
	// compare the absolute error with the desired accuracy "_eps" (to prevent dead-lock)
	else
	{
		if (abs(dss) <= _eps)
		{
			endIntgPoint = int(pow(2, double(endTrapIndex-1))) + 1;
			return true;
		}
	}
	return false;
}

void pdIntgManager::InitMLSCoeff(const double x1, const double x2, const double x3, pdNode * node)
{
	/*

	This function calculate the MLS coefficient for the given trapezoidal (integration) points at
	(x1, x2, x3). Note this is the local coordinate. 

	*/

	// get the contributing node list of this node 
	vector<pdNode*> cbList;
	node->GetContribNodeList(cbList);
	// shift to global coordinate
	double x = x1 + _node_i->GetX1();
	double y = x2 + _node_i->GetX2();
	double z = x3 + _node_i->GetX3();
	int dim = POLYBASISDIM; // basis vector dimension
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
	node->AddMLSCoeff(coeff);
	// delete temporary objects	
	delete lhs;
	delete rhs;
	delete sol;
	delete b;
	delete b2;
}

void pdIntgManager::SetCBNodeDispl(double ** displ)
{
	_cbNodeDispl = displ;
}

void pdIntgManager::SetCurrentRunTime(const double time)
{
	_curRunTime = time;
}
