#include "pdIntgManager.h"

pdIntgManager::pdIntgManager(const double delta, const double dx, pdBond* bond, pdNode* node_i, pdNode* node_j)
: _delta(delta), _dx(dx), _bond(bond), _node_i(node_i), _node_j(node_j)
{
	_x1i = 0.0;
	_x2i = 0.0;
	_x3i = 0.0;
	_x1j = _node_j->GetX1() - _node_i->GetX1();
	_x2j = _node_j->GetX2() - _node_i->GetX2();
	_x3j = _node_j->GetX3() - _node_i->GetX3();
	_u1i = _node_i->GetU1();
	_u2i = _node_i->GetU2();
	_u3i = _node_i->GetU3();
	// turn off initMSLCoeff flag
	_initMLSCoeffOn = false;
	_cbNodeDispl = NULL;
}

pdIntgManager::pdIntgManager(const double delta, const double dx, pdBond* bond, const double* u, pdNode* node_i, pdNode* node_j)
: _delta(delta), _dx(dx), _bond(bond), _node_i(node_i), _node_j(node_j)
{
	_x1i = 0.0;
	_x2i = 0.0;
	_x3i = 0.0;
	_x1j = _node_j->GetX1() - _node_i->GetX1();
	_x2j = _node_j->GetX2() - _node_i->GetX2();
	_x3j = _node_j->GetX3() - _node_i->GetX3();
	// use trial displacement of node_i from RK4 steps for the force calculation
	_u1i = u[0];
	_u2i = u[1];
	_u3i = u[2];
	// turn off initMSLCoeff flag
	_initMLSCoeffOn = false;
	_cbNodeDispl = NULL;
}

pdIntgManager::~pdIntgManager()
{
	_abscis.clear();
	_weights.clear();
	if (_cbNodeDispl)
	{
		delete [] *_cbNodeDispl;
		delete [] _cbNodeDispl;
	}
}

void pdIntgManager::CalcIntg3D(string itg, double a, double b, vector<string> funcName, double* pdForce)
{
	/*

	This function calculates a three-dimensional integration for integrand selected by input string "itg".
	The integration order is not necessarily to be in "X-Y-Z" order.
	In this function, the trapzoidal points in X (outer) direction are generated.
	There are totally three loops and this is the outer loop of the integration.

	Input:
	  1. itg: name of the integrand.
	  2. a: lower integration limit for X (outer) direction.
	  3. b: upper integration limit for X (outer) direction.
	  4. funcName[0]: lower integration limits function name for the middle direction.
      5. funcName[1]: upper integration limits function name for the middle direction.
      6. funcName[2]: lower integration limits function name for the inner direction.
      7. funcName[3]: upper integration limits function name for the inner direction.
      8. funcName[4]: indicates the integration order (from inner to outer).

    Output:
	  1. the final integration result.

     */

	_intgPointCounter = 0;
	_intgPointEnd[0] = 0;
	// Start the integration loop
	if (abs(a-b)<=1.0e-10) // the upper and lower integration limits for the outer direction are equal
	{
		pdForce[0] = 0.0;
		pdForce[1] = 0.0;
		pdForce[2] = 0.0;
	}
	else
	{
		_nrfunc = SelectIntegrand(itg); // save the integrand
		if (_errorControlOn) // Adaptive trapezoidal integration with error control
		{
			double s[JMAXP]; // store the results from each run
			for (int j = _intgPoint; j <= JMAX; ++j)
			{
				CalcTrapzdQuad(&pdIntgManager::CalcIntgX1, a, b, funcName, j, 1, pdForce);
				s[j] = pdForce[0];
				if (j >= (_intgPoint + 1))
				{
					if (j >= SMAX)
					{
						_intgPointEnd[1] = int(pow(2, double(j-1))) + 1;
						pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
						break;
					}
					else
					{
						if (IsConverge(s[j-1], s[j], j, _intgPointEnd[1]))
						{
							// result converges
							pdForce[0] = s[j];
							break;
						}
					}
				}
			}
		}
		else
		{
			// Adaptive trapezoidal integration without error control, grid is not refined
			if (_aiOn)
			{
				_intgPointEnd[1] = int(pow(2, double(_intgPoint-1))) + 1;
				CalcTrapzdQuad(&pdIntgManager::CalcIntgX1, a, b, funcName, _intgPoint, 1, pdForce);
			}
			// Fixed Gaussian integration
			else
			{
				CalcGaussianQuad(&pdIntgManager::CalcIntgX1, a, b, funcName, _intgPoint, 1, pdForce);
			}

		}
	}
}

void pdIntgManager::CalcIntgX1(double x, vector<string> funcName, double* pdForce)
{
	/*

	This function generates the trapzoidal points in Y (middle) direction of a three-dimensional integration.
	This is the middle loop of the integration.

	Input:
	  1. x: value of the trapezoidal point for X (outer) direction.
	  4. funcName[0]: lower integration limits function name for the middle direction.
      5. funcName[1]: upper integration limits function name for the middle direction.
      6. funcName[2]: lower integration limits function name for the inner direction.
      7. funcName[3]: upper integration limits function name for the inner direction.
      8. funcName[4]: indicates the integration order (from inner to outer).

    Output:
	  1. the value of the integrand with trapezoidal points from three loops.

     */

	_xSave = x;
	double a = SelectYLimitFuncValue(x, funcName[0]); // lower Y (middle direction) integration limit
	double b = SelectYLimitFuncValue(x, funcName[1]); // upper Y (middle direction) integration limit

	if (abs(a-b)<=1.0e-10) // If the upper and lower integration limits are equal, return 0
	{
		pdForce[0] = 0.0;
		pdForce[1] = 0.0;
		pdForce[2] = 0.0;
	}
	else
	{
		if (_errorControlOn) // Adaptive trapezoidal integration with error control
		{
			double s[JMAXP]; // store the results from each run
			for (int j=_intgPoint;j<=JMAX;++j)
			{
				CalcTrapzdQuad(&pdIntgManager::CalcIntgX2, a, b, funcName, j, 2, pdForce);
				s[j] = pdForce[1];
				if (j>=(_intgPoint+1))
				{
					if (j >= SMAX)
					{
						_intgPointEnd[2] = int(pow(2, double(j-1))) + 1;
						pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
						break;
					}
					else
					{
						if (IsConverge(s[j-1], s[j], j, _intgPointEnd[2]))
						{
							// result converges
							pdForce[1] = s[j];
							break;
						}
					}
				}
			}
		}
		else
		{
			// Adaptive trapezoidal integration without error control, grid is not refined
			if (_aiOn)
			{
				_intgPointEnd[2] = int(pow(2, double(_intgPoint-1))) + 1;
				CalcTrapzdQuad(&pdIntgManager::CalcIntgX2, a, b, funcName, _intgPoint, 2, pdForce);
			}
			// Fixed Gaussian integration
			else
			{
				CalcGaussianQuad(&pdIntgManager::CalcIntgX2, a, b, funcName, _intgPoint, 2, pdForce);
			}
		}
	}
}

void pdIntgManager::CalcIntgX2(double y, vector<string> funcName, double* pdForce)
{
	/*

	This function generates the trapzoidal points in Z (inner) direction of a three-dimensional integration.
	This is the inner loop of the integration.

	Input:
	  1. y: value of the trapezoidal value for Y (middle) direction.
	  4. funcName[0]: lower integration limits function name for the middle direction.
      5. funcName[1]: upper integration limits function name for the middle direction.
      6. funcName[2]: lower integration limits function name for the inner direction.
      7. funcName[3]: upper integration limits function name for the inner direction.
      8. funcName[4]: indicates the integration order (from inner to outer).

    Output:
	  1. the value of the integrand with trapezoidal points from three loops.

     */

	_ySave = y;
	double a = SelectZLimitFuncValue(_xSave, y, funcName[2]); // lower Z (inner direction) integration limit
	double b = SelectZLimitFuncValue(_xSave, y, funcName[3]); // upper Z (inner direction) integration limit

	if (abs(a-b)<=1.0e-10) // If the upper and lower integration limits are equal, return 0
	{
		pdForce[0] = 0.0;
		pdForce[1] = 0.0;
		pdForce[2] = 0.0;
	}
	else
	{
		if (_errorControlOn) // Adaptive trapezoidal integration with error control
		{
			double s[JMAXP]; // store the results from each run
			for (int j=_intgPoint;j<=JMAX;++j)
			{
				CalcTrapzdQuad(&pdIntgManager::CalcIntgX3, a, b, funcName, j, 3, pdForce);
				s[j] = pdForce[2];
				if (j>=(_intgPoint+1))
				{
					if (j >= SMAX)
					{
						_intgPointEnd[3] = int(pow(2, double(j-1))) + 1;
						pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
						break;
					}
					else
					{
						if (IsConverge(s[j-1], s[j], j, _intgPointEnd[3]))
						{
							// result converges
							pdForce[2] = s[j];
							break;
						}
					}
				}
			}
		}
		else
		{
			// Adaptive trapezoidal integration without error control, grid is not refined
			if (_aiOn)
			{
				_intgPointEnd[3] = int(pow(2, double(_intgPoint-1))) + 1;
				CalcTrapzdQuad(&pdIntgManager::CalcIntgX3, a, b, funcName, _intgPoint, 3, pdForce);
			}
			// Fixed Gaussian integration
			else
			{
				CalcGaussianQuad(&pdIntgManager::CalcIntgX3, a, b, funcName, _intgPoint, 3, pdForce);
			}
		}
	}
}

void pdIntgManager::CalcIntgX3(double z, vector<string> funcName, double* pdForce)
{
	/*

	This function does the following:
	  1. call the integrand and calculates values based on the integration points.
      2. re-arrange the corrdinates of integration points to the right integration order (important to AI method!).
	  3. generate the value of trapezoidal points if _initTrapPointOn flag is true.

	Input:
	  1. z: value of the trapezoidal value for Z (inner) direction.
	  4. funcName[0]: lower integration limits function name for the middle direction.
      5. funcName[1]: upper integration limits function name for the middle direction.
      6. funcName[2]: lower integration limits function name for the inner direction.
      7. funcName[3]: upper integration limits function name for the inner direction.
      8. funcName[4]: indicates the integration order (from inner to outer).

    Output:
	  1. the value of the integrand with integration points from three loops.

     */

    string swapXYZ = funcName[4];
	double xFunc, yFunc, zFunc;
    if (swapXYZ == "XZY")
	{
		xFunc = _xSave;
		yFunc = z;
		zFunc = _ySave;        
	}
    else if (swapXYZ == "ZXY")
	{
		xFunc = _ySave;
		yFunc = z;
		zFunc = _xSave; 
	}
    else if (swapXYZ == "ZYX")
	{
		xFunc = z;
		yFunc = _ySave;
		zFunc = _xSave;        
	}
    else if (swapXYZ == "YZX")
	{
		xFunc = z;
		yFunc = _xSave;
		zFunc = _ySave; 
	}
    else if (swapXYZ == "YXZ")
	{
		xFunc = _ySave;
		yFunc = _xSave;
		zFunc = z; 
	}
	else // swapXYZ == "XYZ"
	{
		xFunc = _xSave;
		yFunc = _ySave;
		zFunc = z; 
	}
	(this->*_nrfunc)(xFunc, yFunc, zFunc, pdForce);
	// increase integration point index by 1
	++_intgPointCounter;
}

void pdIntgManager::CalcTrapzdQuad(FUNC_DV func, double a, double b, vector<string> funcName, int n, int flag, double* pdForce)
{
	/*

	This function computes the nth stage of refinement of an extended trapezoidal rule where n is the trapezoidal index.
	"func" is a pointer to the function to be integrated between limits of a and b.
	When called with n=1, the routine returns the crudest estimate. Subsequent calls with n=2, 3...(in that sequential order)
	will improve the accuracy by adding 2^(n-2) additional interior points.\
	Note: case 2, 3 are required if the initial trapezoidal index is not 1.

	Input:
	  1. func: function pointer to the integrand.
	  2. a, b: lower and upper integration limits.
	  3. funcName: container of function names for integration limits in other two directions.
	  4. n: trapezoidal index. The total number of trapezoidal (integration) points = 2^{n-1}+1 (including two end points of interval).
	  5. flag: indicates the integration direction the integration result belongs to. flag = 1, 2, 3.

    Output:
      1. pdForce: value of the integrand.

    */

	static double s0[4];
	static double s1[4];
	static double s2[4];

	switch (n)
	{
	case 1: // just two end points of the integration interval (2 trapezoidal points)
		{
			double pdForceTemp1[3];
			double pdForceTemp2[3];
			(this->*func)(a, funcName, pdForceTemp1);
			(this->*func)(b, funcName, pdForceTemp2);
			s0[flag] = (b-a)*(0.5*pdForceTemp1[0] + 0.5*pdForceTemp2[0]);
			s1[flag] = (b-a)*(0.5*pdForceTemp1[1] + 0.5*pdForceTemp2[1]);
			s2[flag] = (b-a)*(0.5*pdForceTemp1[2] + 0.5*pdForceTemp2[2]);
			pdForce[0] = s0[flag];
			pdForce[1] = s1[flag];
			pdForce[2] = s2[flag];
			break;
		}
	case 2: // divide the interval into two segments (3 trapezoidal points)
		{
			double pdForceTemp1[3];
			double pdForceTemp2[3];
			double pdForceTemp3[3];
			(this->*func)(a, funcName, pdForceTemp1);
			(this->*func)(0.5*(a+b), funcName, pdForceTemp2);
			(this->*func)(b, funcName, pdForceTemp3);
			s0[flag] = 0.5*(b-a)*(0.5*pdForceTemp1[0] + pdForceTemp2[0] + 0.5*pdForceTemp3[0]);
			s1[flag] = 0.5*(b-a)*(0.5*pdForceTemp1[1] + pdForceTemp2[1] + 0.5*pdForceTemp3[1]);
			s2[flag] = 0.5*(b-a)*(0.5*pdForceTemp1[2] + pdForceTemp2[2] + 0.5*pdForceTemp3[2]);
			pdForce[0] = s0[flag];
			pdForce[1] = s1[flag];
			pdForce[2] = s2[flag];
			break;
		}
	case 3: // divide the interval into 4 segments (5 trapezoidal points)
		{
			double pdForceTemp1[3];
			double pdForceTemp2[3];
			double pdForceTemp3[3];
			double pdForceTemp4[3];
			double pdForceTemp5[3];
			(this->*func)(a, funcName, pdForceTemp1);
			(this->*func)(0.25*(3.0*a+b), funcName, pdForceTemp2);
			(this->*func)(0.5*(a+b), funcName, pdForceTemp3);
			(this->*func)(0.25*(a+3.0*b), funcName, pdForceTemp4);
			(this->*func)(b, funcName, pdForceTemp5);
			s0[flag] = 0.25*(b-a)*(0.5*pdForceTemp1[0] + pdForceTemp2[0] + pdForceTemp3[0] + pdForceTemp4[0] + 0.5*pdForceTemp5[0]);
			s1[flag] = 0.25*(b-a)*(0.5*pdForceTemp1[1] + pdForceTemp2[1] + pdForceTemp3[1] + pdForceTemp4[1] + 0.5*pdForceTemp5[1]);
			s2[flag] = 0.25*(b-a)*(0.5*pdForceTemp1[2] + pdForceTemp2[2] + pdForceTemp3[2] + pdForceTemp4[2] + 0.5*pdForceTemp5[2]);
			pdForce[0] = s0[flag];
			pdForce[1] = s1[flag];
			pdForce[2] = s2[flag];
			break;
		}
	default:
		{
			double x, tnm, sum0, sum1, sum2, del;
			int it, j;
			for (it=1, j=1;j<n-1;j++)
				it <<= 1;
			tnm = it;
			del = (b-a)/tnm; // this is the spacing of the points to be added
			x = a + 0.5*del;
			for (sum0=0.0, sum1=0.0, sum2=0.0,j=1;j<=it;j++, x+=del)
			{
				(this->*func)(x, funcName, pdForce);
				sum0 += pdForce[0];
				sum1 += pdForce[1];
				sum2 += pdForce[2];
			}
			s0[flag] = 0.5*(s0[flag]+(b-a)*sum0/tnm); // this replaces s by its refined value
			s1[flag] = 0.5*(s1[flag]+(b-a)*sum1/tnm);
			s2[flag] = 0.5*(s2[flag]+(b-a)*sum2/tnm);
			// save the results
			pdForce[0] = s0[flag];
			pdForce[1] = s1[flag];
			pdForce[2] = s2[flag];
			break;
		}
	}
}

void pdIntgManager::CalcGaussianQuad(FUNC_DV func, double a, double b, vector<string> funcName, int n, int flag, double* pdForce)
{
	/*

	This function computes the nth stage of refinement using Gaussian-Legendre quadrature where n is the number of Gaussian point.
	"func" is a pointer to the function to be integrated between limits of a and b.
	The location of Gaussian points passed to the integrand is in family (local) coordinate and fixed.

	Input:
	  1. func: function pointer to the integrand.
	  2. a, b: lower and upper integration limits.
	  3. funcName: container of function names for integration limits in other two directions.
	  4. n: number of Gaussian points.
	  5. flag: indicates the integration direction the integration result belongs to. flag = 1, 2, 3.
	  6. pdForce: container to save the peridynamic force.

    */

	static double s0[4];
	static double s1[4];
	static double s2[4];
	s0[flag] = 0.0;
	s1[flag] = 0.0;
	s2[flag] = 0.0;
	// calculate the function value
	for (int j=0;j!=_intgPoint;++j)
	{
		double dx = 0.5*_dx*_abscis[j];
		// select xCenter based on the integration order stored in funcName[4]
		double xCenter;
		if (funcName[4][flag-1] == 'X')
		{
			xCenter = _x1j;// family coordinates
		}
		else if (funcName[4][flag-1] == 'Y')
		{
			xCenter = _x2j;
		}
		else
		{
			xCenter = _x3j;
		}
		(this->*func)(dx + xCenter, funcName, pdForce);
		// use fixed Gaussian points		
		s0[flag] += _weights[j]*pdForce[0];
		s1[flag] += _weights[j]*pdForce[1];
		s2[flag] += _weights[j]*pdForce[2];
	}
	// scale the answer to the integration interval
	s0[flag] *= 0.5*(b-a);
	s1[flag] *= 0.5*(b-a); 
	s2[flag] *= 0.5*(b-a); 
	// save the results
	pdForce[0] = s0[flag];
	pdForce[1] = s1[flag];
	pdForce[2] = s2[flag];
}

double pdIntgManager::SelectYLimitFuncValue(double x, string funcName)
{
	FUNC_S func = SelectYLimitFunc(funcName);
	return (this->*func)(x);
}

double pdIntgManager::SelectZLimitFuncValue(double x, double y, string funcName)
{
	FUNC_D func = SelectZLimitFunc(funcName);
	return (this->*func)(x, y);
}

void pdIntgManager::AdaptiveIntegration(string intergrand, double* pdForce)
{
	// the geometric configurations of node i and j fall into two categoeries
    // char intergrand is the name of the intergrand function

    // type 1 (node j on one of the axis)
    if (_x1j != _x1i && _x2j == _x2i && _x3j == _x3i) // node i and j are sitting on the same axis (X)
    {
        CalcConfigurationType1(intergrand, pdForce);
    }
    else if (_x1j == _x1i && _x2j != _x2i && _x3j == _x3i) // node i and j are sitting on the same axis (Y)
    {
        CalcConfigurationType1(intergrand, pdForce);
    }
    else if (_x1j == _x1i && _x2j == _x2i && _x3j != _x3i) // node i and j are sitting on the same axis (Z)
    {
        CalcConfigurationType1(intergrand, pdForce);
    }
    // type 2
    else
    {
        CalcConfigurationType2(intergrand, pdForce);
    }
}

void pdIntgManager::CalcConfigurationType1(string itg, double* pdForce)
{
    // this function calculates integration for geometric configuration type 1: two nodes on the same axis

    // first find the intersection points with eight edges of node j's cell
    double iPoints[8]; // store the coordinate of intersectinon points
    vector<string> funcName;
    funcName.push_back("YLimitXLow");
    funcName.push_back("YLimitXHigh");
    funcName.push_back("YLimitYLow");
    funcName.push_back("YLimitYHigh");
    funcName.push_back("YLimitZLow");
    funcName.push_back("YLimitZHigh");
    int b = this->FindIntersectionPointOnEightEdges(funcName, iPoints); // binary flag

    double dist = abs(_x1j - _x1i) + abs(_x2j - _x2i) + abs(_x3j - _x3i); // distant between two nodes
    // distance limit for each subtype (from large to small)
    double d1 = _delta + 0.5 * _dx;
    double d2 = sqrt(_delta * _delta - 0.25 * _dx * _dx) + 0.5 * _dx;
    double d3 = sqrt(_delta * _delta - 0.5 * _dx * _dx) + 0.5 * _dx;
    double d4 = _delta - 0.5 * _dx;
    double d5 = sqrt(_delta * _delta - 0.25 * _dx * _dx) - 0.5 * _dx;
    double d6 = sqrt(_delta * _delta - 0.5 * _dx * _dx) - 0.5 * _dx;

    // find the integration limtis and carry out the integration
    // subtype 1: d2<dist<d1
    if (dist >= d2 && dist <= d1)
    {
        if (_x1j != 0.0) // subtype1, node j on X axis
            CalcConfigurationType1Sub1("Z", "Y", "X", itg, dist, pdForce);
        else if (_x2j != 0.0) // subtype1, node j on Y axis
            CalcConfigurationType1Sub1("Z", "X", "Y", itg, dist, pdForce);
        else // (_x3j != 0.0)  subtype1, node j on Z axis
            CalcConfigurationType1Sub1("X", "Y", "Z", itg, dist, pdForce);
    }
    // subtype 2: d3<dist<d2
    else if (dist >= d3 && dist < d2)
    {
        switch (b)
        {
            case 17: // horizon sphere intersects with T1, B1 (actually two intersection points on both edges), node j on +X axis
				{
					CalcConfigurationType1Sub2("Z", "Y", "X", itg, dist, pdForce);
					break;
				}
            case 34: // horizon sphere intersects with T3, B3 (actually two intersection points on both edges), node j on -X axis
				{
					CalcConfigurationType1Sub2("Z", "Y", "X", itg, dist, pdForce);
					break;
				}
            case 68: // horizon sphere intersects with T2, B2 (actually two intersection points on both edges), node j on +Y axis
                {
					CalcConfigurationType1Sub2("Z", "X", "Y", itg, dist, pdForce);
					break;
				}
            case 136: // horizon sphere intersects with T4, B4 (actually two intersection points on both edges), node j on -Y axis
				{
					CalcConfigurationType1Sub2("Z", "X", "Y", itg, dist, pdForce);
					break;
				}
            case 240: // horizon sphere intersects with B1-B4 (actually two intersection points on each four edges), node j on +Z axis
				{
					CalcConfigurationType1Sub2("Y", "X", "Z", itg, dist, pdForce);
					break;
				}
            case 15: // horizon sphere intersects with T1-T4 (actually two intersection points on each four edges), node j on -Z axis
				{
					CalcConfigurationType1Sub2("Y", "X", "Z", itg, dist, pdForce);
					break;
				}
            default:
				{
					pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
					break;
				}
        }
    }
    // subtype 3: d4<dist<d3
    else if (dist >= d4 && dist < d3)
    {
        switch (b)
        {
            case 204: // horizon sphere intersects T2, T4, B2, B4, node j on X axis
				{
					CalcConfigurationType1Sub3("Y", "Z", "X", itg, pdForce);
					break;
				}
            case 51: // horizon sphere intersects T1, T3, B1, B3, node j on Y axis
				{
					CalcConfigurationType1Sub3("X", "Z", "Y", itg, pdForce);
					break;
				}
            case 0: // horizon sphere intersects V1-V4, node j on Z axis
				{
					CalcConfigurationType1Sub3("X", "Y", "Z", itg, pdForce);
					break;
				}
            default:
                {
					pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
					break;
				}
        }
    }
    // subtype 4: d5<dist<d4
    else if (dist >= d5 && dist < d4)
    {
        switch (b)
        {
            case 204: // horizon sphere intersects T2, T4, B2, B4, node j on X axis
				{
					CalcConfigurationType1Sub4("Z", "Y", "X", itg, dist, pdForce);
					break;
				}
            case 51: // horizon sphere intersects T1, T3, B1, B3, node j on Y axis
				{
					CalcConfigurationType1Sub4("Z", "X", "Y", itg, dist, pdForce);
					break;
				}
            case 0: // horizon sphere intersects V1-V4, node j on Z axis
				{
					CalcConfigurationType1Sub4("Y", "X", "Z", itg, dist, pdForce);
					break;
				}
            default:
                {
					pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
					break;
				}
        }
    }
    // subtype 5: d6<dist<d5
    else if (dist >= d6 && dist < d5)
    {
        switch (b)
        {
            case 238: // horizon sphere intersects T2, T3, T4, B2, B3, B4 (actually two intersection points on T3&B3), node j on the +X axis
				{
					CalcConfigurationType1Sub5("Z", "Y", "X", itg, dist, pdForce);
					break;
				}
            case 221: // horizon sphere intersects T1, T2, T4, B1, B2, B4 (actually two intersection points on T1&B1), node j on the -X axis
				{
					CalcConfigurationType1Sub5("Z", "Y", "X", itg, dist, pdForce);
					break;
				}
            case 187: // horizon sphere intersects T1, T3, T4, B1, B3, B4 (actually two intersection points on T4&B4), node j on the +Y axis
				{
					CalcConfigurationType1Sub5("Z", "X", "Y", itg, dist, pdForce);
					break;
				}
            case 119: // horizon sphere intersects T1, T2, T3, B1, B2, B3 (actually two intersection points on T2&B2), node j on the -Y axis
				{
					CalcConfigurationType1Sub5("Z", "X", "Y", itg, dist, pdForce);
					break;
				}
            case 15: // horizon sphere intersects with T1-T4 (actually two intersection points on each four edges), node j on the +Z axis
				{
					CalcConfigurationType1Sub5("Y", "X", "Z", itg, dist, pdForce);
					break;
				}
            case 240: // horizon sphere intersects with B1-B4 (actually two intersection points on each four edges), node j on the -Z axis
				{
					CalcConfigurationType1Sub5("Y", "X", "Z", itg, dist, pdForce);
					break;
				}
            default:
                {
					pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
					break;
				}
        }
    }
    else
    {
        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
    }
}

void pdIntgManager::CalcConfigurationType2(string itg, double* pdForce)
{
    // this function calculates integration for geometric configuration type 2: node j not on any of the axis
    // node i is the source node, and node j is the field node

    // first find the intersection points with eight edges of node j's cell
    double iPoints[8]; // store the coordinate of intersectinon points
    vector<string> funcName;
    funcName.push_back("YLimitXLow");
    funcName.push_back("YLimitXHigh");
    funcName.push_back("YLimitYLow");
    funcName.push_back("YLimitYHigh");
    funcName.push_back("YLimitZLow");
    funcName.push_back("YLimitZHigh");
    int b = this->FindIntersectionPointOnEightEdges(funcName, iPoints); // binary flag

    // find the integration limtis and carry out the integration
    vector<string> flag; // store the integration limit function names
    // the binary b corresponding to the 8 edges of top and bottom surfaces: B4 B2 B3 B1 | T4 T2 T3 T1 (read binary b from right to left)
	double pdForceStep[3];
	pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
    switch (b)
    {
        case 85: // b=0101|0101, horizon sphere intersects with T1, T2, B1, B2
            {
                if (_x1j > 0 && _x2j > 0) // node j in +X +Y +/-Z octant
                    CalcConfigurationType2Sub1("Z", "X", "Y", itg, pdForce);
                else if (_x1j < 0 && _x2j < 0) // node j on -X-Y plane
                    CalcConfigurationType2Sub2("Z", "X", "Y", itg, pdForce);
                else
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
				break;
            }
        case 102: // b=0110|0110, horizon sphere intersects with T2, T3, B2, B3
            {
                if (_x1j < 0 && _x2j > 0)  // node j in  -X +Y +/-Z octant
                    CalcConfigurationType2Sub1("Z", "X", "Y", itg, pdForce);
                else if (_x1j > 0 && _x2j < 0) // node j on +X-Y plane
                    CalcConfigurationType2Sub2("Z", "X", "Y", itg, pdForce);
                else
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
				break;
            }
        case 153: // b=1001|1001, horizon sphere intersects with T1, T4, B1, B4
            {
                if (_x1j > 0 && _x2j < 0)  // node j in +X -Y +/-Z octant
                    CalcConfigurationType2Sub1("Z", "X", "Y", itg, pdForce);
                else if (_x1j < 0 && _x2j > 0) // node j on -X+Y plane
                    CalcConfigurationType2Sub2("Z", "X", "Y", itg, pdForce);
                else
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
				break;
            }
        case 170: // b=1001|1001, horizon sphere intersects with T3, T4, B3, B4
            {
                if (_x1j < 0 && _x2j < 0)   // node j is in -X -Y +/-Z octant
                    CalcConfigurationType2Sub1("Z", "X", "Y", itg, pdForce);
                else if (_x1j > 0 && _x2j > 0)  // node j on +X+Y plane
                    CalcConfigurationType2Sub2("Z", "X", "Y", itg, pdForce);
                else
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
				break;
            }
        case 12: // b=0000|1100, horizon sphere intersects with T2, T4 only
            {
                // check the intersection points with the other four edges
                double nPoints[4]; // store the coordinate of intersectinon points
                int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                if (nb == 6) // horizon sphere also intersects with V2, V3
                {
                    if (_x1j > 0 && _x3j > 0) // node j on +X+Z plane
                        CalcConfigurationType2Sub2("Y", "X", "Z", itg, pdForce);
                    else if (_x1j < 0 && _x3j < 0)  // node j on -X+/-Y-Z octant
                        CalcConfigurationType2Sub1("Y", "X", "Z", itg, pdForce);
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else if (nb == 9) // horizon sphere also intersects with V1, V4
                {
                    if (_x1j < 0 && _x3j > 0) // node j on -X+Z plane
                        CalcConfigurationType2Sub2("Y", "X", "Z", itg, pdForce);
                    else if (_x1j > 0 && _x3j < 0)  // node j on +X+/-Y-Z octant
                        CalcConfigurationType2Sub1("Y", "X", "Z", itg, pdForce);
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
				break;
            }
        case 192: // b=1100|0000, horizon sphere intersects with B2, B4 only
            {
                // check the intersection points with the other four edges
                double nPoints[4]; // store the coordinate of intersectinon points
                int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                if (nb == 6) // horizon sphere also intersects with V2, V3
                {
                    if (_x1j > 0 && _x3j < 0)  // node j on +X-Z plane
                        CalcConfigurationType2Sub2("Y", "X", "Z", itg, pdForce);
                    else if (_x1j < 0 && _x3j > 0)  // node j on -X+/-Y+Z octant
                        CalcConfigurationType2Sub1("Y", "X", "Z", itg, pdForce);
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else if (nb == 9) // horizon sphere also intersects with V1, V4
                {
                    if (_x1j < 0 && _x3j < 0)  // node j on -X-Z plane
                        CalcConfigurationType2Sub2("Y", "X", "Z", itg, pdForce);
                    else if (_x1j > 0 && _x3j > 0)  // node j on +X+/-Y+Z octant
                        CalcConfigurationType2Sub1("Y", "X", "Z", itg, pdForce);
                      else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else
                      pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
				break;
            }
        case 3: // b=0000|0011, horizon sphere intersects with T1, T3 only
            {
                // check the intersection points with the other four edges
                double nPoints[4]; // store the coordinate of intersectinon points
                int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                if (nb == 12) // horizon sphere also intersects with V3, V4
                {
                    if (_x2j > 0 && _x3j > 0) // node j on +Y+Z plane
                        CalcConfigurationType2Sub2("X", "Z", "Y", itg, pdForce);
                    else if (_x2j < 0 && _x3j < 0) // node j on +/-X-Y-Z octant
                        CalcConfigurationType2Sub1("X", "Z", "Y", itg, pdForce);
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else if (nb == 3) // horizon sphere also intersects with V1, V2
                {
                    if (_x2j < 0 && _x3j > 0)  // node j on -Y+Z plane
                        CalcConfigurationType2Sub2("X", "Z", "Y", itg, pdForce);
                    else if (_x2j > 0 && _x3j < 0)  // node j on +/-X+Y-Z octant
                        CalcConfigurationType2Sub1("X", "Z", "Y", itg, pdForce);
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
				break;
            }
        case 48: // b=0011|0000, horizon sphere intersects with B1, B3 only
            {
                // check the intersection points with the other four edges
                double nPoints[4]; // store the coordinate of intersectinon points
                int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                if (nb == 12) // horizon sphere also intersects with V3, V4
                {
                    if (_x2j > 0 && _x3j < 0)  // node j on +Y-Z plane
                        CalcConfigurationType2Sub2("X", "Z", "Y", itg, pdForce);
                    else if (_x2j < 0 && _x3j > 0)  // node j on +/-X-Y+Z ocatnt
                        CalcConfigurationType2Sub1("X", "Z", "Y", itg, pdForce);
                      else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else if (nb == 3) // horizon sphere also intersects with V1, V2
                {
                    if (_x2j < 0 && _x3j < 0)  // node j on -Y-Z plane
                        CalcConfigurationType2Sub2("X", "Z", "Y", itg, pdForce);
                    else if (_x2j > 0 && _x3j > 0)  // node j on +/-X+Y+Z octant
                        CalcConfigurationType2Sub1("X", "Z", "Y", itg, pdForce);
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
				break;
            }
        case 165: // b=1010|0101, horizon sphere intersects with T1, T2, B4, B3
            {
                if (iPoints[2] < iPoints[7])
                {
                    if (_x1j > 0 && _x2j > 0 && _x3j > 0) // node j in +X +Y +Z octant
                    {
                        // integrate Z first, then project to X-Y
                        flag.push_back("YLimitHighYXPos");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[2], iPoints[7], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitLowYXPos");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitHighYXPos");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j < 0 && _x2j < 0 && _x3j < 0) // node j in -X -Y -Z octant
                    {
                        // integrate Z first, then project to X-Y
                        flag.push_back("YLimitHighYXNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[2], iPoints[7], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitLowYXNeg");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowYXNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else // (iPoints[2] >= iPoints[7])
                {
                    if (_x1j > 0 && _x2j > 0 && _x3j > 0) // node j in +X +Y +Z octant
                    {
                        // integrate Z first, then project to X-Y
                        flag.push_back("YLimitHighYXPos");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[7], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitHighYXPos");
                        flag.push_back("YLimitLowYXPos");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[7], iPoints[2], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitLowYXPos");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[2], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitHighYXPos");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j < 0 && _x2j < 0 && _x3j < 0) // node j in -X -Y -Z octant
                    {
                        // integrate Z first, then project to X-Y
                        flag.push_back("YLimitHighYXNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[7], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitHighYXNeg");
                        flag.push_back("YLimitLowYXNeg");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[7], iPoints[2], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitLowYXNeg");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[2], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowYXNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
				break;
            }
        case 90: // b=0101|1010, horizon sphere intersects with T3, T4, B1, B2
            {
                if (iPoints[3] < iPoints[6])
                {
                    if (_x1j > 0 && _x2j > 0 && _x3j < 0) // node j in +X +Y -Z octant
                    {
                        // integrate Y first, then project to X-Z
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitHighZXNeg");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[3], iPoints[6], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowZXNeg");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitHighZXNeg");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j < 0 && _x2j < 0 && _x3j > 0) // node j in -X -Y +Z octant
                    {
                        // integrate Y first, then project to X-Z
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitHighZXPos");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[3], iPoints[6], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowZXPos");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitLowZXPos");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else // (iPoints[3] >= iPoints[6])
                {
                    if (_x1j > 0 && _x2j > 0 && _x3j < 0) // node j in +X +Y -Z octant
                    {
                        // integrate Y first, then project to X-Z
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitHighZXNeg");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[6], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowZXNeg");
                        flag.push_back("YLimitHighZXNeg");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[6], iPoints[3], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowZXNeg");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[3], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitHighZXNeg");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j < 0 && _x2j < 0 && _x3j > 0) // node j in -X -Y +Z octant
                    {
                        // integrate Y first, then project to X-Z
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitHighZXPos");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[6], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowZXPos");
                        flag.push_back("YLimitHighZXPos");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[6], iPoints[3], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowZXPos");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[3], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitLowZXPos");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
				break;
            }
        case 150: // b=1001|0110, horizon sphere intersects with T2, T3, B1, B4
            {
                if (iPoints[2] < iPoints[7])
                {
                    if (_x1j < 0 && _x2j > 0 && _x3j > 0)
                    {
                        // node j in -X +Y +Z octant
                        flag.push_back("YLimitHighZXPos");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[2], iPoints[7], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitLowZXPos");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitHighZXPos");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j > 0 && _x2j < 0 && _x3j < 0)
                    {
                        // node j in +X -Y -Z octant
                        flag.push_back("YLimitHighZXNeg");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[2], iPoints[7], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitLowZXNeg");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowZXNeg");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else // (iPoints[2] >= iPoints[7])
                {
                    if (_x1j < 0 && _x2j > 0 && _x3j > 0)
                    {
                        // node j in -X +Y +Z octant
                        flag.push_back("YLimitHighZXPos");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[2], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitHighZXPos");
                        flag.push_back("YLimitLowZXPos");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[7], iPoints[2], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitLowZXPos");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[7], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitHighZXPos");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j > 0 && _x2j < 0 && _x3j < 0)
                    {
                        // node j in +X -Y -Z octant
                        flag.push_back("YLimitHighZXNeg");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[2], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitHighZXNeg");
                        flag.push_back("YLimitLowZXNeg");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[7], iPoints[2], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitLowZXNeg");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[7], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowZXNeg");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
				break;
            }
        case 105: // b=0110|1001, horizon sphere intersects with T1, T4, B2, B3
            {
                if (iPoints[3] < iPoints[6])
                {
                    if (_x1j < 0 && _x2j > 0 && _x3j < 0) // node j in -X +Y -Z octant
                    {
                        // integrate Z first, then project to X-Y
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitHighYXPos");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[3], iPoints[6], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowYXPos");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitLowYXPos");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j > 0 && _x2j < 0 && _x3j > 0) // node j in +X -Y +Z octant
                    {
                        // integrate Z first, then project to X-Y
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitHighYXNeg");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[3], iPoints[6], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowYXNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitHighYXNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else // (iPoints[3] >= iPoints[6])
                {
                    if (_x1j < 0 && _x2j > 0 && _x3j < 0) // node j in -X +Y -Z octant
                    {
                        // integrate Z first, then project to X-Y
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitHighYXPos");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[6], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowYXPos");
                        flag.push_back("YLimitHighYXPos");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[6], iPoints[3], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowYXPos");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[3], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitLowYXPos");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j > 0 && _x2j < 0 && _x3j > 0) // node j in +X -Y +Z octant
                    {
                        // integrate Z first, then project to X-Y
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitHighYXNeg");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[6], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowYXNeg");
                        flag.push_back("YLimitHighYXNeg");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[6], iPoints[3], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowYXNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[3], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitHighYXNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
				break;
            }
        case 10: // b=0000|1010, horizon sphere intersects with T4, T3 only
            {
                // check the intersection points with the other four edges
                double nPoints[4]; // store the coordinate of intersectinon points
                int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                if (nb == 4) // horizon sphere also intersects with V3
                {
                    if (_x1j > 0 && _x2j > 0 && _x3j > 0)
                    {
                        // node j in +X+Y+Z octant
                        flag.push_back("YLimitHighZXPos");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[3], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitHighZXPos");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[3], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j < 0 && _x2j < 0 && _x3j < 0)
                    {
                        // node j in -X-Y-Z octant
                        flag.push_back("YLimitHighZXNeg");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[3], _x1j + 0.5 * _dx, flag, pdForce);
                    }
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else if (nb == 11) // horizon sphere also intersects with V1, V2, V4
                {
                    if (_x1j > 0 && _x2j > 0 && _x3j < 0)
                    {
                        // node j in +X+Y-Z octant
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitHighYXPos");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[3], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j < 0 && _x2j < 0 && _x3j > 0)
                    {
                        // node j in -X-Y+Z octant
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitHighYXNeg");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[3], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitHighYXNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[3], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else
                {
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
				break;
            }
        case 160: // b=1010|0000, horizon sphere intersects with B4, B3 only
            {
                // check the intersection points with the other four edges
                double nPoints[4]; // store the coordinate of intersectinon points
                int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                if (nb == 4) // horizon sphere also intersects with V3
                {
                    if (_x1j > 0 && _x2j > 0 && _x3j < 0)
                    {
                        // node j in +X+Y-Z octant
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitHighZXNeg");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitHighZXNeg");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[7], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j < 0 && _x2j < 0 && _x3j > 0)
                    {
                        // node j in -X-Y+Z octant
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitHighZXPos");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForce);
                    }
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else if (nb == 11) // horizon sphere also intersects with V1, V2, V4
                {
                    if (_x1j > 0 && _x2j > 0 && _x3j > 0)
                    {
                        // node j in +X+Y+Z octant
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitLowYXPos");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[7], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j < 0 && _x2j < 0 && _x3j < 0)
                    {
                        // node j in -X-Y-Z octant
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitLowYXNeg");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[7], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowYXNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else
                {
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
				break;
            }
        case 9: // b=0000|1001, horizon sphere intersects with T1, T4 only
            {
                // check the intersection points with the other four edges
                double nPoints[4]; // store the coordinate of intersectinon points
                int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                if (nb == 8) // horizon sphere also intersects with V4
                {
                    if (_x1j < 0 && _x2j > 0 && _x3j > 0)
                    {
                        // node j in -X+Y+Z octant
                        flag.push_back("YLimitHighZXPos");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitHighZXPos");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[3], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j > 0 && _x2j < 0 && _x3j < 0)
                    {
                        // node j in +X-Y-Z octant
                        flag.push_back("YLimitHighZXNeg");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForce);
                    }
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else if (nb == 7) // horizon sphere also intersects with V1, V2, V3
                {
                    if (_x1j < 0 && _x2j > 0 && _x3j < 0)
                    {
                        // node j in -X+Y-Z octant
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitHighYXPos");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[3], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j > 0 && _x2j < 0 && _x3j > 0)
                    {
                        // node j in +X-Y+Z octant
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitHighYXNeg");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[3], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitHighYXNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else
                {
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
				break;
            }
        case 144: // b=1001|0000, horizon sphere intersects with B1, B4 only
            {
                // check the intersection points with the other four edges
                double nPoints[4]; // store the coordinate of intersectinon points
                int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                if (nb == 8) // horizon sphere also intersects with V4
                {
                    if (_x1j < 0 && _x2j > 0 && _x3j < 0)
                    {
                        // node j in -X+Y-Z octant
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitHighZXNeg");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[7], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitHighZXNeg");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[7], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j > 0 && _x2j < 0 && _x3j > 0)
                    {
                        // node j in +X-Y+Z octant
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitHighZXPos");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[7], flag, pdForce);
                    }
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else if (nb == 7) // horizon sphere also intersects with V1, V2, V3
                {
                    if (_x1j < 0 && _x2j > 0 && _x3j > 0)
                    {
                        // node j in -X+Y+Z octant
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitLowYXPos");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[7], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j > 0 && _x2j < 0 && _x3j < 0)
                    {
                        // node j in +X-Y-Z octant
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitLowYXNeg");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[7], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowYXNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[7], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else
                {
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
				break;
            }
        case 6: // b=0000|0110, horizon sphere intersects with T2, T3 only
            {
                // check the intersection points with the other four edges
                double nPoints[4]; // store the coordinate of intersectinon points
                int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                if (nb == 2) // horizon sphere also intersects with V2
                {
                    if (_x1j > 0 && _x2j < 0 && _x3j > 0)
                    {
                        // node j in +X-Y+Z octant
                        flag.push_back("YLimitLowZXPos");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[2], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitLowZXPos");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[2], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j < 0 && _x2j > 0 && _x3j < 0)
                    {
                        // node j in -X+Y-Z octant
                        flag.push_back("YLimitLowZXNeg");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[2], _x1j + 0.5 * _dx, flag, pdForce);
                    }
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else if (nb == 13) // horizon sphere also intersects with V1, V3, V4
                {
                    if (_x1j > 0 && _x2j < 0 && _x3j < 0)
                    {
                        // node j in +X-Y-Z octant
                        flag.push_back("YLimitHighYXNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[2], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j < 0 && _x2j > 0 && _x3j > 0)
                    {
                        // node j in -X+Y+Z octant
                        flag.push_back("YLimitHighYXPos");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[2], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitHighYXPos");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[2], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else
                {
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
				break;
            }
        case 5: // b=0000|0101,  horizon sphere intersects with T1, T2 only
            {
                // check the intersection points with the other four edges
                double nPoints[4]; // store the coordinate of intersectinon points
                int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                if (nb == 1) // horizon sphere also intersects with V1
                {
                    if (_x1j < 0 && _x2j < 0 && _x3j > 0)
                    {
                        // node j in -X-Y+Z octant
                        flag.push_back("YLimitLowZXPos");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitLowZXPos");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[2], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j > 0 && _x2j > 0 && _x3j < 0)
                    {
                        // node j in +X+Y-Z octant
                        flag.push_back("YLimitLowZXNeg");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForce);
                    }
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else if (nb == 14) // horizon sphere also intersects with V2, V3, V4
                {
                    if (_x1j < 0 && _x2j < 0 && _x3j < 0)
                    {
                        // node j in -X-Y-Z octant
                        flag.push_back("YLimitHighYXNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[2], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j > 0 && _x2j > 0 && _x3j > 0)
                    {
                        // node j in +X+Y+Z octant
                        flag.push_back("YLimitHighYXPos");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[2], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitHighYXPos");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else
                {
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
				break;
            }
        case 96: // b=0110|0000,  horizon sphere intersects with B2, B3 only
            {
                // check the intersection points with the other four edges
                double nPoints[4]; // store the coordinate of intersectinon points
                int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                if (nb == 2) // horizon sphere also intersects with V2
                {
                    if (_x1j > 0 && _x2j < 0 && _x3j < 0)
                    {
                        // node j in +X-Y-Z octant
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitLowZXNeg");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowZXNeg");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[6], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j < 0 && _x2j > 0 && _x3j > 0)
                    {
                        // node j in -X+Y+Z octant
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitLowZXPos");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForce);
                    }
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else if (nb == 13) // horizon sphere also intersects with V1, V3, V4
                {
                    if (_x1j > 0 && _x2j < 0 && _x3j > 0)
                    {
                        // node j in +X-Y+Z octant
                        flag.push_back("YLimitLowYXNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[6], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j < 0 && _x2j > 0 && _x3j < 0)
                    {
                        // node j in -X+Y-Z octant
                        flag.push_back("YLimitLowYXPos");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[6], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitLowYXPos");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else
                {
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
				break;
            }
        case 80: // b=0101|0000, horizon sphere intersects with B1, B2 only
            {
                // check the intersection points with the other four edges
                double nPoints[4]; // store the coordinate of intersectinon points
                int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                if (nb == 1) // horizon sphere also intersects with V1
                {
                    if (_x1j < 0 && _x2j < 0 && _x3j < 0)
                    {
                        // node j in -X-Y-Z octant
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitLowZXNeg");
                        flag.push_back("ZLimitYSphereNeg");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[6], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowZXNeg");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[6], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitZHigh");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYHigh");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j > 0 && _x2j > 0 && _x3j > 0)
                    {
                        // node j in +X+Y+Z octant
                        flag.push_back("YLimitZLow");
                        flag.push_back("YLimitLowZXPos");
                        flag.push_back("ZLimitYLow");
                        flag.push_back("ZLimitYSpherePos");
                        flag.push_back("XZY");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[6], flag, pdForce);
                    }
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else if (nb == 14) // horizon sphere also intersects with V2, V3, V4
                {
                    if (_x1j < 0 && _x2j < 0 && _x3j > 0)
                    {
                        // node j in -X-Y+Z octant
                        flag.push_back("YLimitLowYXNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[6], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else if (_x1j > 0 && _x2j > 0 && _x3j < 0)
                    {
                        // node j in +X+Y-Z octant
                        flag.push_back("YLimitLowYXPos");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[6], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitLowYXPos");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[6], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
                else
                {
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
				break;
            }
        case 53: // b=0011|0101, horizon sphere intersects with T1, T2, B1, B3
            {
                if (_x1j > 0 && _x2j > 0 && _x3j > 0)
                {
                    // node j in +X+Y+Z octant, horizon sphere also intersects with V2
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitLowZXPos");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYSpherePos");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[2], _x1j + 0.5 * _dx, flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYSpherePos");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
                    ArrayAddition(pdForceStep, pdForce, 3);
                }
                else if (_x1j < 0 && _x2j < 0 && _x3j < 0)
                {
                    // node j in -X-Y-Z octant, horizon sphere also intersects with V2
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitLowZXNeg");
                    flag.push_back("ZLimitYSphereNeg");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[2], _x1j + 0.5 * _dx, flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYSphereNeg");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitLowZXNeg");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[2], _x1j + 0.5 * _dx, flag, pdForceStep);
                    ArrayAddition(pdForceStep, pdForce, 3);
                }
                else
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
				break;
            }
        case 83: // b=0101|0011, horizon sphere intersects with T1, T3, B1, B2
            {
                if (_x1j > 0 && _x2j > 0 && _x3j < 0)
                {
                    // node j in +X+Y-Z octant, horizon sphere also intersects with V2
                    flag.push_back("YLimitLowZXNeg");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYSpherePos");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYSpherePos");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[6], flag, pdForceStep);
                    ArrayAddition(pdForceStep, pdForce, 3);
                }
                else if (_x1j < 0 && _x2j < 0 && _x3j > 0)
                {
                    // node j in -X-Y+Z octant, horizon sphere also intersects with V2
                    flag.push_back("YLimitLowZXPos");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYSphereNeg");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYSphereNeg");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[6], flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitLowZXPos");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
                    ArrayAddition(pdForceStep, pdForce, 3);
                }
                else
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
				break;
            }
        case 54: // b=0011|0110, horizon sphere intersects with T2, T3, B1, B3
            {
                if (_x1j < 0 && _x2j > 0 && _x3j > 0)
                {
                    // node j in -X+Y+Z octant, horizon sphere also intersects with V1
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitLowZXPos");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYSpherePos");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYSpherePos");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[2], _x1j + 0.5 * _dx, flag, pdForceStep);
                    ArrayAddition(pdForceStep, pdForce, 3);
                }
                else if (_x1j > 0 && _x2j < 0 && _x3j < 0)
                {
                    // node j in +X-Y-Z octant, horizon sphere also intersects with V1
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitLowZXNeg");
                    flag.push_back("ZLimitYSphereNeg");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYSphereNeg");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[2], _x1j + 0.5 * _dx, flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitLowZXNeg");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[2], flag, pdForceStep);
                    ArrayAddition(pdForceStep, pdForce, 3);
                }
                else
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
				break;
            }
        case 99: // b=0110|0011, horizon sphere intersects with T1, T3, B2, B3
            {
                if (_x1j < 0 && _x2j > 0 && _x3j < 0)
                {
                    // node j in -X+Y-Z octant, horizon sphere also intersects with V1
                    flag.push_back("YLimitLowZXNeg");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYSpherePos");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[6], flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYSpherePos");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
                    ArrayAddition(pdForceStep, pdForce, 3);
                }
                else if (_x1j > 0 && _x2j < 0 && _x3j > 0)
                {
                    // node j in +X-Y+Z octant, horizon sphere also intersects with V1
                    flag.push_back("YLimitLowZXPos");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYSphereNeg");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[6], flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYSphereNeg");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[6], _x1j + 0.5 * _dx, flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitLowZXPos");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[6], flag, pdForceStep);
                    ArrayAddition(pdForceStep, pdForce, 3);
                }
                else
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
				break;
            }
        case 57: // b=0011|1001, horizon sphere intersects with T1, T4, B1, B3
            {
                if (_x1j > 0 && _x2j < 0 && _x3j > 0)
                {
                    // node j in +X-Y+Z octant, horizon sphere also intersects with V3
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitHighZXPos");
                    flag.push_back("ZLimitYSphereNeg");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[3], _x1j + 0.5 * _dx, flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYSphereNeg");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
                    ArrayAddition(pdForceStep, pdForce, 3);
                }
                else if (_x1j < 0 && _x2j > 0 && _x3j < 0)
                {
                    // node j in -X+Y-Z octant, horizon sphere also intersects with V3
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitHighZXNeg");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYSpherePos");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[3], _x1j + 0.5 * _dx, flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYSpherePos");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitHighZXNeg");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[3], _x1j + 0.5 * _dx, flag, pdForceStep);
                    ArrayAddition(pdForceStep, pdForce, 3);
                }
                else
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
				break;
            }
        case 147: // b=1001|0011, horizon sphere intersects with T1, T3, B1, B4
            {
                if (_x1j > 0 && _x2j < 0 && _x3j < 0)
                {
                    // node j in +X-Y-Z octant, horizon sphere also intersects with V3
                    flag.push_back("YLimitHighZXNeg");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYSphereNeg");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYSphereNeg");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[7], flag, pdForceStep);
                    ArrayAddition(pdForceStep, pdForce, 3);
                }
                else if (_x1j < 0 && _x2j > 0 && _x3j > 0)
                {
                    // node j in -X+Y+Z octant, horizon sphere also intersects with V3
                    flag.push_back("YLimitHighZXPos");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYSpherePos");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYSpherePos");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[7], flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitHighZXPos");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
                    ArrayAddition(pdForceStep, pdForce, 3);
                }
                else
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
				break;
            }
        case 58: // b=0011|1010, horizon sphere intersects with T3, T4, B1, B3
            {
                if (_x1j < 0 && _x2j < 0 && _x3j > 0)
                {
                    // node j in -X-Y+Z octant, horizon sphere also intersects with V4
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitHighZXPos");
                    flag.push_back("ZLimitYSphereNeg");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYSphereNeg");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[3], _x1j + 0.5 * _dx, flag, pdForceStep);
                    ArrayAddition(pdForceStep, pdForce, 3);
                }
                else if (_x1j > 0 && _x2j > 0 && _x3j < 0)
                {
                    // node j in +X+Y-Z octant, horizon sphere also intersects with V4
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitHighZXNeg");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYSpherePos");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYSpherePos");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[3], _x1j + 0.5 * _dx, flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitHighZXNeg");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[3], flag, pdForceStep);
                    ArrayAddition(pdForceStep, pdForce, 3);
                }
                else
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
				break;
            }
        case 163: // b=1010|0011, horizon sphere intersects with T1, T3, B4, B3
            {
                if (_x1j < 0 && _x2j < 0 && _x3j < 0)
                {
                    // node j in -X-Y-Z octant, horizon sphere also intersects with V4
                    flag.push_back("YLimitHighZXNeg");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYSphereNeg");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[7], flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYSphereNeg");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
                    ArrayAddition(pdForceStep, pdForce, 3);
                }
                else if (_x1j > 0 && _x2j > 0 && _x3j > 0)
                {
                    // node j in +X+Y+Z octant, horizon sphere also intersects with V4
                    flag.push_back("YLimitHighZXPos");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYSpherePos");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[7], flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYSpherePos");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, iPoints[7], _x1j + 0.5 * _dx, flag, pdForceStep);
					ArrayAddition(pdForceStep, pdForce, 3);
                    flag.clear();
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitHighZXPos");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, iPoints[7], flag, pdForceStep);
                    ArrayAddition(pdForceStep, pdForce, 3);
                }
                else
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
				break;
            }
        case 204: // b=1100|1100, horizon sphere intersects with T2, T4, B2, B4
            {
                if (_x1j > 0)
                {
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitXLow");
                    flag.push_back("ZLimitXSpherePos");
                    flag.push_back("YZX");
                    CalcIntg3D(itg, _x2j - 0.5 * _dx, _x2j + 0.5 * _dx, flag, pdForce);
                }
                else
                {
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitXSphereNeg");
                    flag.push_back("ZLimitXHigh");
                    flag.push_back("YZX");
                    CalcIntg3D(itg, _x2j - 0.5 * _dx, _x2j + 0.5 * _dx, flag, pdForce);
                }
				break;
            }
        case 51: // b=0011|0011, horizon sphere intersects with T1, T3, B1, B3
            {
                if (_x2j > 0)
                {
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYLow");
                    flag.push_back("ZLimitYSpherePos");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, _x1j + 0.5 * _dx, flag, pdForce);
                }
                else
                {
                    flag.push_back("YLimitZLow");
                    flag.push_back("YLimitZHigh");
                    flag.push_back("ZLimitYSphereNeg");
                    flag.push_back("ZLimitYHigh");
                    flag.push_back("XZY");
                    CalcIntg3D(itg, _x1j - 0.5 * _dx, _x1j + 0.5 * _dx, flag, pdForce);
                }
				break;
            }
        case 0: // b=0000|0000, horizon sphere does not intersects with any eight edges
            {
                // check the intersection points with the other four edges
                double nPoints[4]; // store the coordinate of intersectinon points
                int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                if (nb == 15)
                {
                    // horizon sphere intersects with V1-V4
                    if (_x3j > 0)
                    {
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZLow");
                        flag.push_back("ZLimitZSpherePos");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, _x1j + 0.5 * _dx, flag, pdForce);
                    }
                    else
                    {
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitZSphereNeg");
                        flag.push_back("ZLimitZHigh");
                        flag.push_back("XYZ");
                        CalcIntg3D(itg, _x1j - 0.5 * _dx, _x1j + 0.5 * _dx, flag, pdForce);
                    }
                }
                else
                {
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
				break;
            }
        case 197: // b=1100|0101, horizon sphere intersects with T1, T2, B2, B4
            {
                if (_x1j > 0 && _x2j > 0 && _x3j > 0)
                {
                    // node j in +X+Y+Z octant
                    // check the intersection points with the other four edges
                    double nPoints[4]; // store the coordinate of intersectinon points
                    int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                    if (nb == 8)
                    {
                        // horizon sphere also intersects with V4
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitLowYZPos");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXSpherePos");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[3], _x3j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXSpherePos");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[3], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else if (_x1j < 0 && _x2j < 0 && _x3j < 0)
                {
                    // node j in -X-Y-Z octant
                    // check the intersection points with the other four edges
                    double nPoints[4]; // store the coordinate of intersectinon points
                    int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                    if (nb == 8)
                    {
                        // horizon sphere also intersects with V4
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitLowYZNeg");
                        flag.push_back("ZLimitXSphereNeg");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[3], _x3j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXSphereNeg");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[3], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowYZNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[3], _x3j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else
                {
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
				break;
            }
        case 92: // b=0101|1100, horizon sphere intersects with T2, T4, B1, B2
            {
                if (_x1j > 0 && _x2j > 0 && _x3j < 0)
                {
                    // node j in +X+Y-Z octant
                    // check the intersection points with the other four edges
                    double nPoints[4]; // store the coordinate of intersectinon points
                    int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                    if (nb == 8)
                    {
                        // horizon sphere also intersects with V4
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitLowYZPos");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXSpherePos");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[3], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXSpherePos");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[3], _x3j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else if (_x1j < 0 && _x2j < 0 && _x3j > 0)
                {
                    // node j in -X-Y+Z octant
                    // check the intersection points with the other four edges
                    double nPoints[4]; // store the coordinate of intersectinon points
                    int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                    if (nb == 8)
                    {
                        // horizon sphere also intersects with V4
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitLowYZNeg");
                        flag.push_back("ZLimitXSphereNeg");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[3], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXSphereNeg");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[3], _x3j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitLowYZNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[3], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else
                {
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
				break;
            }
        case 201: // b=1100|1001, horizon sphere intersects with T1, T4, B2, B4
            {
                if (_x1j > 0 && _x2j < 0 && _x3j > 0)
                {
                    // node j in +X-Y+Z octant
                    // check the intersection points with the other four edges
                    double nPoints[4]; // store the coordinate of intersectinon points
                    int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                    if (nb == 1)
                    {
                        // horizon sphere also intersects with V1
                        flag.push_back("YLimitLowYZNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXSpherePos");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[0], _x3j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXSpherePos");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[0], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else if (_x1j < 0 && _x2j > 0 && _x3j < 0)
                {
                    // node j in -X+Y-Z octant
                    // check the intersection points with the other four edges
                    double nPoints[4]; // store the coordinate of intersectinon points
                    int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                    if (nb == 1)
                    {
                        // horizon sphere also intersects with V1
                        flag.push_back("YLimitLowYZPos");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXSphereNeg");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[0], _x3j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXSphereNeg");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[0], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitLowYZPos");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[0], _x3j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else
                {
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
				break;
            }
        case 156: // b=1001|1100, horizon sphere intersects with T2, T4, B1, B4
            {
                if (_x1j > 0 && _x2j < 0 && _x3j < 0)
                {
                    // node j in +X-Y-Z octant
                    // check the intersection points with the other four edges
                    double nPoints[4]; // store the coordinate of intersectinon points
                    int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                    if (nb == 1)
                    {
                        // horizon sphere also intersects with V1
                        flag.push_back("YLimitLowYZNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXSpherePos");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[0], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXSpherePos");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[0], _x3j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else if (_x1j < 0 && _x2j > 0 && _x3j > 0)
                {
                    // node j in -X+Y+Z octant
                    // check the intersection points with the other four edges
                    double nPoints[4]; // store the coordinate of intersectinon points
                    int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                    if (nb == 1)
                    {
                        // horizon sphere also intersects with V1
                        flag.push_back("YLimitLowYZPos");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXSphereNeg");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[0], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXSphereNeg");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[0], _x3j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitLowYZPos");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[0], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else
                {
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
				break;
            }
        case 198: // b=1100|0110, horizon sphere intersects with T2, T3, B2, B4
            {
                if (_x1j < 0 && _x2j > 0 && _x3j > 0)
                {
                    // node j in -X+Y+Z octant
                    // check the intersection points with the other four edges
                    double nPoints[4]; // store the coordinate of intersectinon points
                    int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                    if (nb == 4)
                    {
                        // horizon sphere also intersects with V3
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitHighYZPos");
                        flag.push_back("ZLimitXSphereNeg");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[2], _x3j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXSphereNeg");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[2], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else if (_x1j > 0 && _x2j < 0 && _x3j < 0)
                {
                    // node j in +X-Y-Z octant
                    // check the intersection points with the other four edges
                    double nPoints[4]; // store the coordinate of intersectinon points
                    int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                    if (nb == 4)
                    {
                        // horizon sphere also intersects with V3
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitHighYZNeg");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXSpherePos");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[2], _x3j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXSpherePos");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[2], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitHighYZNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[2], _x3j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else
                {
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
				break;
            }
        case 108: // b=0110|1100, horizon sphere intersects with T2, T4, B2, B3
            {
                if (_x1j < 0 && _x2j > 0 && _x3j < 0)
                {
                    // node j in -X+Y-Z octant
                    // check the intersection points with the other four edges
                    double nPoints[4]; // store the coordinate of intersectinon points
                    int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                    if (nb == 4)
                    {
                        // horizon sphere also intersects with V3
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitHighYZPos");
                        flag.push_back("ZLimitXSphereNeg");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[2], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXSphereNeg");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[2], _x3j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else if (_x1j > 0 && _x2j < 0 && _x3j > 0)
                {
                    // node j in +X-Y+Z octant
                    // check the intersection points with the other four edges
                    double nPoints[4]; // store the coordinate of intersectinon points
                    int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                    if (nb == 4)
                    {
                        // horizon sphere also intersects with V3
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitHighYZNeg");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXSpherePos");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[2], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXSpherePos");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[2], _x3j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitHighYZNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[2], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else
                {
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
				break;
            }
        case 202: // b=1100|1010, horizon sphere intersects with T3, T4, B2, B4
            {
                if (_x1j < 0 && _x2j < 0 && _x3j > 0)
                {
                    // node j in -X-Y+Z octant
                    // check the intersection points with the other four edges
                    double nPoints[4]; // store the coordinate of intersectinon points
                    int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                    if (nb == 2)
                    {
                        // horizon sphere also intersects with V2
                        flag.push_back("YLimitHighYZNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXSphereNeg");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[1], _x3j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXSphereNeg");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[1], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else if (_x1j > 0 && _x2j > 0 && _x3j < 0)
                {
                    // node j in +X+Y-Z octant
                    // check the intersection points with the other four edges
                    double nPoints[4]; // store the coordinate of intersectinon points
                    int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                    if (nb == 2)
                    {
                        // horizon sphere also intersects with V2
                        flag.push_back("YLimitHighYZPos");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXSpherePos");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[1], _x3j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXSpherePos");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[1], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitHighYZPos");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[1], _x3j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else
                {
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
				break;
            }
        case 172: // b=1010|1100, horizon sphere intersects with T2, T4, B3, B4
            {
                if (_x1j < 0 && _x2j < 0 && _x3j < 0)
                {
                    // node j in -X-Y-Z octant
                    // check the intersection points with the other four edges
                    double nPoints[4]; // store the coordinate of intersectinon points
                    int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                    if (nb == 2)
                    {
                        // horizon sphere also intersects with V2
                        flag.push_back("YLimitHighYZNeg");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXSphereNeg");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[1], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXSphereNeg");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[1], _x3j + 0.5 * _dx, flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else if (_x1j > 0 && _x2j > 0 && _x3j > 0)
                {
                    // node j in +X+Y+Z octant
                    // check the intersection points with the other four edges
                    double nPoints[4]; // store the coordinate of intersectinon points
                    int nb = this->FindIntersectionPointOnFourEdges(funcName, nPoints);
                    if (nb == 2)
                    {
                        // horizon sphere also intersects with V2
                        flag.push_back("YLimitHighYZPos");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXSpherePos");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[1], flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitYHigh");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXSpherePos");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, nPoints[1], _x3j + 0.5 * _dx, flag, pdForceStep);
						ArrayAddition(pdForceStep, pdForce, 3);
                        flag.clear();
                        flag.push_back("YLimitYLow");
                        flag.push_back("YLimitHighYZPos");
                        flag.push_back("ZLimitXLow");
                        flag.push_back("ZLimitXHigh");
                        flag.push_back("ZYX");
                        CalcIntg3D(itg, _x3j - 0.5 * _dx, nPoints[1], flag, pdForceStep);
                        ArrayAddition(pdForceStep, pdForce, 3);
                    }
                    else
                    {
                        pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                    }
                }
                else
                {
                    pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
                }
				break;
            }
        default: // any other unfound possibilities
            {
                pdForce[0] = pdForce[1] = pdForce[2] = 0.0;
				break;
            }
    }
}

void pdIntgManager::CalcConfigurationType3(string itg, double* pdForce)
{
    // this function calculates integration for geometric configuration when the family node cell is fully in the horizon sphere
    vector<string> flag; // table index for limit function
    // integrating over the whole cell
    flag.push_back("YLimitYLow");
    flag.push_back("YLimitYHigh");
    flag.push_back("ZLimitZLow");
    flag.push_back("ZLimitZHigh");
    flag.push_back("XYZ");
    return CalcIntg3D(itg, _x1j - 0.5 * _dx, _x1j + 0.5 * _dx, flag, pdForce);
}









