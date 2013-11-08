#include "pdDataManager.h"
using namespace std;

void pdDataManager::ReadInputFile()
{
	// This function read in all data from the input file and initial all containers.

	ofstream outfile(_outFilePath.c_str(), ios::app);
	string tempKeyword;

	/* 1. basic problem setups 
	*/
	// read problem title (optional)
	tempKeyword = "problem_title";
    if (FindKeyword(tempKeyword, OPT))
	{
        ReadKeyword(tempKeyword, 1, _title);
	}
    else
	{
        _title = "Peridynamic test problem";
	}
	// read grid center (optional)
	tempKeyword = "grid_center";
    if (FindKeyword(tempKeyword, OPT))
	{
        ReadKeyword(tempKeyword, 1, _gridCenterX1);
		ReadKeyword(tempKeyword, 2, _gridCenterX2);
		ReadKeyword(tempKeyword, 3, _gridCenterX3);
	}
    else
	{
        _gridCenterX1 = 0.0;
		_gridCenterX2 = 0.0;
		_gridCenterX3 = 0.0;
	}       
	// read grid dimension (required)
	tempKeyword = "grid_dimension";
	if (FindKeyword(tempKeyword, REQ))
	{
		ReadKeyword(tempKeyword, 1, _gridDimensionX1);
		ReadKeyword(tempKeyword, 2, _gridDimensionX2);
		ReadKeyword(tempKeyword, 3, _gridDimensionX3);
	}
	// read global grid spacing (required)
	tempKeyword = "grid_spacing";
	if (FindKeyword(tempKeyword, REQ))
	{
		ReadKeyword(tempKeyword, 1, _gridSpacing);
	}
	// read horizon (optional)
	tempKeyword = "horizon";
	if (FindKeyword(tempKeyword, OPT))
	{
		ReadKeyword("horizon", 1, _horizon);
	}
	else
	{
		_horizon = 3.015 * _gridSpacing; // default
	}
	// read total run time and total run steps (required)
	tempKeyword = "total_time";
	if (FindKeyword(tempKeyword, REQ))
	{
		ReadKeyword(tempKeyword, 1, _totalTime);
	}
	tempKeyword = "total_step";
	if (FindKeyword(tempKeyword, REQ))
	{
		ReadKeyword(tempKeyword, 1, _totalStep);
	}
	// read plot dump frequency (required)
	tempKeyword = "step_dump_frequency";
	if (FindKeyword(tempKeyword, REQ))
	{
		ReadKeyword(tempKeyword, 1, _stepDumpFrequency);
	}
	// read safety factor (optional)
	tempKeyword = "safety_factor";
	if (FindKeyword(tempKeyword, OPT))
	{
		ReadKeyword(tempKeyword, 1, _safeFactor);
	}
	else
	{
		_safeFactor = 0.8; // default
	}

	if (_myRank==0)
	{
		outfile << "ReadInputFile: Problem title = " << _title << endl
			<< "ReadInputFile: Grid center: x1 = " << _gridCenterX1 << ", x2 = " << _gridCenterX2 << ", x3 = " << _gridCenterX3 << endl
			<< "ReadInputFile: Grid dimension = " << _gridDimensionX1 << " * " << _gridDimensionX2 << " * " << _gridDimensionX3 << endl
			<< "ReadInputFile: Global grid spacing = " << _gridSpacing << endl
			<< "ReadInputFile: Global grid horizon = " << _horizon << endl
			<< "ReadInputFile: Total run time = " << _totalTime << endl
			<< "ReadInputFile: Total run steps = " << _totalStep << endl
			<< "ReadInputFile: Step dump frequency = " << _stepDumpFrequency << endl;
	}
	// read fixed step length (optional), using stable time step as default
	_fixedTimeStepOn = false;// default
	tempKeyword = "fixed_time_step";
	if (FindKeyword(tempKeyword, OPT))
	{
		_fixedTimeStepOn = true;
		ReadKeyword(tempKeyword, 1, _fixedTimeStep);
		if (_myRank==0)
		{
			outfile << "ReadInputFile: Fixed time step = " << _fixedTimeStep << endl;
		}
	}

	/* 2. material
	*/
	// read material data (required), id starts from 0
	// also check if state-based material is used
	int stateSum = 0;
	_stateMatOn = false;
	tempKeyword = "number_of_materials";
	if (FindKeyword(tempKeyword, REQ))
	{
		ReadKeyword(tempKeyword, 1, _numMats);
		string keywd = "material_";
		for (int i=0;i!=_numMats;++i)
		{
			keywd = KeywordAppend(keywd, i+1);
			// get material type
			int mat_type;
			ReadKeyword(keywd, 1, mat_type);		
			switch(mat_type)
			{
			case BOND: // type 1 for bond-based material
				{
					double ymod, sigmay, ecrit, denst;
					ReadKeyword(keywd, 2, ymod);
					ReadKeyword(keywd, 3, sigmay);
					ReadKeyword(keywd, 4, ecrit);
					ReadKeyword(keywd, 5, denst);
					_mats.push_back(pdHandle<pdMaterial>(new pdMaterialBond(i, mat_type, ymod, sigmay, ecrit, denst)));
					break;
				}
			case STATE: // type 2 for state-based material
				{
					double ymod, sigmay, nu, ecrit, denst;
					ReadKeyword(keywd, 2, ymod);
					ReadKeyword(keywd, 3, sigmay);
					ReadKeyword(keywd, 4, nu);
					ReadKeyword(keywd, 5, ecrit);
					ReadKeyword(keywd, 6, denst);
					_mats.push_back(pdHandle<pdMaterial>(new pdMaterialState(i, mat_type, ymod, sigmay, nu, ecrit, denst)));
					stateSum++;
					break;
				}
			default:
				{
				}
			}
			// read body force (gavity) data (optional)
			string gravString = KeywordAppend("gravity_", i+1);
			if (FindKeyword(gravString, OPT))
			{
				double grav1, grav2, grav3;
				ReadKeyword(gravString, 1, grav1);
				ReadKeyword(gravString, 2, grav2);
				ReadKeyword(gravString, 3, grav3);
				_mats[i]->SetGravity(grav1, grav2, grav3);
			}
			else
			{
				_mats[i]->SetGravity(0.0, 0.0, 0.0); // default
			}
			keywd = "material_";
		}
		if (stateSum>0)
		{
			_stateMatOn = true;
		}
	}

	/* 3. boundary condition 
	*/
	// read boundary condition data (required), id starts from 0
	tempKeyword = "number_of_boundary_conditions";
	if (FindKeyword(tempKeyword, REQ))
	{
		ReadKeyword(tempKeyword, 1, _numBCs);
		string keywd = "boundary_condition_";
		for (int i=0;i!=_numBCs;++i)
		{
			keywd = KeywordAppend(keywd, i+1);
			// get boundary condition type
			int bc_type;
			ReadKeyword(keywd, 1, bc_type);
			switch (bc_type)
			{		
				case DISPL:	// type 1 for prescribed displacement boundary condition
				{
					int dir_x1, dir_x2, dir_x3;
					ReadKeyword(keywd, 2, dir_x1);
					ReadKeyword(keywd, 3, dir_x2);
					ReadKeyword(keywd, 4, dir_x3);
					double val[3];
					ReadKeyword(keywd, 5, val[0]);
					ReadKeyword(keywd, 6, val[1]);
					ReadKeyword(keywd, 7, val[2]);
					vector<double> value(val, val+sizeof(val)/sizeof(double));
					double end_time;
					ReadKeyword(keywd, 8, end_time);
					double load_type;
					ReadKeyword(keywd, 9, load_type);
					_BCs.push_back(pdHandle<pdBdryCondition>(new pdBdryConditionDispl(i, dir_x1, dir_x2, dir_x3, value, end_time, load_type)));
					break;
				}		
				case VELOC:	// type 2 for prescibed veloctiy boundary condition
				{
					int dir_x1, dir_x2, dir_x3;
					ReadKeyword(keywd, 2, dir_x1);
					ReadKeyword(keywd, 3, dir_x2);
					ReadKeyword(keywd, 4, dir_x3);
					double val[3];
					ReadKeyword(keywd, 5, val[0]);
					ReadKeyword(keywd, 6, val[1]);
					ReadKeyword(keywd, 7, val[2]);
					vector<double> value(val, val+sizeof(val)/sizeof(double));
					double end_time;
					ReadKeyword(keywd, 8, end_time);
					double load_type;
					ReadKeyword(keywd, 9, load_type);
					_BCs.push_back(pdHandle<pdBdryCondition>(new pdBdryConditionVeloc(i, dir_x1, dir_x2, dir_x3, value, end_time, load_type)));
					break;
				}		
				case DISPLGRAD:	// type 3 for prescribed displacement gradient boundary condition (12/06/09)
				{
					int dir_x1, dir_x2, dir_x3;
					ReadKeyword(keywd, 2, dir_x1);
					ReadKeyword(keywd, 3, dir_x2);
					ReadKeyword(keywd, 4, dir_x3);
					double val[9];
					ReadKeyword(keywd, 5, val[0]);
					ReadKeyword(keywd, 6, val[1]);
					ReadKeyword(keywd, 7, val[2]);
					ReadKeyword(keywd, 8, val[3]);
					ReadKeyword(keywd, 9, val[4]);
					ReadKeyword(keywd, 10, val[5]);
					ReadKeyword(keywd, 11, val[6]);
					ReadKeyword(keywd, 12, val[7]);
					ReadKeyword(keywd, 13, val[8]);
					vector<double> value(val, val+sizeof(val)/sizeof(double));
					double end_time;
					ReadKeyword(keywd, 14, end_time);
					double load_type;
					ReadKeyword(keywd, 15, load_type);
					_BCs.push_back(pdHandle<pdBdryCondition>(new pdBdryConditionDisplGrad(i, dir_x1, dir_x2, dir_x3, value, end_time, load_type)));
					break;
				}
				case TRACTION:	// type 4 for traction force boundary condition
				{
					int dir_x1, dir_x2, dir_x3;
					ReadKeyword(keywd, 2, dir_x1);
					ReadKeyword(keywd, 3, dir_x2);
					ReadKeyword(keywd, 4, dir_x3);
					double val[3];
					ReadKeyword(keywd, 5, val[0]);
					ReadKeyword(keywd, 6, val[1]);
					ReadKeyword(keywd, 7, val[2]);
					vector<double> value(val, val+sizeof(val)/sizeof(double));
					double end_time;
					ReadKeyword(keywd, 8, end_time);
					double load_type;
					ReadKeyword(keywd, 9, load_type);
					_BCs.push_back(pdHandle<pdBdryCondition>(new pdBdryConditionTract(i, dir_x1, dir_x2, dir_x3, value, end_time, load_type)));
					break;
				}		
				default:
				{
				}
			}
			keywd = "boundary_condition_";
		}
		// print boundary condition data
		if (_myRank==0)
		{
			outfile << endl << "ReadInputFile: Total boundary condition number = " << _BCs.size() << endl;
			for (vector<pdHandle<pdBdryCondition>>::iterator itor=_BCs.begin();itor!=_BCs.end();++itor)
			{
				outfile << *(*itor).GetPtr() << endl;
			}
		}
	}

	/* 4. material spaces
	*/
	// read material spaces data (required), id starts from 0
	tempKeyword = "number_of_material_spaces";
	if (FindKeyword(tempKeyword, REQ))
	{
		ReadKeyword(tempKeyword, 1, _numMatsSpaces);
		string keywd = "material_space_";
		for (int i=0;i!=_numMatsSpaces;++i)
		{
			keywd = KeywordAppend(keywd, i+1);
			// get space geometry type
			int geo_type, mid;
			ReadKeyword(keywd, 1, geo_type);
			switch(geo_type)
			{
			case CUBOID: // type 1 for cuboidal space
				{
					// get id of the material associated with this material space
					ReadKeyword(keywd, 2, mid);
					mid = mid - 1; // material id in code starts from zero
					// get space boundary
					double reg_x1lo, reg_x1hi, reg_x2lo, reg_x2hi, reg_x3lo, reg_x3hi;
					ReadKeyword(keywd, 3, reg_x1lo);
					ReadKeyword(keywd, 4, reg_x1hi);
					ReadKeyword(keywd, 5, reg_x2lo);
					ReadKeyword(keywd, 6, reg_x2hi);
					ReadKeyword(keywd, 7, reg_x3lo);
					ReadKeyword(keywd, 8, reg_x3hi);
					pdSpaceMaterial* matspace = new pdSpaceMaterial(i, geo_type, GetMaterial(mid));
					matspace->SetCuboidSpaceBdry(reg_x1lo, reg_x1hi, reg_x2lo, reg_x2hi, reg_x3lo, reg_x3hi);
					_matSpaces.push_back(matspace);
					break;
				}
			case CYLINDER: // type 2 for cylindrical space
				{
					// get material id for this space
					ReadKeyword(keywd, 2, mid);
					mid = mid - 1;
					// get space boundary
					double reg_rad, reg_x1cen, reg_x2cen, reg_x3lo, reg_x3hi;
					ReadKeyword(keywd, 3, reg_rad);
					ReadKeyword(keywd, 4, reg_x1cen);
					ReadKeyword(keywd, 5, reg_x2cen);
					ReadKeyword(keywd, 6, reg_x3lo);
					ReadKeyword(keywd, 7, reg_x3hi);
					pdSpaceMaterial* matspace = new pdSpaceMaterial(i, geo_type, GetMaterial(mid));
					matspace->SetCylinderSpaceBdry(reg_rad, reg_x1cen, reg_x2cen, reg_x3lo, reg_x3hi);
					_matSpaces.push_back(matspace);
					break;
				}
			default:
				{
				}
			}
			keywd = "material_space_";
		}
		// print material space data
		if (_myRank==0)
		{
			outfile << endl << "ReadInputFile: Total material space number = " << _matSpaces.size() << endl;
			for (vector<pdSpaceMaterial*>::iterator itor=_matSpaces.begin();itor!=_matSpaces.end();++itor)
			{
				outfile << *(*itor) << endl;
			}
		}
	}

	/* 5. deletion spaces
	*/
	// read deletion spaces data (optional), id starts from 0
	tempKeyword = "number_of_deletion_spaces";
	if (FindKeyword(tempKeyword, OPT))
	{
		ReadKeyword(tempKeyword, 1, _numDelSpaces);
		string keywd = "deletion_space_";
		for (int i=0;i<_numDelSpaces;++i)
		{
			keywd = KeywordAppend(keywd, i+1);
			// Get space geometry type
			int geo_type;
			ReadKeyword(keywd, 1, geo_type);
			switch(geo_type)
			{
			case CUBOID: // type 1 for cuboidal space
				{
					// Get space boundary
					double reg_x1lo, reg_x1hi, reg_x2lo, reg_x2hi, reg_x3lo, reg_x3hi;
					ReadKeyword(keywd, 2, reg_x1lo);
					ReadKeyword(keywd, 3, reg_x1hi);
					ReadKeyword(keywd, 4, reg_x2lo);
					ReadKeyword(keywd, 5, reg_x2hi);
					ReadKeyword(keywd, 6, reg_x3lo);
					ReadKeyword(keywd, 7, reg_x3hi);
					pdSpaceDeletion* delspace = new pdSpaceDeletion(i, geo_type);
					delspace->SetCuboidSpaceBdry(reg_x1lo, reg_x1hi, reg_x2lo, reg_x2hi, reg_x3lo, reg_x3hi);
					_delSpaces.push_back(delspace);
					break;
				}
			case CYLINDER: // type 2 for cylindrical space
				{
					// Get space boundary
					double reg_rad, reg_x1cen, reg_x2cen, reg_x3lo, reg_x3hi;
					ReadKeyword(keywd, 2, reg_rad);
					ReadKeyword(keywd, 3, reg_x1cen);
					ReadKeyword(keywd, 4, reg_x2cen);
					ReadKeyword(keywd, 5, reg_x3lo);
					ReadKeyword(keywd, 6, reg_x3hi);
					pdSpaceDeletion* delspace = new pdSpaceDeletion(i, geo_type);
					delspace->SetCylinderSpaceBdry(reg_rad, reg_x1cen, reg_x2cen, reg_x3lo, reg_x3hi);
					_delSpaces.push_back(delspace);
					break;
				}
			default:
				{
				}
			}
			keywd = "deletion_space_";
		}
		// print deletion space data
		if (_myRank==0)
			{
			outfile << endl << "ReadInputFile: Total deletion space number = " << _delSpaces.size() << endl;
			for (vector<pdSpaceDeletion*>::iterator itor=_delSpaces.begin();itor!=_delSpaces.end();++itor)
			{
				outfile << *(*itor)<< endl;
			}
		}
	}
	else
	{
		_numDelSpaces = 0;
	}

	/* 6. boundary condition space 
	*/
	// read boundary condition spaces data (required), id starts from 0
	tempKeyword = "number_of_boundary_spaces";
	if (FindKeyword(tempKeyword, REQ))
	{
		ReadKeyword(tempKeyword, 1, _numBdSpaces);
		string keywd = "boundary_space_";
		for (int i=0;i!=_numBdSpaces;++i)
		{
			keywd = KeywordAppend(keywd, i+1);
			// Get space geometry type
			int geo_type, bcid;
			ReadKeyword(keywd, 1, geo_type);
			switch (geo_type)
			{
			case CUBOID: // type 1 for cuboidal space
				{
					// Get id of the boundary condition associated with this boundary condition space
					ReadKeyword(keywd, 2, bcid);
					bcid = bcid - 1; // boundary condition id in code starts from zero
					// Get space boundary
					double reg_x1lo, reg_x1hi, reg_x2lo, reg_x2hi, reg_x3lo, reg_x3hi;
					ReadKeyword(keywd, 3, reg_x1lo);
					ReadKeyword(keywd, 4, reg_x1hi);
					ReadKeyword(keywd, 5, reg_x2lo);
					ReadKeyword(keywd, 6, reg_x2hi);
					ReadKeyword(keywd, 7, reg_x3lo);
					ReadKeyword(keywd, 8, reg_x3hi);
					pdBdryCondition* bc = GetBC(bcid);
					pdSpaceBdryCondition* bdspace = new pdSpaceBdryCondition(i, geo_type, bc);
					bdspace->SetCuboidSpaceBdry(reg_x1lo, reg_x1hi, reg_x2lo, reg_x2hi, reg_x3lo, reg_x3hi);
					_bdrySpaces.push_back(bdspace);
					/* 
					
					Store the space limis and strain field of each traction boundary condition (07/19/2011).
					Current assumptions:
					  1. Uniaxial traction, no tangential component.
					  2. Single material in whole body.
					  3. The body is a regular cuboid.
					  4. The surface applied with traction force is the whole boundary surface.

					*/
					if (bc->GetType()==TRACTION)
					{
						double e11, e22, e33;
						// get the body dimension
						double x1Length = _gridDimensionX1*_gridSpacing;
						double x2Length = _gridDimensionX2*_gridSpacing;
						double x3Length = _gridDimensionX3*_gridSpacing;
						// get the Young's modulus and Poisson's ratio
						double ymod = _mats[0]->GetYoungModulus();
						double nu = _mats[0]->GetNu();
						// get the traction boundary condition data
						int dirFlagX1 = bc->GetDirFlagX1();
						int dirFlagX2 = bc->GetDirFlagX2();
						int dirFlagX3 = bc->GetDirFlagX3();
						if (dirFlagX1)
						{
							double f = bc->GetValue(1);							
							double crsArea = x2Length*x3Length;
							e11 = f/(ymod*crsArea);
							e22 = -nu*e11;
							e33 = -nu*e11;
						}
						else if (dirFlagX2)
						{
							double f = bc->GetValue(2);							
							double crsArea = x1Length*x3Length;
							e22 = f/(ymod*crsArea);
							e11 = -nu*e22;							
							e33 = -nu*e22;
						}
						else // dirFlagX3==1
						{
							double f = bc->GetValue(3);							
							double crsArea = x1Length*x2Length;
							e33 = f/(ymod*crsArea);
							e11 = -nu*e33;
							e22 = -nu*e33;							
						}
						double val[12] = {dirFlagX1, dirFlagX2, dirFlagX3, e11, e22, e33, reg_x1lo, reg_x1hi, reg_x2lo, reg_x2hi, reg_x3lo, reg_x3hi};
						vector<double> value(val, val+sizeof(val)/sizeof(double));
						_tractSpaces.push_back(value);
					}
					break;
				}
			case CYLINDER: // type 2 for cylindrical space
				{
					// Get boundary condition id for this space
					ReadKeyword(keywd, 2, bcid);
					bcid = bcid - 1;
					// Get space boundary
					double reg_rad, reg_x1cen, reg_x2cen, reg_x3lo, reg_x3hi;
					ReadKeyword(keywd, 3, reg_rad);
					ReadKeyword(keywd, 4, reg_x1cen);
					ReadKeyword(keywd, 5, reg_x2cen);
					ReadKeyword(keywd, 6, reg_x3lo);
					ReadKeyword(keywd, 7, reg_x3hi);
					pdSpaceBdryCondition* bdspace = new pdSpaceBdryCondition(i, geo_type, GetBC(bcid));
					bdspace->SetCylinderSpaceBdry(reg_rad, reg_x1cen, reg_x2cen, reg_x3lo, reg_x3hi);
					_bdrySpaces.push_back(bdspace);
					break;
				}
			default:
				{
				}
			}
			keywd = "boundary_space_";
		}
		// print boundary space data
		if (_myRank==0)
		{
			outfile << endl << "ReadInputFile: Total boundary space number = " << _bdrySpaces.size() << endl;
			for (vector<pdSpaceBdryCondition*>::iterator itor=_bdrySpaces.begin();itor!=_bdrySpaces.end();++itor)
			{
				outfile << *(*itor)<< endl;
			}
		}
	}

	/* 7. other parameters and control flags 
	*/
	// read linear damping coefficient (required)
	tempKeyword = "damping_coeff";
	if (FindKeyword(tempKeyword, REQ))
	{
        ReadKeyword(tempKeyword, 1, _dampCoeff);
	}
    else
	{
        _dampCoeff = 0.0; // default
	}
	// read force normalization flag (required)
	tempKeyword = "bdryEff_on";
	if (FindKeyword(tempKeyword, REQ))
	{
        ReadKeyword(tempKeyword, 1, _bdryEffOn);
	}
	else
	{
		_bdryEffOn = true; // default
	}
	// Set the global tolerance
	_tolerance = 1.0e-4*_gridSpacing;

	/*
	 Integration methods current in use:
	   1. Integration with fixed Gaussian points (FGI).
	     a. If 1*1*1 Gaussian point is used, moving least square approximation is turned off (_mlsOn=false).
		 b. If more Gaussian points are used, moving least square approximation is turned on (_mlsOn=true).
	   2. Adaptive integration (AI). with or without error control
	     a. If _mlsOn=false (default), then nodal displacments are prescribed.
		 b. If _mlsOn=true, then moving least square approximation is used to calculate the current nodal displacements.
		 c. Error control flag (error_control_on) can be turned on/off. 
	   3. One-point integration, or cubic-cell integration (CCI).
	     a. If _cciMdf=false (default), it is original CCI method where nodes outside the horizon are excluded from the calculation.
		 b. If _cciMdf=true, it is the modified CCI method to include all nodes which have intersecting volume with the horizon.
	*/

	// first check CCI flag (required)
	tempKeyword = "cci_on";
	if (FindKeyword(tempKeyword, REQ))
	{
		ReadKeyword(tempKeyword, 1, _cciOn);
        if (_cciOn)
		{
			// turn off other flags
			_fgiOn = false;
			_mlsOn = false;
			_aiOn = false;		
			_errorControlOn = false;
			_eps = 0.0;
			_intgPoint = 0;
			// read if this CCI is modified (same as !*1*1 Gaussian) or original
			ReadKeyword(tempKeyword, 2, _cciMdf);
		}
		else
		{
			// then check FGI flag (required)
			tempKeyword = "fgi_on";
			if (FindKeyword(tempKeyword, REQ))
			{
				ReadKeyword(tempKeyword, 1, _fgiOn);
				if (_fgiOn)
				{
					// read the number of integration points (Gaussian points) in each direction (optional)
					tempKeyword = "intg_point";
					if (FindKeyword(tempKeyword, OPT))
					{
						ReadKeyword(tempKeyword, 1, _intgPoint);
						if (_intgPoint == 1) 
						{
							// turn off moving least square when using 1*1*1 Gaussian point
							_mlsOn = false;
						}
						else
						{
							_mlsOn = true;
						}
						// read support domain size (optional)
						tempKeyword = "support_domain_size";
						if (FindKeyword(tempKeyword, OPT))
						{
							ReadKeyword(tempKeyword, 1, _rw);
						}
						else
						{
							_rw = 3.0 * _gridSpacing; // default
						}
					}
					else
					{
						_intgPoint = 3; // default
						_mlsOn = true;				
						// read support domain size (optional)
						tempKeyword = "support_domain_size";
						if (FindKeyword(tempKeyword, OPT))
						{
							ReadKeyword(tempKeyword, 1, _rw);
						}
						else
						{
							_rw = 3.0 * _gridSpacing; // default
						}
					}
					// turn off other flags
					_aiOn = false;		
					_errorControlOn = false;
					_eps = 0.0;
					_cciMdf = false;
				}
				else
				{
					_cciMdf = false;
					// at last check AI flag (required)
					tempKeyword = "ai_on";
					if (FindKeyword(tempKeyword, REQ))
					{
						ReadKeyword(tempKeyword, 1, _aiOn);
						if (_aiOn)
						{
							// check the moving least square flag (required)
							ReadKeyword(tempKeyword, 2, _mlsOn);
							// read support domain size (optional)
							tempKeyword = "support_domain_size";
							if (FindKeyword(tempKeyword, OPT))
							{
								ReadKeyword(tempKeyword, 1, _rw);
							}
							else
							{
								_rw = 3.0 * _gridSpacing; // default
							}
							// read trapezoidal index in each direction (optional)
							tempKeyword = "intg_point";
							if (FindKeyword(tempKeyword, OPT))
							{
								// trapezoidal index != number of trapezoidal (integration) points in AI method
								ReadKeyword(tempKeyword, 1, _intgPoint);
							}
							else
							{
								_intgPoint = 3; // default
							}
							// read error control flag (optional)
							tempKeyword = "error_control_on";
							if (FindKeyword(tempKeyword, OPT))
							{
								ReadKeyword(tempKeyword, 1, _errorControlOn);
								if (_errorControlOn)
								{
									// read desired error control accuracy (optional)
									tempKeyword = "error_control_eps";
									if (FindKeyword(tempKeyword, OPT))
									{
										ReadKeyword(tempKeyword, 1, _eps);
									}
									else
									{
										_eps = 0.01; // default
									}
								}
								else
								{
									_eps = 0.0;
								}
							}
							else
							{
								_errorControlOn = false;
								_eps = 0.0;
							}
						}
						else
						{
							// fatal error, since no integration method is used
							cerr << "ReadInputFile: no integration method is used! \n";
							outfile << "ReadInputFile: no integration method is used! \n";
							exit(0);
						}
					}
				}
			}
		}
	}


	// print parameters and control flags data
	outfile << endl << "ReadInputFile: Linear damping coefficient = " << _dampCoeff << endl;
	outfile << endl << "ReadInputFile: Boundary effect compensation factor flag = " << _bdryEffOn << endl;
	outfile << endl << "ReadInputFile: CCI integration method flag = " << _cciOn << endl;
	if (_cciOn)
	{
		outfile << endl << "ReadInputFile: Modified CCI beta flag = " << _cciMdf << endl;
	}
	outfile << endl << "ReadInputFile: Integration with fixed Gaussian points method flag = " << _fgiOn << endl;
	if (_fgiOn)
	{
		outfile << endl << "ReadInputFile: Starting Gaussian points = " << _intgPoint << endl;
		outfile << endl << "ReadInputFile: Moving least square flag = " << _mlsOn << endl;
	}
	outfile << endl << "ReadInputFile: Adaptive integration method flag = " << _aiOn << endl;
	if (_aiOn)
	{
		outfile << endl << "ReadInputFile: Moving least square flag = " << _mlsOn << endl;
		outfile << endl << "ReadInputFile: Error control method flag = " << _errorControlOn << endl;
		outfile << endl << "ReadInputFile: Error control accuracy = " << _eps << endl;
		outfile << endl << "ReadInputFile: Starting trapezoidal index = " << _intgPoint << endl;
	}

	// close input file
	outfile.close();
}

void pdDataManager::InitMPI() 
{
	// This function initial mpi environment ( Kebing 08/10/09)
	
	/*MPI::Init();
	_myRank = MPI::COMM_WORLD.Get_rank();
	_numProcs = MPI::COMM_WORLD.Get_size();*/
	_myRank = 0;
	_numProcs = 1;
	// get working folder path
	SetWkDir();
	ofstream outfile(_outFilePath.c_str(), ios::app);		
	//// read processor data (required)
	//ReadKeyword("number_of_processors", 1, _proc1);
	//ReadKeyword("number_of_processors", 2, _proc2);
	//ReadKeyword("number_of_processors", 3, _proc3);
	//if (_numProcs < _proc1*_proc2*_proc3)
	//{
	//	if (_myRank==0)
	//	{
	//		cerr << "InitMPI: wrong allocation of processors!" << endl;
	//		outfile << "InitMPI: wrong allocation of processors!" << endl;
	//	}
	//	exit(0);
	//}
	outfile.close();
}

void pdDataManager::InitLocalBonds()
{
	// This function allocates local bonds to each processor
	// try to assign equal number of bonds when the total number of bonds can not be divided by number of processors
	// example: total bonds=30, total procs=4, will allocate 0-7(8), 8-15(8), 16-22(7), 23-29(7) bonds to proc0, 1, 2, 3

	int s, e; // start and end of bond number
	for (int ip=0;ip!=_numProcs;++ip)
	{
		int incre = _numBonds/_numProcs;
		int rmd = _numBonds%_numProcs; // the remainder
		s = ip*incre + min(ip, rmd);
		if ( ip< rmd)
		{
			++incre;
		}
		e = s + incre -1;
		if (ip==_myRank)
		{
			_myBonds.resize(incre);
			copy(_bonds.begin()+s, _bonds.begin()+e+1, _myBonds.begin());
		}
	}
}

void pdDataManager::InitNode()
{
	/* 

	This function initials global nodes, including:
     1. Generating node coordinates and find total number of nodes after deletion.
     2. Initial material properties associated with each node.
     3. Initial boundary conditions associated with each node.
     4. Initial nodal displacements, velocities and other variables.
	 
	 */

	ofstream outfile(_outFilePath.c_str(), ios::app);
	double tol = _tolerance;
	// read nodal data from external grid file (optional)
	_extGridFileOn = false;
	string fpath;
	// check if keyword "grid_file" is in the input file
	if (FindKeyword("grid_file", OPT))
	{
		_extGridFileOn = true;
		ReadKeyword("grid_file", 1, fpath);
		ifstream ext_input(fpath.c_str(), ios::in);
		ReadExternalGridFile(ext_input);
		if (_myRank==0)
		{
			outfile << endl << "ReadInputFile: external grid file " << fpath << " is read" << endl;
		}
	}
	else
	{
		// generate node coordinates if no external grid file.
		int id = 0;
		for (int n=1;n!=_gridDimensionX3+1;++n) // z direction
		{
			for (int m=1;m!=_gridDimensionX2+1;++m) // y direction
			{
				for (int l=1;l!=_gridDimensionX1+1;++l) // x direction
				{
					// for cubic lattice only 
					double x1i = (l - (_gridDimensionX1 + 1) * 0.5) * _gridSpacing + _gridCenterX1;
					double x2i = (m - (_gridDimensionX2 + 1) * 0.5) * _gridSpacing + _gridCenterX2;
					double x3i = (n - (_gridDimensionX3 + 1) * 0.5) * _gridSpacing + _gridCenterX3;
					// check if this node is in a deletion space
					bool idel = false;
					for (int d=0;d!=_numDelSpaces;++d)
					{
						// find deletion space geometry type
						pdSpaceDeletion* dreg = GetDelSpace(d);
						int dgeo_type = dreg->GetGeomType();
						switch (dgeo_type)
						{
						case CUBOID:
							{
								double x1lo = dreg->GetX1Low();
								double x1hi = dreg->GetX1High();
								double x2lo = dreg->GetX2Low();
								double x2hi = dreg->GetX2High();
								double x3lo = dreg->GetX3Low();
								double x3hi = dreg->GetX3High();
								if (x1i>=x1lo-tol && x1i<=x1hi+tol && x2i>=x2lo-tol && x2i<=x2hi+tol
									&& x3i>=x3lo-tol && x3i<=x3hi+tol)
								{
									idel = true;
								}
								break;
							}
						case CYLINDER:
							{
								double rad   = dreg->GetRadius();
								double x1cen = dreg->GetX1Center();
								double x2cen = dreg->GetX2Center();
								double x3lo  = dreg->GetX3Low();
								double x3hi  = dreg->GetX3High();
								double dist  = sqrt((x1i-x1cen)*(x1i-x1cen)+(x2i-x2cen)*(x2i-x2cen));
								if (dist<=rad+tol && x3i>=x3lo-tol && x3i<=x3hi+tol)
								{
									idel = true;
								}
								break;
							}
						default:
							{
							}
						}
					}
					if (!idel)
					{
						pdNode* node_i = new pdNode(id++); // node id starts from 0
						node_i->SetX1(x1i);
						node_i->SetX2(x2i);
						node_i->SetX3(x3i);
						_nodes.push_back(node_i);						
					}					
				}
			}
		}
		_numNodes = int(_nodes.size());
		if (_myRank==0)
		{
			outfile << endl << "InitNode: total node number = " << _numNodes << endl;
		}
	}
	// find material associated with each node by checking if node falls into one of the material spaces.
	// for overlap space, material properties of the latter material space will be applied.
	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{
		pdNode* node_i = *itor;
		double x1i = node_i->GetX1();
		double x2i = node_i->GetX2();
		double x3i = node_i->GetX3();
		// loop over material spaces
		for (vector<pdSpaceMaterial*>::const_iterator itor_m=_matSpaces.begin();itor_m!=_matSpaces.end();++itor_m)
		{
			// find space geometry type and material id
			pdSpaceMaterial* mreg = *itor_m;
			int mgeo_type = mreg->GetGeomType();
			pdMaterial* mat = mreg->GetMaterial();
			switch (mgeo_type)
			{
			case CUBOID:
				{
					double x1lo = mreg->GetX1Low();
					double x1hi = mreg->GetX1High();
					double x2lo = mreg->GetX2Low();
					double x2hi = mreg->GetX2High();
					double x3lo = mreg->GetX3Low();
					double x3hi = mreg->GetX3High();
					if (x1i>=x1lo-tol && x1i<=x1hi+tol && x2i>=x2lo-tol && x2i<=x2hi+tol
						&& x3i>=x3lo-tol && x3i<=x3hi+tol)
					{
						node_i->SetMaterial(mat);
					}
					else
					{
						node_i->SetMaterial(0);						
					}
					break;
				}
			case CYLINDER:
				{
					double rad   = mreg->GetRadius();
					double x1cen = mreg->GetX1Center();
					double x2cen = mreg->GetX2Center();
					double x3lo  = mreg->GetX3Low();
					double x3hi  = mreg->GetX3High();
					double dist  = sqrt((x1i-x1cen)*(x1i-x1cen)+(x2i-x2cen)*(x2i-x2cen));
					if (dist<=rad+tol && x3i>=x3lo-tol && x3i<=x3hi+tol)
					{
						node_i->SetMaterial(mat);
					}
					else
					{
						node_i->SetMaterial(0);						
					}
					break;
				}
			default:
				{
				}
			}
		}
	}
	// check if all nodes have material associated with it
	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{
		pdNode* node_i = *itor;
		if (node_i->GetMaterial()==0) // fatal error
		{
			if (_myRank==0)
			{
				outfile << endl << "InitNode: node " << node_i->GetID() << " not in any material space!" << endl
					<< "InitNode: node coordinate: x1 = " << node_i->GetX1()
					<< ", x2 = " << node_i->GetX2() << ", x3 = " << node_i->GetX3() << endl;
				cerr << endl << "InitNode: node " << (*itor)->GetID() << " not in any material space!" << endl
					<< "InitNode: node coordinate: x1 = " << node_i->GetX1()
					<< ", x2 = " << node_i->GetX2() << ", x3 = " << node_i->GetX3() << endl;
			}
			exit(0);
		}
	}
	if (_myRank==0)
	{
		outfile << endl << "InitNode: Material association done " << endl;
	}

	// find boundary condition associated with each node by checking if node falls into one of 
	// the boundary condition spaces.
	// for overlap space, the boundary condition of latter boundary space will be applied.
	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{		
		pdNode* node_i = *itor;
		double x1i = node_i->GetX1();
		double x2i = node_i->GetX2();
		double x3i = node_i->GetX3();
		// default: not in any boundary condition space (set to 0)
		node_i->SetBC(0);
		// loop over boundary spaces
		for (vector<pdSpaceBdryCondition*>::const_iterator itor_b=_bdrySpaces.begin();itor_b!=_bdrySpaces.end();++itor_b)
		{
			// find space geometry type and boundary condition id
			pdSpaceBdryCondition* breg = *itor_b;
			int bgeo_type = breg->GetGeomType();
			pdBdryCondition* bc = breg->GetBC();
			switch (bgeo_type)
			{
			case CUBOID:
				{
					double x1lo = breg->GetX1Low();
					double x1hi = breg->GetX1High();
					double x2lo = breg->GetX2Low();
					double x2hi = breg->GetX2High();
					double x3lo = breg->GetX3Low();
					double x3hi = breg->GetX3High();
					if (x1i>=x1lo-tol && x1i<=x1hi+tol && x2i>=x2lo-tol && x2i<=x2hi+tol
						&& x3i>=x3lo-tol && x3i<=x3hi+tol)
					{
						node_i->SetBC(bc);
					}
					break;
				}
			case CYLINDER:
				{
					double rad   = breg->GetRadius();
					double x1cen = breg->GetX1Center();
					double x2cen = breg->GetX2Center();
					double x3lo  = breg->GetX3Low();
					double x3hi  = breg->GetX3High();
					double dist  = sqrt((x1i-x1cen)*(x1i-x1cen)+(x2i-x2cen)*(x2i-x2cen));
					if (dist<=rad+tol && x3i>=x3lo-tol && x3i<=x3hi+tol)
					{
						node_i->SetBC(bc);
					}
					break;
				}
			default:
				{
				}
			}
		}
	}
	if (_myRank==0)
	{
		outfile << endl << "InitNode: boundary condition association done " << endl;
	}


	// initial other variables
	int nrow = 0; // the row number for 1st dof of each node, starting from 0
	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{
		pdNode* node_i = *itor;
		// set initial displacement and velocity to zero, initial velocity buff to zero
		node_i->SetU1(0.0);
		node_i->SetU2(0.0);
		node_i->SetU3(0.0);
		node_i->SetV1(0.0);
		node_i->SetV2(0.0);
		node_i->SetV3(0.0);
		node_i->SetV1Buffer(0.0);
		node_i->SetV2Buffer(0.0);
		node_i->SetV3Buffer(0.0);
		// initial force normalization factor to 1
		node_i->SetBdryEffFactor(1.0);
		// initial bond damage, yield faction to zero
		node_i->SetNumDmgBond(0);
		node_i->SetNumYldBond(0);
		// set node radius and volume, initial force normalization factor to 1 (default)
		node_i->SetNodeVolume(pow(_gridSpacing, 3.0));
		// set support domain size
        node_i->SetSupportDomainSize(_rw);
		// store the row number for 1st dof of each node
		_rowOfNode.push_back(nrow);
		nrow += 3;
	}
	if (_myRank==0)
	{
		outfile << endl << "InitNode: nodes initialization done " << endl;
	}
	
	// initial the total work done through boundary and elastic energy to zero
	_externalWorkTotal = 0.0;
	_elasticEnergyTotal = 0.0;
	_externalWorkStep = 0.0;
	outfile.close();
}

void pdDataManager::InitGridBoundary()
{
	// This function initialize the grid boundry. 
	// It must be called after InitNode. Also initial force-through-plane grid.

	ofstream outfile(_outFilePath.c_str(), ios::app);
	double x1a, x1b, x2a, x2b, x3a, x3b;
	if (_extGridFileOn)
	{
		x1a = pow(10.0, 10.0);
		x1b = -pow(10.0, 10.0);
		x2a = pow(10.0, 10.0);
		x2b = -pow(10.0, 10.0);
		x3a = pow(10.0, 10.0);
		x3b = -pow(10.0, 10.0);
		for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
		{
			pdNode* node_i = *itor;
			double x1 = node_i->GetX1() + node_i->GetU1();
			double x2 = node_i->GetX2() + node_i->GetU2();
			double x3 = node_i->GetX3() + node_i->GetU3();
			x1a = min(x1a, x1);
			x1b = max(x1b, x1);
			x2a = min(x2a, x2);
			x2b = max(x2b, x2);
			x3a = min(x3a, x3);
			x3b = max(x3b, x3);
		}
	}
	else
	{
		x1a = -(_gridDimensionX1-1)*0.5*_gridSpacing + _gridCenterX1;
		x1b =  (_gridDimensionX1-1)*0.5*_gridSpacing + _gridCenterX1;
		x2a = -(_gridDimensionX2-1)*0.5*_gridSpacing + _gridCenterX2;
		x2b =  (_gridDimensionX2-1)*0.5*_gridSpacing + _gridCenterX2;
		x3a = -(_gridDimensionX3-1)*0.5*_gridSpacing + _gridCenterX3;
		x3b =  (_gridDimensionX3-1)*0.5*_gridSpacing + _gridCenterX3;
	}
	// set the grid boundary for the default deformed grid.
	// this is the outer bounding.
	double deformedBdry = 50.0 * _gridSpacing; // a default deformed grid boundary 
	_gridBoundaryX1Low = x1a - deformedBdry;
	_gridBoundaryX1High = x1b + deformedBdry;
	_gridBoundaryX2Low = x2a - deformedBdry;
	_gridBoundaryX2High = x2b + deformedBdry;
	_gridBoundaryX3Low = x3a - deformedBdry;
	_gridBoundaryX3High = x3b + deformedBdry;
	// read locations of planes (optional) for total force calculation
	if (FindKeyword("total_force_plane", OPT))
	{
		ReadKeyword("total_force_plane", 1, _forcePlaneX1);
		ReadKeyword("total_force_plane", 2, _forcePlaneX2);
		ReadKeyword("total_force_plane", 3, _forcePlaneX3);		
	}
	else
	{
		// default position is xy, yz, zx plane
		_forcePlaneX1 = 0.5*(x1a+x1b);
		_forcePlaneX2 = 0.5*(x2a+x2b);
		_forcePlaneX3 = 0.5*(x3a+x3b);
	}
	if (_myRank==0)
	{
		outfile << endl << "InitGridBoundary: grid boundaries: x1a = " << x1a << ", x1b = " << x1b << endl
			<< "InitGridBoundary: grid boundaries: x2a = " << x2a << ", x2b = " << x2b << endl
			<< "InitGridBoundary: grid boundaries: x3a = " << x3a << ", x3b = " << x3b << endl
			<< endl << "InitGridBoundary: total force planes: x1= " << _forcePlaneX1
			<< ", x2= " << _forcePlaneX2 << ", x3= " << _forcePlaneX3 << endl
			<< endl << "InitGridBoundary: model space boundaries: x1clo = " << _gridBoundaryX1Low << ", x1chi = " << _gridBoundaryX1High << endl
			<< "InitGridBoundary: model space boundaries: x2clo = " << _gridBoundaryX2Low << ", x2chi = " << _gridBoundaryX2High << endl
			<< "InitGridBoundary: model space boundaries: x3clo = " << _gridBoundaryX3Low << ", x3chi = " << _gridBoundaryX3High << endl;
	}
	//// initial force-through-plane grid (x-z plane current only 10/06/08)
	//int id = 0;
	//for (int j=0;j!=_gridDimensionX1;++j)
	//{
	//	for (int k=0;k!=_gridDimensionX3;++k)
	//	{
	//		double x1l = -(0.5*_gridDimensionX1-j)*_gridSpacing;
	//		double x1h = -(0.5*_gridDimensionX1-j)*_gridSpacing + _gridSpacing;
	//		double x2l = 0.;
	//		double x2h = 0.;
	//		double x3l = -(0.5*_gridDimensionX3-k)*_gridSpacing;
	//		double x3h = -(0.5*_gridDimensionX3-k)*_gridSpacing + _gridSpacing;
	//		pdSpacePlane* plane = new pdSpacePlane(id, 1);
	//		plane->SetCuboidSpaceBdry(x1l, x1h, x2l, x2h, x3l, x3h);
	//		_planes.push_back(plane);
	//		++id;
	//	}
	//}
	//// loop through each bond to find how many planes it passes (10/13/08)
	//// this is done once before dynamic calculation begins
	//for (vector<pdBond*>::const_iterator itor=_bonds.begin();itor!=_bonds.end();++itor)
	//{
	//	// loops through all force-through-planes to see if this bond is pass through
	//	// If so, calculate the force normal to the plane
	//	// i and j is two ends of a bond
	//	pdBond* bond_ij = *itor;
	//	pdNode* node_i = bond_ij->GetNodeI();
	//	pdNode* node_j = bond_ij->GetNodeJ();
	//	double x1i = node_i->GetX1(); //Get coord of node i and j
	//	double x2i = node_i->GetX2();
	//	double x3i = node_i->GetX3();
	//	double x1j = node_j->GetX1(); 
	//	double x2j = node_j->GetX2();
	//	double x3j = node_j->GetX3();
	//	for (vector<pdSpacePlane*>::const_iterator itor_p=_planes.begin();itor_p!=_planes.end();++itor_p)
	//	{
	//		pdSpacePlane* i_plane = *itor_p;
	//		// get plane id and position
	//		double x1lo = i_plane->GetX1Low();
	//		double x1hi = i_plane->GetX1High();
	//		double x2lo = i_plane->GetX2Low();
	//		double x2hi = i_plane->GetX2High();
	//		double x3lo = i_plane->GetX3Low();
	//		double x3hi = i_plane->GetX3High();
	//		// get plane normal direction
	//		double xi, xj, xplane;
	//		int dir = abs(i_plane->GetNorm());
	//		switch (dir)
	//		{
	//		case 1:
	//			{
	//				xi = x1i;
	//				xj = x1j;
	//				xplane = x1lo;
	//				break;
	//			}
	//		case 2:
	//			{
	//				xi = x2i;
	//				xj = x2j;
	//				xplane = x2lo;
	//				break;
	//			}
	//		case 3:
	//			{
	//				xi = x3i;
	//				xj = x3j;
	//				xplane = x3lo;
	//				break;
	//			}
	//		default:
	//			{
	//			}
	//		}
	//		// first check if node i and j are not at the same side of the plane
	//		if ((xi-xplane)*(xj-xplane)<0.)
	//		{
	//			// then check if this bond is passing through the plane
	//			// p1, p2, p3 is the coordinates of intersection
	//			double p1, p2, p3;
	//			switch (dir)
	//			{
	//			case 1:
	//				{
	//					p1 = xplane;
	//					p2 = x2i + (x2j-x2i)*(x1lo-x1i)/(x1j-x1i);
	//					p3 = x3i + (x3j-x3i)*(x1lo-x1i)/(x1j-x1i);
	//					if ( (p2-x2lo)*(p2-x2hi)<=0. && (p3-x3lo)*(p3-x3hi)<=0. )
	//					{
	//						bond_ij->AddPlane(i_plane);
	//					}
	//					break;
	//				}
	//			case 2:
	//				{
	//					p1 = x1i + (x1j-x1i)*(x2lo-x2i)/(x2j-x2i);
	//					p2 = xplane;
	//					p3 = x3i + (x3j-x3i)*(x2lo-x2i)/(x2j-x2i);
	//					if ( (p1-x1lo)*(p1-x1hi)<=0. && (p3-x3lo)*(p3-x3hi)<=0. )
	//					{
	//						bond_ij->AddPlane(i_plane);
	//					}
	//					break;
	//				}
	//			case 3:
	//				{
	//					p1 = x1i + (x1j-x1i)*(x3lo-x3i)/(x3j-x3i);
	//					p2 = x2i + (x2j-x2i)*(x3lo-x3i)/(x3j-x3i);
	//					p3 = xplane;
	//					if ( (p2-x2lo)*(p2-x2hi)<=0. && (p1-x1lo)*(p1-x1hi)<=0. )
	//					{
	//						bond_ij->AddPlane(i_plane);
	//					}
	//					break;
	//				}
	//			default:
	//				{
	//				}
	//			}
	//		}
	//	}
	//}
	outfile.close();
}

void pdDataManager::InitMaterial()
{
	// This function calculates the spring constant for each material. It must be called after InitNode.
	// Also find maximum spring constant and minimum denstiy.

	ofstream outfile(_outFilePath.c_str(), ios::app);
	if (_myRank==0)
	{
		outfile << endl << "InitMatieral: Total material number = " << _mats.size() << endl;
	}
	for (vector<pdHandle<pdMaterial>>::iterator itor=_mats.begin();itor!=_mats.end();++itor)
	{
		pdMaterial* mat = (*itor).GetPtr();
		int mat_type = mat->GetType();
		// type 1 for bond-based material.
		if (mat_type==BOND)
		{
			double spring = 18*mat->GetBulkModulus()/(PI*pow(_horizon, 4.0));
			mat->SetSpringConstant(spring);
		}
		// type 2 for state-based material.
		else if (mat_type==STATE)
		{
			mat->SetSpringConstant(0.0); // sprint constant is not necessary in state-based material
		}
		// print material data
		if (_myRank==0)
		{
			outfile << *mat << endl;
		}
	}
	// find maximum spring constant and minimum denstiy
	// set max spring constant to zero for cycle 1 stalbe time step calculation
	_minMaterialDensity = 1.0e10;
	_maxBulkSoundSpeed = 0.0;
	for (vector<pdHandle<pdMaterial>>::iterator itor=_mats.begin();itor!=_mats.end();++itor)
	{
		_minMaterialDensity = min(_minMaterialDensity, (*itor)->GetDensity());
		_maxBulkSoundSpeed = max(_maxBulkSoundSpeed, (*itor)->GetBulkSoundSpeed());
	}		
	outfile.close();
}

void pdDataManager::InitBdryEffectFactor()
{
	// This function calculates the boundary effect compensation factor for each node.
	// It has to be called after InitMaterial() since the calculation needs to know the value of spring constant.

	if (_bdryEffOn)
	{
		// 1. find the total number of family nodes in a full horizon sphere
		// create a virtual grid center at the origin
		int vd = 2 * (int(_horizon/_gridSpacing) + 1); 
		pdNode* node_i = new pdNode(-1); // virtual source node
		node_i->SetX1(0.0);
		node_i->SetX2(0.0);
		node_i->SetX3(0.0);
		node_i->SetU1(0.0);
		node_i->SetU2(0.0);
		node_i->SetU3(0.0);
		_numMaxFamily = 0;
		for (int n=1;n!=vd+1;++n) // z direction
		{
			for (int m=1;m!=vd+1;++m) // y direction
			{
				for (int l=1;l!=vd+1;++l) // x direction
				{
					// generate the coordinate of the virtual nodes in the grid
					double x1i = (l-vd*0.5)*_gridSpacing;
					double x2i = (m-vd*0.5)*_gridSpacing;
					double x3i = (n-vd*0.5)*_gridSpacing;
					// create virtual node j
					pdNode* node_j = new pdNode(-1);
					node_j->SetX1(x1i);
					node_j->SetX2(x2i);
					node_j->SetX3(x3i);
					node_j->SetU1(0.0);
					node_j->SetU2(0.0);
					node_j->SetU3(0.0);
					// check if it is a family node of the source node (0, 0, 0)
					double dist_ij = sqrt(x1i*x1i + x2i*x2i + x3i*x3i);
					// find the shortest and longest distance from node i to the cubic volume of node j
					double dist_min, dist_max;
					CalcDistanceToCube(node_i, node_j, dist_min, dist_max);
					if (dist_min < _horizon + _tolerance) // node j's cell has intersection volume with horizon sphere
						++_numMaxFamily;
					delete node_j;
				}
			}
		}
		delete node_i;

		// 2. use virtual 3*3*3 Gaussian points to calculate the compensation factor
		// backup the real number of integration points
		int intgBuffer = _intgPoint;
		_intgPoint = 3;
		// because this function is called before InitIntegrationPoint()
		// must set up a virtual 3*3*3 Gaussian points array for each node
		SetGaussianPoint();
		// find the real body boundaries (accurate for cubic grid)
		double x1Low = -_gridDimensionX1*0.5*_gridSpacing + _gridCenterX1;
		double x1High = _gridDimensionX1*0.5*_gridSpacing + _gridCenterX1;
		double x2Low = -_gridDimensionX2*0.5*_gridSpacing + _gridCenterX2;
		double x2High = _gridDimensionX2*0.5*_gridSpacing + _gridCenterX2;
		double x3Low = -_gridDimensionX3*0.5*_gridSpacing + _gridCenterX3;
		double x3High = _gridDimensionX3*0.5*_gridSpacing + _gridCenterX3;
		double ns = (int(_horizon/_gridSpacing)) * _gridSpacing;
		for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
		{
			pdNode* node = *itor;
			double x1i = node->GetX1();
			double x2i = node->GetX2();
			double x3i = node->GetX3();
			pdMaterial* mat = node->GetMaterial();
			// 1) nodes with full horizon sphere
			if (node->GetNumFamily()==_numMaxFamily) 
			{
				node->SetBdryEffFactor(1.0);
			}
			else
			{
				// 2) boundary nodes whose horizon spheres intersecting with only one surface (type one)
				// X1
				if (x1i>=x1High-ns && x1i<=x1High && x2i>=x2Low+ns && x2i<=x2High-ns
					&& x3i>=x3Low+ns && x3i<=x3High-ns)
				{				
					double h = x1High - x1i;
					node->SetBdryEffFactor(this->CalcBdryEffFactorTypeOne(h, mat));
				}
				else if (x1i>=x1Low && x1i<=x1Low+ns && x2i>=x2Low+ns && x2i<=x2High-ns
					&& x3i>=x3Low+ns && x3i<=x3High-ns)
				{				
					double h = x1i - x1Low;
					node->SetBdryEffFactor(this->CalcBdryEffFactorTypeOne(h, mat));
				}
				// X2
				else if (x2i>=x2High-ns && x2i<=x2High && x1i>=x1Low+ns && x1i<=x1High-ns
					&& x3i>=x3Low+ns && x3i<=x3High-ns)
				{				
					double h = x2High - x2i;
					node->SetBdryEffFactor(this->CalcBdryEffFactorTypeOne(h, mat));
				}
				else if (x2i>=x2Low && x2i<=x2Low+ns && x1i>=x1Low+ns && x1i<=x1High-ns
					&& x3i>=x3Low+ns && x3i<=x3High-ns)
				{				
					double h = x2i - x2Low;
					node->SetBdryEffFactor(this->CalcBdryEffFactorTypeOne(h, mat));
				}
				// X3
				else if (x3i>=x3High-ns && x3i<=x3High && x1i>=x1Low+ns && x1i<=x1High-ns
					&& x2i>=x2Low+ns && x2i<=x2High-ns)
				{				
					double h = x3High - x3i;
					node->SetBdryEffFactor(this->CalcBdryEffFactorTypeOne(h, mat));
				}
				else if (x3i>=x3Low && x3i<=x3Low+ns && x1i>=x1Low+ns && x1i<=x1High-ns
					&& x2i>=x2Low+ns && x2i<=x2High-ns)
				{				
					double h = x3i - x3Low;
					node->SetBdryEffFactor(this->CalcBdryEffFactorTypeOne(h, mat));
				}
				else
				{
					// 3) any other nodes
					CalcBdryEffFactor(node);
				}
			}
		}
		// restore the initial number of integration points 
		// and clear the Gaussian points list for each node
		_intgPoint = intgBuffer;
		for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
		{
			(*itor)->ClearGaussianPointList();	
		}
	}	
}

//void pdDataManager::CalcStateBasedMotion(const double dt)
//{
//	// This function updates the node velocity due to peridynamic state interactions.
//
//	double tol = _tolerance;
//	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
//	{
//		pdNode* node_i = *itor;
//		pdMaterial* mat = node_i->GetMaterial();
//		if (mat!=0 && mat->GetType()==2) // a state-based material
//		{
//			// find the unbroken bonds, total nodes and total velocity
//			int num_nodes_in_state = 0;
//			int nbond_not_broken = 0;
//			double sum_v1 = 0.;
//			double sum_v2 = 0.;
//			double sum_v3 = 0.;
//			int numfam = node_i->GetNumFamily();
//			// loop over family nodes
//			for (int n=0;n!=numfam;++n)
//			{
//				pdNode* node_j = node_i->GetFamilyNode(n);
//				pdBond* bond_ij = node_i->GetFamilyBond(n);
//				if (!bond_ij->IsBreak_state())
//				{
//					++nbond_not_broken;
//					++num_nodes_in_state;
//					sum_v1 += node_j->GetV1();
//					sum_v2 += node_j->GetV2();
//					sum_v3 += node_j->GetV3();
//				}
//			}
//			// find the damage and mean velocity at each node
//			if (nbond_not_broken>0)
//			{
//				vel_family[0] = sum_v1/nbond_not_broken;
//				vel_family[1] = sum_v2/nbond_not_broken;
//				vel_family[2] = sum_v3/nbond_not_broken;
//				double damage_state = 1. - nbond_not_broken/(num_nodes_in_state+1.0e-10); // a local damage_state ??
//			}
//			else
//			{
//				vel_family[0] = 0.;
//				vel_family[1] = 0.;
//				vel_family[2] = 0.;
//				double damage_state = 1.; // a local damage_state ??
//			}
//			// find force density with its family nodes
//			CalcStateBasedForce(node_i, dt);
//			// apply damping to node mi to bring its velocity closer to the mean velocity within its family.
//			// this is done only after some bonds are broken.
//			if (num_nodes_in_state>0)
//			{
//				double damage_state = 1. - nbond_not_broken/(num_nodes_in_state+1.0e-10);
//				double fac = 0.5*damage_state;
//				double v1 = (1. - fac)*node_i->GetV1Buffer() + fac*vel_family[0];
//				double v2 = (1. - fac)*node_i->GetV2Buffer() + fac*vel_family[1];
//				double v3 = (1. - fac)*node_i->GetV3Buffer() + fac*vel_family[2];
//				node_i->SetV1Buffer(v1);
//				node_i->SetV2Buffer(v2);
//				node_i->SetV3Buffer(v3);
//			}
//			// sum force density at each node due to all state interactions.
//			// this includes terms due the force states both at node mi and mj.
//			double x1i = node_i->GetX1();
//			double x2i = node_i->GetX2();
//			double x3i = node_i->GetX3();
//			double voli = node_i->GetNodeVolume();
//			for (int n=0;n!=numfam;++n)
//			{
//				pdNode* node_j = node_i->GetFamilyNode(n);
//				pdBond* bond_ij = node_i->GetFamilyBond(n);
//				if (!bond_ij->IsBreak_state())
//				{
//					double volj = node_j->GetNodeVolume();
//					double force_dens_1i =  bond_ij->GetForce_state_1()*volj;
//					double force_dens_2i =  bond_ij->GetForce_state_2()*volj;
//					double force_dens_3i =  bond_ij->GetForce_state_3()*volj;
//					double force_dens_1j = -bond_ij->GetForce_state_1()*voli;
//					double force_dens_2j = -bond_ij->GetForce_state_2()*voli;
//					double force_dens_3j = -bond_ij->GetForce_state_3()*voli;
//					// update the force density for node mi and mj
//					node_i->AddForce(force_dens_1i, force_dens_2i, force_dens_3i);
//					node_j->AddForce(force_dens_1j, force_dens_2j, force_dens_3j);
//					// find total normal force through force planes
//					double x1j = node_j->GetX1();
//					double x2j = node_j->GetX2();
//					double x3j = node_j->GetX3();
//					if (x1j>_forcePlaneX1-tol && x1i<=_forcePlaneX1+tol)
//					{
//						_forceX1 += bond_ij->GetForce_state_1()*volj*voli;
//					}
//					else if (x1i>_forcePlaneX1-tol && x1j<=_forcePlaneX1+tol)
//					{
//						_forceX1 -= bond_ij->GetForce_state_1()*volj*voli;
//					}
//					if (x2j>_forcePlaneX2-tol && x2i<=_forcePlaneX2+tol)
//					{
//						_forceX2 += bond_ij->GetForce_state_2()*volj*voli;
//					}
//					else if (x2i>_forcePlaneX2-tol && x2j<=_forcePlaneX2+tol)
//					{
//						_forceX2 -= bond_ij->GetForce_state_2()*volj*voli;
//					}
//					if (x3j>_forcePlaneX3-tol && x3i<=_forcePlaneX3+tol)
//					{
//						_forceX3 += bond_ij->GetForce_state_3()*volj*voli;
//					}
//					else if (x3i>_forcePlaneX3-tol && x3j<=_forcePlaneX3+tol)
//					{
//						_forceX3 -= bond_ij->GetForce_state_3()*volj*voli;
//					}
//				}
//			}
//		}
//	}	
//	// update velocity based on the force densities due to state interactions.
//	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
//	{
//		pdNode* node_i = *itor;
//		double denst = node_i->GetMaterial()->GetDensity();
//		double v1 = node_i->GetV1Buffer() + dt*node_i->GetForce1()/denst;
//		double v2 = node_i->GetV2Buffer() + dt*node_i->GetForce2()/denst;
//		double v3 = node_i->GetV3Buffer() + dt*node_i->GetForce3()/denst;
//		node_i->SetV1Buffer(v1);
//		node_i->SetV2Buffer(v2);
//		node_i->SetV3Buffer(v3);		
//	}
//}

void pdDataManager::CalcBondBasedMotion(const double dt, const double iter)
{
	/*

	This function updates the node velocity due to peridynamice force.
	Input:
	  1. dt: the current time step.
	  2. iter: the current iteration.

	*/

	double tol = _tolerance;
	double vBuffer[3]; // store update node velocity
	double pdForce[3]; // store peridynamic force between node i and j
	double bodyForce[3]; // store body force on node i
	int numIntg;

	// for MPI
	//MPI::Request send_req_i, send_req_j;
	//double vsend_i[3], vsend_j[3]; 
	
	//// first calculate the damping coefficient at current time step
	//// also get the fictitious diagonal mass matrix
	//pdVector* mMatrix = new pdVector(_numNodes*3); // only diagonal terms are stored
	//double dampCoeff = CalcDampingCoeff(dt, mMatrix);
	double pdForceSum[3] = {0.0, 0.0, 0.0};
	//string ofpath = _workDir + "pds_debug.dat";
	//ofstream fout(ofpath.c_str(), ios::app);

	/*** loop over nodes starts ***/
	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{
		pdNode* node_i = *itor;
		pdMaterial* mat = node_i->GetMaterial();
		double denst = mat->GetDensity();
		// get body force
		bodyForce[0] = mat->GetGrav1()*denst;
		bodyForce[1] = mat->GetGrav2()*denst;
		bodyForce[2] = mat->GetGrav3()*denst;
		// loop over family of node i
		int numFam = node_i->GetNumFamily();
		for (int j=0;j!=numFam;++j)
		{
			pdNode* node_j = node_i->GetFamilyNode(j);
			pdBond* bond_ij = node_i->GetFamilyBond(j);
			if (!bond_ij->IsBreak())  // check if this bond is already broken
			{
				// force calculation
				if (_cciOn)
				{
					/*** CCI method ***/
					CalcBondForceCCI(bond_ij, node_i, node_j, pdForce, numIntg);			
				}
				else
				{
					/*** Adaptive trapezoidal integration or Gaussian quadrature integration ***/
					CalcBondForceIntg(bond_ij, node_i, node_j, pdForce, numIntg); 
				}
				// get the update to node velocity buffer from this bond					
				for (int k=0;k!=3;++k)
				{
					vBuffer[k] = (pdForce[k] + bodyForce[k]) * dt / denst;						
				}
				// update node velocity buffer
				node_i->AddV1Buffer(vBuffer[0]);
				node_i->AddV2Buffer(vBuffer[1]);
				node_i->AddV3Buffer(vBuffer[2]);
			}
			else
			{
				pdForce[0] = 0.0;
				pdForce[1] = 0.0;
				pdForce[2] = 0.0;
				node_i->AddNumDmgBond(1);
				node_i->AddNumYldBond(1);
				node_i->AddElastEnergyDensity(0.0);
			}

			//// add the forces from this bond to the total forces
			//pdForceSum[0] += pdForce[0];
			//pdForceSum[1] += pdForce[1];
			//pdForceSum[2] += pdForce[2];			
		}
		//// get the update to node velocity buffer from this bond					
		//if (iter==1) // assuming v0 == 0.0
		//{
		//	// Calculate node velocity at the first half time step
		//	for (int k=0;k!=3;++k)
		//	{
		//		vBuffer[k] = dt*(pdForceSum[k]+bodyForce[k])/(0.5*denst);
		//	}
		//}
		//else
		//{
		//	// Get node velocity at the previous half time step
		//	double vPrev[3] = {node_i->GetV1(), node_i->GetV2(), node_i->GetV3()};
		//	// Calculate node velocity at current half time step
		//	for (int k=0;k!=3;++k)
		//	{
		//		vBuffer[k] = ((2.0*denst-_dampCoeff*dt)*vPrev[k] + 2.0*dt*(pdForceSum[k]+bodyForce[k]))
		//			/(2.0*denst + _dampCoeff*dt);
		//	}
		//}
		//node_i->AddV1Buffer(vBuffer[0]);
		//node_i->AddV2Buffer(vBuffer[1]);
		//node_i->AddV3Buffer(vBuffer[2]);
		//pdForceSum[0] = pdForceSum[1] = pdForceSum[2] = 0.0;

		// apply linear viscous damping
		node_i->AddV1Buffer(-_dampCoeff * node_i->GetV1() * dt / denst);
		node_i->AddV2Buffer(-_dampCoeff * node_i->GetV2() * dt / denst);
		node_i->AddV3Buffer(-_dampCoeff * node_i->GetV3() * dt / denst);		

		//// Update node velocity buffer (Adaptive dynamic relaxation)
		//// Get the row number of node i in three dof
		//int irowx = _rowOfNode[node_i->GetID()];
		//int irowy = irowx + 1;
		//int irowz = irowx + 2;
		//// Get the elements in the inverse mass matrix corresponding to this node
		//double invMass[3] = {1.0/mMatrix->GetCoeff(irowx), 1.0/mMatrix->GetCoeff(irowy),
		//	1.0/mMatrix->GetCoeff(irowz)};
		//if (iter==1) // assuming v0 == 0.0
		//{
		//	// Calculate node velocity at the first half time step
		//	for (int k=0;k!=3;++k)
		//	{
		//		vBuffer[k] = dt*invMass[k]*pdForceSum[k]*0.5;
		//	}
		//}
		//else
		//{
		//	// Get node velocity at the previous half time step
		//	double vPrev[3] = {node_i->GetV1(), node_i->GetV2(), node_i->GetV3()};
		//	// Calculate node velocity at current half time step
		//	for (int k=0;k!=3;++k)
		//	{
		//		vBuffer[k] = ((2.0-dampCoeff*dt)*vPrev[k] + 2.0*dt*invMass[k]*pdForceSum[k])
		//			/(2.0+dampCoeff*dt);
		//	}
		//}
		//node_i->AddV1Buffer(vBuffer[0]);
		//node_i->AddV2Buffer(vBuffer[1]);
		//node_i->AddV3Buffer(vBuffer[2]);
		///*if (iter>1)
		//{
		//	fout << node_i->GetID() << "\t"
		//		<< node_i->GetX1() << "\t" << node_i->GetX2() << "\t" << node_i->GetX3() << "\t"
		//		<< pdForceSum[0] << "\t" << pdForceSum[1] << "\t" << pdForceSum[2] << "\t"
		//		<< invMass[0] << "\t" << invMass[1] << "\t" << invMass[2] << "\t"
		//		<< vBuffer[0] << "\t" << vBuffer[1] << "\t" << vBuffer[2] << "\n";
		//}*/
		//pdForceSum[0] = pdForceSum[1] = pdForceSum[2] = 0.0;
	}
	/*** loop over nodes ends ***/

	/*if (iter>1)
	{
		fout << iter << "\t" << dampCoeff << endl;
	}
	delete mMatrix;*/
}

void pdDataManager::CalcBondBasedMotionRK4(const double dt, const double time, const vector<double>& displStep, 
										   vector<double>& k)
{
	/*

	This function updates the node velocity using trial displacements 'displStep' from RK4 steps. 
	The nodal displacements are updated after the force calculation is done.

	Input:
	  1. dt: the current time step.
	  2. time: the current run time till current step.
	  3. displStep: the trial nodal displacements. 
	
	Outupt:
	  1. k: the four K in the RK4 steps.

	*/

	double tol = _tolerance;
	double pdForce[3]; // store peridynamic force between node i and j
	double bodyForce[3]; // store body force on node i
	int numIntg = 0;

	/*string ofpath = _workDir + "pds_debug.dat";
	ofstream fout(ofpath.c_str(), ios::app);*/

	// first update nodal velocity based on the peridynamic force integration over volume 
	/*** loop over nodes starts ***/
	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{
		pdNode* node_i = *itor;
		// Get the row number of node i in three dof
		int irowx = _rowOfNode[node_i->GetID()];
		int irowy = irowx + 1;
		int irowz = irowx + 2;
		pdMaterial* mat_i = node_i->GetMaterial();
		double denst = mat_i->GetDensity();
		// get body force
		bodyForce[0] = mat_i->GetGrav1()*denst;
		bodyForce[1] = mat_i->GetGrav2()*denst;
		bodyForce[2] = mat_i->GetGrav3()*denst;
		// loop over family of node i
		int numFam = node_i->GetNumFamily();
		for (int j=0;j!=numFam;++j)
		{
			pdNode* node_j = node_i->GetFamilyNode(j);
			pdMaterial* mat_j = node_j->GetMaterial();
			pdBond* bond_ij = node_i->GetFamilyBond(j);
			if (!bond_ij->IsBreak())  // check if this bond is already broken
			{
				// force calculation
				if (_cciOn)
				{
					/*** CCI method ***/
					double xi1 = node_j->GetX1() - node_i->GetX1();
					double xi2 = node_j->GetX2() - node_i->GetX2();
					double xi3 = node_j->GetX3() - node_i->GetX3();
					double xi[3] = {xi1, xi2, xi3};
					// Get the row number of node j in three dof
					int jrowx = _rowOfNode[node_j->GetID()];
					int jrowy = jrowx + 1;
					int jrowz = jrowx + 2;
					// calculate eta
					double u1i, u2i, u3i, u1j, u2j, u3j;
					if (node_i->GetBC())
					{
						u1i = node_i->GetU1();
						u2i = node_i->GetU2();
						u3i = node_i->GetU3();
					}
					else
					{
						u1i = displStep[irowx];
						u2i = displStep[irowy];
						u3i = displStep[irowz];
					}
					if (node_j->GetBC())
					{
						u1j = node_j->GetU1();
						u2j = node_j->GetU2();
						u3j = node_j->GetU3();
					}
					else
					{
						u1j = displStep[jrowx];
						u2j = displStep[jrowy];
						u3j = displStep[jrowz];
					}
					double eta1 = u1j - u1i;
					double eta2 = u2j - u2i;
					double eta3 = u3j - u3i;
					double eta[3] = {eta1, eta2, eta3};
					// get the average critical strecth, spring constant and boundary effect compensation factor
					double ecrit = min(mat_i->GetCriticalStretch(), mat_j->GetCriticalStretch());
					double spring = min(mat_i->GetSpringConstant(), mat_j->GetSpringConstant());
					double bdryEffij = 0.5*(node_i->GetBdryEffFactor() + node_j->GetBdryEffFactor());
					double elast;
					CalcBondForceCCI(xi, eta, ecrit, spring, bdryEffij, 2.0, pdForce, elast);
					node_i->AddElastEnergyDensity(elast);
				}
				else
				{
					/*** Adaptive trapezoidal integration or Gaussian quadrature integration ***/
					// Note: the displacement of node_j is not used in the integration.
					// But the displacesment of node_i must be updated according to the trial nodal displacements in RK procedure.
					double u1i, u2i, u3i;
					if (node_i->GetBC())
					{
						u1i = node_i->GetU1();
						u2i = node_i->GetU2();
						u3i = node_i->GetU3();
					}
					else
					{
						u1i = displStep[irowx];
						u2i = displStep[irowy];
						u3i = displStep[irowz];
					}
					double u[3] = { u1i, u2i, u3i };
					// (1) Integration scheme 
					// works for both AI and fixed Gaussian integration
					if (_aiOn)
					{
						// get the trial displacements of contributing nodes for family node_j
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
							// Get the row number of each contributing node in three dof
							int jrowx = _rowOfNode[cbList[j]->GetID()];
							int jrowy = jrowx + 1;
							int jrowz = jrowx + 2;
							double u1j = displStep[jrowx];
							double u2j = displStep[jrowy];
							double u3j = displStep[jrowz];
							cbDispl[j][0] = u1j;
							cbDispl[j][1] = u2j;
							cbDispl[j][2] = u3j;
						}
						CalcBondForceIntg(bond_ij, u, node_i, cbDispl, node_j, time, pdForce, numIntg);
					}
					else // fixed Gaussian integration
					{
						// (2) Summation scheme 
						// only works for fixed Gaussian integration since the integration limits are always the same
						// but it is much faster than the integration scheme
						double x1i = node_i->GetX1();
						double x2i = node_i->GetX2();
						double x3i = node_i->GetX3();					
						double ecrit = min(mat_i->GetCriticalStretch(), mat_j->GetCriticalStretch());
						double spring = min(mat_i->GetSpringConstant(), mat_j->GetSpringConstant());
						double bdryEffij = 0.5*(node_i->GetBdryEffFactor() + node_j->GetBdryEffFactor());
						double elast;
						double pdForcePoint[3]; // force from each Gaussian point
						pdForce[0] = pdForce[1] = pdForce[2] = 0.0; // zero the summation 
						int numGauss = node_j->GetGaussianPointListSize();
						for (int i=0;i!=numGauss;++i)
						{
							pdGauss* gPoint = node_j->GetGaussianPoint(i);
							// get relative position xi
							double xi1 = gPoint->GetX1() - x1i;
							double xi2 = gPoint->GetX2() - x2i;
							double xi3 = gPoint->GetX3() - x3i;
							double xi[3] = {xi1, xi2, xi3};
							// get relative displacement eta	
							double eta1 = gPoint->GetU1() - u1i;
							double eta2 = gPoint->GetU2() - u2i;
							double eta3 = gPoint->GetU3() - u3i;
							double eta[3] = {eta1, eta2, eta3};
							// get the effective weight
							double w1 = gPoint->GetWeightX1(); // weights of Gaussian point in three directions
							double w2 = gPoint->GetWeightX2();
							double w3 = gPoint->GetWeightX3();
							double effWeight = pow(w1 * w2 * w3, 1.0 / 3.0); // take the geometric average
							CalcBondForceCCI(xi, eta, ecrit, spring, bdryEffij, effWeight, pdForcePoint, elast);
							node_i->AddElastEnergyDensity(elast);
							// Summation of the force between node_i and the Gaussian points in node_j
							pdForce[0] += (pow(_gridSpacing, 3)/8.0)*w1*w2*w3*pdForcePoint[0];
							pdForce[1] += (pow(_gridSpacing, 3)/8.0)*w1*w2*w3*pdForcePoint[1];
							pdForce[2] += (pow(_gridSpacing, 3)/8.0)*w1*w2*w3*pdForcePoint[2];
						}
					}
				}
				// update node velocity buffer
				node_i->AddV1Buffer((pdForce[0] + bodyForce[0]) * dt / denst);
				node_i->AddV2Buffer((pdForce[1] + bodyForce[1]) * dt / denst);
				node_i->AddV3Buffer((pdForce[2] + bodyForce[2]) * dt / denst);
			}
			else
			{
				pdForce[0] = 0.0;
				pdForce[1] = 0.0;
				pdForce[2] = 0.0;
			}	
		}
	}
	/*** loop over nodes ends ***/
	
	// next apply boundary condition to affected nodes
	UpdateNodeStatusWithBC(time);
	/*** loop over nodes starts ***/
	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{
		pdNode* node_i = *itor;
		double vol = node_i->GetNodeVolume();
		double denst = node_i->GetMaterial()->GetDensity();
		// updated node velocity buffer at current time step
		double v1 = node_i->GetV1Buffer();
		double v2 = node_i->GetV2Buffer();
		double v3 = node_i->GetV3Buffer();	
		// store the result of k
		k.push_back(v1*dt);
		k.push_back(v2*dt);
		k.push_back(v3*dt);
	}
	/*** loop over nodes ends ***/
}

void pdDataManager::RK4(const double dt, const double time)
{
	/*

	This function implements an explicit RK4 time integration scheme for the dynamic solution.
	For each iteration:
	  1. For the first step, use the updated displacements of integration points (from pervious iteration) to calculate the force.
	  2. For the 2nd, 3rd and 4th step, first update displacements of integration points in CalcIntgPointDisplRK4(.) function. Then call
	     CalcBondBasedMotionRK4(.) to calculate the force. After that, the nodal displacements are updated in CalcBondBasedMotionRK4(.).
	  3. After the 4th step, use all the K to accumulate displacement increments with proper weights. Then update displacements of integration points.
	     Also update node energy and total energy.
      4. Note: the nodes/integration points with or without boundary condtions are treated differently.

   Input:
     1. dt: current time step.
	 2. time: current run time.

	*/

	vector<double> u1, u234, k1, k24, k3;
	double dtHalf = 0.5*dt;
	double dt6 = 1.0/6.0;

	// fill the starting displacements with displacements of nodes from last step
	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{
		pdNode* node_i = *itor;
		u1.push_back(node_i->GetU1());
		u1.push_back(node_i->GetU2());
		u1.push_back(node_i->GetU3());
	}

	/*string ofpath = _workDir + "pds_debug.dat";
	ofstream fout(ofpath.c_str(), ios::app);*/

	// step 1
	//fout << endl<< "***********STEP 1***************" << endl;
	CalcBondBasedMotionRK4(dtHalf, time, u1, k1);
	// step 2
	//fout << endl<< "***********STEP 2***************" << endl;
	for (int i=0;i!=_numNodes*3;++i)
	{
		u234.push_back(u1[i] + 0.5*k1[i]);
	}
	CalcIntgPointDisplRK4(time, u234);
	CalcBondBasedMotionRK4(dtHalf, time, u234, k24);
	// step 3
	//fout << endl<< "***********STEP 3***************" << endl;
	for (int i=0;i!=_numNodes*3;++i)
	{
		u234[i] = u1[i] + 0.5*k24[i];
	}
	CalcIntgPointDisplRK4(time, u234);
	CalcBondBasedMotionRK4(dt, time, u234, k3);
	// add results of step 2 to step 3, clear step 2 container for step 4
	//fout << endl<< "***********STEP 4***************" << endl;
	for (int i=0;i!=_numNodes*3;++i)
	{
		u234[i] = u1[i] + k3[i];
		k3[i] += k24[i];
	}
	k24.clear();
	CalcIntgPointDisplRK4(time, u234);
	CalcBondBasedMotionRK4(dt, time, u234, k24);
	// accumulate increments with proper weights
	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{
		pdNode* node_i = *itor;
		double denst = node_i->GetMaterial()->GetDensity();
		// Get the row number of node i in three dof
		int irowx = _rowOfNode[node_i->GetID()];
		int irowy = irowx + 1;
		int irowz = irowx + 2;
		// node displacement at previous time step
		double u1 = node_i->GetU1();
		double u2 = node_i->GetU2();
		double u3 = node_i->GetU3();		
		// store nodal displacements
		if (node_i->GetBC())
		{
			node_i->SetU1(u1 + dt6*(k1[irowx] + k24[irowx] + 2.0*k3[irowx]));
			node_i->SetU2(u2 + dt6*(k1[irowy] + k24[irowy] + 2.0*k3[irowy]));
			node_i->SetU3(u3 + dt6*(k1[irowz] + k24[irowz] + 2.0*k3[irowz]));
		}
		else
		{
			// apply linear viscous damping on nodal displacement
			node_i->SetU1(u1 + dt6*(k1[irowx] + k24[irowx] + 2.0*k3[irowx]) - dt*_dampCoeff*node_i->GetV1()*dt/denst);
			node_i->SetU2(u2 + dt6*(k1[irowy] + k24[irowy] + 2.0*k3[irowy]) - dt*_dampCoeff*node_i->GetV2()*dt/denst);
			node_i->SetU3(u3 + dt6*(k1[irowz] + k24[irowz] + 2.0*k3[irowz]) - dt*_dampCoeff*node_i->GetV3()*dt/denst);
		}
		// apply linear viscous damping on nodal velocity
		node_i->AddV1Buffer(-_dampCoeff * node_i->GetV1() * dt / denst);
		node_i->AddV2Buffer(-_dampCoeff * node_i->GetV2() * dt / denst);
		node_i->AddV3Buffer(-_dampCoeff * node_i->GetV3() * dt / denst);
		// store nodal velocity
		node_i->SetV1(node_i->GetV1Buffer());
		node_i->SetV2(node_i->GetV2Buffer());
		node_i->SetV3(node_i->GetV3Buffer());	
	}

	// update node energy and total energy
	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{
		pdNode* node_i = *itor;
		// find total strain energy, dissipated energy, momentum and kinetic energy
		double v1i = node_i->GetV1();
		double v2i = node_i->GetV2();
		double v3i = node_i->GetV3();
		double velsq = v1i*v1i + v2i*v2i + v3i*v3i;
		double vol_i = node_i->GetNodeVolume();
		double denst_i = node_i->GetMaterial()->GetDensity();
		// this is the total elastic energy not elastic energy density
		_elasticEnergyStep += vol_i*node_i->GetElastEnergyDensity();
		_kineticEnergyStep += 0.5*denst_i*velsq*vol_i;
		// find kinetic energy of each node.
		node_i->SetKinetEnergy(0.5*denst_i*vol_i*velsq);
	}
	// end of loop over nodes
	// update energy
	_elasticEnergyTotal += _elasticEnergyStep;
	_externalWorkTotal += _externalWorkStep;

	// finally update approximated displacments at the integration points
	u1.clear();
	// fill the updated nodal displacements (current step)
	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{
		pdNode* node_i = *itor;
		u1.push_back(node_i->GetU1());
		u1.push_back(node_i->GetU2());
		u1.push_back(node_i->GetU3());
	}
	CalcIntgPointDisplRK4(time, u1);

	u1.clear();
	u234.clear();
	k1.clear();
	k24.clear();
	k3.clear();
}

void pdDataManager::UpdateNodeStatus(const double dt, const double time)
{
	/* 

	This function updates status for each node:
	  1. apply boundary condition to affected nodes, update node displacement and velocity.
	  2. update nodes energy and total energy.
	  3. calculate approximated displacments at Gaussian points.

	Input: 
	  1. dt: the current time step.
	  2. time: the current total run time.

    */

	//// MPI
	//if (_numProcs>1)
	//{
	//	// communication
	//	// need sychrinize before receive
	//	MPI::Request recv_req_i, recv_req_j;
	//	MPI::Status stat;
	//	double vadd_i[2], vadd_j[2];
	//	vector<pdBond*> listed;
	//	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	//	{
	//		pdNode* node_i = *itor;
	//		for (int n=0;n!=node_i->GetNumFamily();++n)
	//		{
	//			pdBond* bond_ij = node_i->GetFamilyBond(n);
	//			vector<pdBond*>::iterator itor_b = find(_myBonds.begin(), _myBonds.end(), bond_ij);
	//			vector<pdBond*>::iterator itor_l = find(listed.begin(), listed.end(), bond_ij);
	//			if (itor_b==_myBonds.end() && itor_l==listed.end())
	//			{
	//				recv_req_i = MPI::COMM_WORLD.Irecv(&vadd_i, 3, MPI::DOUBLE, MPI::ANY_SOURCE, bond_ij->GetNodeI()->GetID());
	//				recv_req_j = MPI::COMM_WORLD.Irecv(&vadd_j, 3, MPI::DOUBLE, MPI::ANY_SOURCE, bond_ij->GetNodeJ()->GetID());
	//				recv_req_i.Wait(stat);
	//				recv_req_j.Wait(stat);
	//				// update node velocity complete
	//				bond_ij->GetNodeI()->AddV1Buffer(vadd_i[0]);
	//				bond_ij->GetNodeI()->AddV2Buffer(vadd_i[1]);
	//				bond_ij->GetNodeI()->AddV3Buffer(vadd_i[2]);
	//				bond_ij->GetNodeJ()->AddV1Buffer(vadd_j[0]);
	//				bond_ij->GetNodeJ()->AddV2Buffer(vadd_j[1]);
	//				bond_ij->GetNodeJ()->AddV3Buffer(vadd_j[2]);
	//				// save ptr to this bond to listed container to prevent duplicate update
	//				listed.push_back(bond_ij);
	//			}
	//		}
	//	}
	//}

	// apply boundary condition to affected nodes
	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{
		pdNode* node_i = *itor;
		double vol = node_i->GetNodeVolume();
		double denst = node_i->GetMaterial()->GetDensity();
		// node velocity at current time step
		double v1 = node_i->GetV1Buffer();
		double v2 = node_i->GetV2Buffer();
		double v3 = node_i->GetV3Buffer();
		double f_int_1 = denst*(v1-node_i->GetV1())/dt; // for calculation of external work done on boundary
		double f_int_2 = denst*(v2-node_i->GetV2())/dt;
		double f_int_3 = denst*(v3-node_i->GetV3())/dt;
		// node displacement at previous time step
		double u1 = node_i->GetU1();
		double u2 = node_i->GetU2();
		double u3 = node_i->GetU3();
		// apply boundary conditions. If no BC applied, do nothing
		if (node_i->GetBC()!=0)
		{
			pdBdryCondition* bc = node_i->GetBC();
			// find boundary condition type, direction flag and time-dependant factor
			pdBdryCondition::BCType bc_type = bc->GetType();
			int dir_x1 = bc->GetDirFlagX1();
			int dir_x2 = bc->GetDirFlagX2();
			int dir_x3 = bc->GetDirFlagX3();
			double fac = bc->GetTimeFactor(time);
			// type 1: prescribed displacement boundary conditon
			if (bc_type==DISPL)
			{
				if (dir_x1==1)
				{
					u1 = bc->GetValue(1)*fac;
					v1 = 0.0;
					// also update the displacement of Gaussian points inside the cell
					int numGauss = node_i->GetGaussianPointListSize();
					for (int g=0;g!=numGauss;g++)
						node_i->GetGaussianPoint(g)->SetU1(u1);
				}
				if (dir_x2==1)
				{
					u2 = bc->GetValue(2)*fac;
					v2 = 0.0;
					int numGauss = node_i->GetGaussianPointListSize();
					for (int g=0;g!=numGauss;g++)
						node_i->GetGaussianPoint(g)->SetU2(u2);
				}
				if (dir_x3==1)
				{
					u3 = bc->GetValue(3)*fac;
					v3 = 0.0;
					int numGauss = node_i->GetGaussianPointListSize();
					for (int g=0;g!=numGauss;g++)
						node_i->GetGaussianPoint(g)->SetU3(u3);
				}
			}
			// type 2: prescribed velocity boundary condition
			else if (bc_type==VELOC)
			{
				if (dir_x1==1)
				{
					v1 = bc->GetValue(1)*fac;
					_externalWorkStep -= vol*f_int_1*v1*dt;
				}
				if (dir_x2==1)
				{
					v2 = bc->GetValue(2)*fac;
					_externalWorkStep -= vol*f_int_2*v2*dt;
				}
				if (dir_x3==1)
				{
					v3 = bc->GetValue(3)*fac;
					_externalWorkStep -= vol*f_int_3*v3*dt;
				}
			}				
			// type 3: prescribed displacement gradient (12/06/09)
			else if (bc_type == DISPLGRAD)
			{
				if (dir_x1==1)
				{
					u1 = (bc->GetValue(1, 1)*node_i->GetX1()
						+ bc->GetValue(1, 2)*node_i->GetX2()
						+ bc->GetValue(1, 3)*node_i->GetX3())*fac;
					v1 = 0.0;
					// also update the displacement of Gaussian points inside the cell
					int numGauss = node_i->GetGaussianPointListSize();
					for (int g=0;g!=numGauss;g++)
					{
						pdGauss* gPoint = node_i->GetGaussianPoint(g);
						double x1g = gPoint->GetX1();
						double x2g = gPoint->GetX2();
						double x3g = gPoint->GetX3();
						double u1g = (bc->GetValue(1, 1)*x1g + bc->GetValue(1, 2)*x2g + bc->GetValue(1, 3)*x3g)*fac;
						gPoint->SetU1(u1g);
					}
				}
				if (dir_x2==1)
				{
					u2 = (bc->GetValue(2, 1)*node_i->GetX1()
						+ bc->GetValue(2, 2)*node_i->GetX2()
						+ bc->GetValue(2, 3)*node_i->GetX3())*fac;
					v2 = 0.0;
					int numGauss = node_i->GetGaussianPointListSize();
					for (int g=0;g!=numGauss;g++)
					{
						pdGauss* gPoint = node_i->GetGaussianPoint(g);
						double x1g = gPoint->GetX1();
						double x2g = gPoint->GetX2();
						double x3g = gPoint->GetX3();
						double u2g = (bc->GetValue(2, 1)*x1g + bc->GetValue(2, 2)*x2g + bc->GetValue(2, 3)*x3g)*fac;
						gPoint->SetU2(u2g);
					}
				}
				if (dir_x3==1)
				{
					u3 = (bc->GetValue(3, 1)*node_i->GetX1()
						+ bc->GetValue(3, 2)*node_i->GetX2()
						+ bc->GetValue(3, 3)*node_i->GetX3())*fac;
					v3 = 0.0;
					int numGauss = node_i->GetGaussianPointListSize();
					for (int g=0;g!=numGauss;g++)
					{
						pdGauss* gPoint = node_i->GetGaussianPoint(g);
						double x1g = gPoint->GetX1();
						double x2g = gPoint->GetX2();
						double x3g = gPoint->GetX3();
						double u3g = (bc->GetValue(3, 1)*x1g + bc->GetValue(3, 2)*x2g + bc->GetValue(3, 3)*x3g)*fac;
						gPoint->SetU3(u3g);
					}
				}
			}
		}
		// get the node displacement for the current time step
		u1 += v1*dt;
		u2 += v2*dt;
		u3 += v3*dt;
		// store the current node displacement and velocity
		node_i->SetU1(u1);
		node_i->SetU2(u2);
		node_i->SetU3(u3);
		node_i->SetV1(v1);
		node_i->SetV2(v2);
		node_i->SetV3(v3);		
	}
	// end of loop over nodes

	// update node energy and total energy
	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{
		pdNode* node_i = *itor;
		// find total strain energy, dissipated energy, momentum and kinetic energy
		double v1i = node_i->GetV1();
		double v2i = node_i->GetV2();
		double v3i = node_i->GetV3();
		double velsq = v1i*v1i + v2i*v2i + v3i*v3i;
		double vol_i = node_i->GetNodeVolume();
		double denst_i = node_i->GetMaterial()->GetDensity();
		// this is the total elastic energy not elastic energy density
		_elasticEnergyStep += vol_i*node_i->GetElastEnergyDensity();
		_kineticEnergyStep += 0.5*denst_i*velsq*vol_i;
		// find kinetic energy of each node.
		node_i->SetKinetEnergy(0.5*denst_i*vol_i*velsq);
	}
	// end of loop over nodes
	// update energy
	_elasticEnergyTotal += _elasticEnergyStep;
	_externalWorkTotal += _externalWorkStep;

	// update approximated displacments at the Gaussian points
	if (_fgiOn)
	{
		if (_mlsOn) // initial Gaussian points are greater than one
		{		
			for (vector<pdNode*>::iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
			{
				pdNode* node = *itor;
				CalcApproxDispl(node);
			}
		}
		else // one Gaussian point, so mls is turned off
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

	/*// for debug
	string ofpath = _workDir + "pds_debug.dat";
	ofstream fout(ofpath.c_str(), ios::app);
	fout << "current time = " << time << endl;
	for (vector<pdNode*>::const_iterator itor_i=_nodes.begin();itor_i!=_nodes.end();++itor_i)
	{
		pdNode* node = *itor_i;
		fout << "node " << node->GetID() 
			<< "\t" << node->GetX1() << "\t" << node->GetX2() << "\t" << node->GetX3()
			<< "\t" << node->GetU1() << "\t" << node->GetU2() << "\t" << node->GetU3();
		int gPointNum = node->GetGaussianPointListSize();
		for (int g=0;g!=gPointNum;++g)
		{
			pdGauss* gPoint = node->GetGaussianPoint(g);
			fout << "\t" << g << "\t" << gPoint->GetU1() << "\t" << gPoint->GetU2() << "\t" << gPoint->GetU3();
		}
		fout << endl;
	}
	fout << endl << endl;*/
}


void pdDataManager::CalcIntgPointDisplRK4(const double time, const vector<double>& displStep)
{
	/* 

	This function updates the displacements of the integration points with the trial displacemens from RK4 steps.
	Integration points with or without boundary conditions are treated differently.
	  1. adaptive integration method
	  2. fixed Gaussian integration method

	Input:
	  1. time: current run time till this time step.
	  2. displStep: trial nodal displacements.

    */

	// first update Gaussian points with their nodes associated with boundary conditions.
	if (_fgiOn)
	{
		UpdateIntgPointStatusWithBC(time);
	}
	else if (_aiOn)
	{
		// update is done in pdIntgManager::CalcBondBasedForce(.) function
	}
	else
	{
	}

	// next update rest of the Gaussian points (no Boundary condition).
	if (_fgiOn)
	{
		if (_mlsOn) // initial Gaussian points are greater than one
		{		
			for (vector<pdNode*>::iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
			{
				pdNode* node = *itor;
				if (node->GetBC()==0)
				{
					int cbNumber = node->GetContribNodeListSize();
					int gPointNum = node->GetGaussianPointListSize();
					for (int g=0;g!=gPointNum;++g)
					{
						pdGauss* gPoint = node->GetGaussianPoint(g);
						pdVector* coeff = gPoint->GetCoeff();
						// get the approximated displacement at this Gaussian point
						double appDisplX1 = 0.0;
						double appDisplX2 = 0.0;
						double appDisplX3 = 0.0;
						for (int i = 0; i != cbNumber; ++i)
						{
							// Get the row number of node i in three dof
							int irowx = _rowOfNode[node->GetContribNode(i)->GetID()];
							int irowy = irowx + 1;
							int irowz = irowx + 2;
							// instead of using nodal displacement, use the trial displacements
							appDisplX1 += displStep[irowx] * coeff->GetCoeff(i);
							appDisplX2 += displStep[irowy] * coeff->GetCoeff(i);
							appDisplX3 += displStep[irowz] * coeff->GetCoeff(i);
						}
						// save the approximate displacment
						gPoint->SetU1(appDisplX1);
						gPoint->SetU2(appDisplX2);
						gPoint->SetU3(appDisplX3);					
					}
				}
			}
		}
		else // one Gaussian point, so mls is turned off
		{
			// copy the nodal displacement to its Gaussian point since they are the same
			for (vector<pdNode*>::iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
			{
				if ((*itor)->GetBC()==0)
				{
					// Get the row number of node i in three dof
					int irowx = _rowOfNode[(*itor)->GetID()];
					int irowy = irowx + 1;
					int irowz = irowx + 2;
					(*itor)->GetGaussianPoint(0)->SetU1(displStep[irowx]);
					(*itor)->GetGaussianPoint(0)->SetU2(displStep[irowy]);
					(*itor)->GetGaussianPoint(0)->SetU3(displStep[irowz]);
				}
			}
		}
	}
	else if (_aiOn)
	{
		// update is done in pdIntgManager::CalcBondBasedForce(.) function
	}
	else
	{
		// implement other integration methods
	}
}

void pdDataManager::UpdateNodeStatusWithBC(const double time)
{
	/*

	This function updates nodal displacements and veclocities with associated boundary conditions.
	No node displacement increment is done in this function.
	Nodes with or without boundary conditions are treated differently.

	Input:
	  1. time: current run time till this step.

	*/

	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{
		pdNode* node_i = *itor;
		double vol = node_i->GetNodeVolume();
		double denst = node_i->GetMaterial()->GetDensity();
		// node velocity at current time step
		double v1 = node_i->GetV1Buffer();
		double v2 = node_i->GetV2Buffer();
		double v3 = node_i->GetV3Buffer();
		// node displacement at previous time step
		double u1 = node_i->GetU1();
		double u2 = node_i->GetU2();
		double u3 = node_i->GetU3();
		// apply boundary conditions. If no BC applied, do nothing
		if (node_i->GetBC())
		{
			pdBdryCondition* bc = node_i->GetBC();
			// find boundary condition type, direction flag and time-dependant factor
			pdBdryCondition::BCType bc_type = bc->GetType();
			int dir_x1 = bc->GetDirFlagX1();
			int dir_x2 = bc->GetDirFlagX2();
			int dir_x3 = bc->GetDirFlagX3();
			double fac = bc->GetTimeFactor(time);
			// type 1: prescribed displacement boundary conditon
			if (bc_type==DISPL)
			{
				if (dir_x1==1)
				{
					u1 = bc->GetValue(1)*fac;
					v1 = 0.0;
					}
				if (dir_x2==1)
				{
					u2 = bc->GetValue(2)*fac;
					v2 = 0.0;
				}
				if (dir_x3==1)
				{
					u3 = bc->GetValue(3)*fac;
					v3 = 0.0;
				}
			}
			// type 2: prescribed velocity boundary condition
			else if (bc_type==VELOC)
			{
				if (dir_x1==1)
				{
					v1 = bc->GetValue(1)*fac;
				}
				if (dir_x2==1)
				{
					v2 = bc->GetValue(2)*fac;
				}
				if (dir_x3==1)
				{
					v3 = bc->GetValue(3)*fac;
				}
			}				
			// type 3: prescribed displacement gradient (12/06/09)
			else if (bc_type == DISPLGRAD)
			{
				if (dir_x1==1)
				{
					u1 = (bc->GetValue(1, 1)*node_i->GetX1()
						+ bc->GetValue(1, 2)*node_i->GetX2()
						+ bc->GetValue(1, 3)*node_i->GetX3())*fac;
					v1 = 0.0;
				}
				if (dir_x2==1)
				{
					u2 = (bc->GetValue(2, 1)*node_i->GetX1()
						+ bc->GetValue(2, 2)*node_i->GetX2()
						+ bc->GetValue(2, 3)*node_i->GetX3())*fac;
					v2 = 0.0;
				}
				if (dir_x3==1)
				{
					u3 = (bc->GetValue(3, 1)*node_i->GetX1()
						+ bc->GetValue(3, 2)*node_i->GetX2()
						+ bc->GetValue(3, 3)*node_i->GetX3())*fac;
					v3 = 0.0;
				}
			}
		}
		// update nodal displacement and velocity
		node_i->SetV1Buffer(v1);
		node_i->SetV2Buffer(v2);
		node_i->SetV3Buffer(v3);
		node_i->SetU1(u1);
		node_i->SetU2(u2);
		node_i->SetU3(u3);
	}
}

void pdDataManager::UpdateIntgPointStatusWithBC(const double time)
{
	/*

	This function updates the displacements of integration points using the boundary conditions associated with its nodes.

	Input:
	  1. time: current run time till this step.

	*/

	for (vector<pdNode*>::const_iterator itor=_nodes.begin();itor!=_nodes.end();++itor)
	{
		pdNode* node_i = *itor;
		// apply boundary conditions. If no BC applied, do nothing
		if (node_i->GetBC())
		{
			pdBdryCondition* bc = node_i->GetBC();
			// find boundary condition type, direction flag and time-dependant factor
			pdBdryCondition::BCType bc_type = bc->GetType();
			int dir_x1 = bc->GetDirFlagX1();
			int dir_x2 = bc->GetDirFlagX2();
			int dir_x3 = bc->GetDirFlagX3();
			double fac = bc->GetTimeFactor(time);
			int numGauss = node_i->GetGaussianPointListSize();
			// type 1: prescribed displacement boundary conditon
			if (bc_type==DISPL)
			{
				if (dir_x1==1)
				{					
					for (int g=0;g!=numGauss;g++)
						node_i->GetGaussianPoint(g)->SetU1(bc->GetValue(1)*fac);
				}
				if (dir_x2==1)
				{
					for (int g=0;g!=numGauss;g++)
						node_i->GetGaussianPoint(g)->SetU2(bc->GetValue(2)*fac);
				}
				if (dir_x3==1)
				{
					for (int g=0;g!=numGauss;g++)
						node_i->GetGaussianPoint(g)->SetU3(bc->GetValue(3)*fac);
				}
			}
			// type 3: prescribed displacement gradient (12/06/09)
			else if (bc_type == DISPLGRAD)
			{
				if (dir_x1==1)
				{
					double u1 = (bc->GetValue(1, 1)*node_i->GetX1()
						+ bc->GetValue(1, 2)*node_i->GetX2()
						+ bc->GetValue(1, 3)*node_i->GetX3())*fac;
					for (int g=0;g!=numGauss;g++)
					{
						pdGauss* gPoint = node_i->GetGaussianPoint(g);
						double x1g = gPoint->GetX1();
						double x2g = gPoint->GetX2();
						double x3g = gPoint->GetX3();
						double u1g = (bc->GetValue(1, 1)*x1g + bc->GetValue(1, 2)*x2g + bc->GetValue(1, 3)*x3g)*fac;
						gPoint->SetU1(u1g);
					}
				}
				if (dir_x2==1)
				{
					double u2 = (bc->GetValue(2, 1)*node_i->GetX1()
						+ bc->GetValue(2, 2)*node_i->GetX2()
						+ bc->GetValue(2, 3)*node_i->GetX3())*fac;
					for (int g=0;g!=numGauss;g++)
					{
						pdGauss* gPoint = node_i->GetGaussianPoint(g);
						double x1g = gPoint->GetX1();
						double x2g = gPoint->GetX2();
						double x3g = gPoint->GetX3();
						double u2g = (bc->GetValue(2, 1)*x1g + bc->GetValue(2, 2)*x2g + bc->GetValue(2, 3)*x3g)*fac;
						gPoint->SetU2(u2g);
					}
				}
				if (dir_x3==1)
				{
					double u3 = (bc->GetValue(3, 1)*node_i->GetX1()
						+ bc->GetValue(3, 2)*node_i->GetX2()
						+ bc->GetValue(3, 3)*node_i->GetX3())*fac;
					for (int g=0;g!=numGauss;g++)
					{
						pdGauss* gPoint = node_i->GetGaussianPoint(g);
						double x1g = gPoint->GetX1();
						double x2g = gPoint->GetX2();
						double x3g = gPoint->GetX3();
						double u3g = (bc->GetValue(3, 1)*x1g + bc->GetValue(3, 2)*x2g + bc->GetValue(3, 3)*x3g)*fac;
						gPoint->SetU3(u3g);
					}
				}
			}
			else
			{
				// implement other boundary condition types
			}
		}
	}
}

void pdDataManager::CalcApproxDispl(pdNode* node)
{
	// This function is called every time step to calculate the interpolated displacement 
	// of the Gaussian points for this node using moving least square approximation.
	// the coefficients of each contributing nodes are calculated before the iteration starts
	// since it is only related to the referrence coordinates of the nodes.

	// get the contributing node list of this node 
	vector<pdNode*> cbList;
	node->GetContribNodeList(cbList);
	int ndNumber = node->GetContribNodeListSize();
	int gPointNum = node->GetGaussianPointListSize();
	for (int g=0;g!=gPointNum;++g)
	{
		pdGauss* gPoint = node->GetGaussianPoint(g);
		pdVector* coeff = gPoint->GetCoeff();
		// get the approximated displacement at this Gaussian point
		double appDisplX1 = 0.0;
		double appDisplX2 = 0.0;
		double appDisplX3 = 0.0;		
		for (int i = 0; i != ndNumber; ++i)
		{
			appDisplX1 += cbList[i]->GetU1() * coeff->GetCoeff(i);
			appDisplX2 += cbList[i]->GetU2() * coeff->GetCoeff(i);
			appDisplX3 += cbList[i]->GetU3() * coeff->GetCoeff(i);
		}
		// save the approximate displacment
		if (node->GetBC()==0)
		{
			gPoint->SetU1(appDisplX1);
			gPoint->SetU2(appDisplX2);
			gPoint->SetU3(appDisplX3);
		}
		else
		{
			pdBdryCondition* bc = node->GetBC();
			int dir_x1 = bc->GetDirFlagX1();
			int dir_x2 = bc->GetDirFlagX2();
			int dir_x3 = bc->GetDirFlagX3();
			if (dir_x1==0)
				gPoint->SetU1(appDisplX1);
			if (dir_x2==0)
				gPoint->SetU2(appDisplX2);
			if (dir_x3==0)
				gPoint->SetU3(appDisplX3);
		}
	}
}



