#include "math.h"
#include "pdSolve.h"
#include "pdDataManager.h"
using namespace std;

pdSolve::pdSolve() {}

pdSolve::~pdSolve() {}

void pdSolve::InitialProblem(pdDataManager& dat)
{
	// This function read data from input file and initialize all necessary data for solving process.

	dat.ReadInputFile();
	dat.InitNode();
	dat.InitBond();
	dat.InitMaterial();
	dat.InitBdryEffectFactor();
	dat.InitGridBoundary();	
	dat.InitIntegrationPoint();
//	dat.InitLocalBonds();
}

void pdSolve::SolveDynamic(pdDataManager& dat)
{
	// This function solves problem using dynamic method and writes out data file.

	// open the general info output file
	string outputFileKeyword = dat.GetWkDir() + "pds_output.dat";
	ofstream outputFile(outputFileKeyword.c_str(), ios::app);
	outputFile.setf(ios::scientific);
	outputFile.precision(5);

	// open the step-wise output file
	string stepFileKeyword = dat.GetWkDir() + "pds_steps.dat";
	ofstream stepFile(stepFileKeyword.c_str(), ios::out);
	stepFile << "Data at given node for each time step output mode\n\n";
	stepFile << "Step Time ID X1 X2 X3 U1 U2 U3 V1 V2 V3 MatID BCID bdryEff DmgFrac \n";

	//string ofpath = dat.GetWkDir() + "pds_debug.dat";
	//ofstream fout(ofpath.c_str(), ios::out);

	// Get total time, total step and dump frequency
	// initial current time and time step
	int totalStep = dat.GetTotalStep();
	int lastStep = totalStep;
	int dump = dat.GetStepDumpFrequency();
	double totalTime = dat.GetTotalTime();
	double stableTimeStep = dat.CalcStableTimeStep(); // stable time step
	double currentTime = 0.0; // current time
	if (dat.GetRank()==0)
	{
		outputFile << endl << "** Dynamic solution: starting from cycle 1 **" << endl;
		cout << endl << "** Dynamic solution: starting from cycle 1 **" << endl;
	}
	
	// iteration starts
	double currentTimeStep = stableTimeStep;
	for (int iter=1;iter!=totalStep+1;++iter)
	{		
		currentTime += currentTimeStep;
		if (currentTime <= totalTime)
		{			
			// *** beginning of dynamic solution ***
			dat.InitDynamicSolve();
			/*if (dat.IsStateOn())
			{
				dat.CalcStateBasedMotion(currentTimeStep);
			}
			else
			{
				// bond-based force calculation
			}*/
			//if (_rk4On)
			{
				dat.RK4(currentTimeStep, currentTime);
			}
			//else
			{
				//dat.CalcBondBasedMotion(currentTimeStep, iter);
				//dat.UpdateNodeStatus(currentTimeStep, currentTime);
			}			
			// *** end of dynamic solution ***

			if (dat.GetRank()==0)
			{
				outputFile << iter << ": dt= " << currentTimeStep << " t= " << currentTime 
					<< " KE= " << dat.GetKinetTotal() << " EE= " << dat.GetElastStep() << endl;
				cout << iter << ": dt= " << currentTimeStep << " t= " << currentTime 
					<< " KE= " << dat.GetKinetTotal() << " EE= " << dat.GetElastStep() << endl;
			}
			// output nodal data for first step, last step and each dump frequency
			if (iter==1 || iter==lastStep || iter%dump==0)
			{
				//WriteDetail(dat, iter, currentTime);
				WriteTecplot(dat, iter, currentTime);
				/*if (iter==lastStep)
				{
					string path = dat.GetWkDir() + "pds_plane.dat";
					ofstream ofile(path.c_str(), ios::out);
					dat.PrintForceThroughPlane(ofile, iter);
					ofile.close();
				}*/
			}
			// output nodal data at given node for each time step
			if (1)
			{				
				WriteEachStep(dat, iter, currentTime, stepFile);
			}
			// find the last step if total run time reaches the end
			if (currentTime + stableTimeStep > totalTime && currentTime != totalTime)
			{
				lastStep = iter+1;
				currentTimeStep = totalTime - currentTime;
			}
		}
		else
		{
			break;
		}
	}
	/*if (dat.GetRank()==0)
	{
		outputFile << endl << "Total work done through boundary = " << dat.GetWorkBdryTotal()
			<< endl << "Total elastic energy = " << dat.GetElastTotal()
			<< endl << "Total kinetic energy = " << dat.GetKinetTotal() << endl;
	}*/
	outputFile.close();
	stepFile.close();
}

void pdSolve::SolveStatic(pdDataManager& dat)
{
	// This function solves problem using static method and writes out data file.

	// get total time, total step and dump frequency
	// initial current time and time step
	int total_step = dat.GetTotalStep();
	int last_step = total_step;
	int dump = dat.GetStepDumpFrequency();
	double tend = dat.GetTotalTime();
	double time = 0.;
	string fpath = dat.GetWkDir() + "pds_output.dat";
	ofstream outfile(fpath.c_str(), ios::app);
	if (dat.GetRank()==0)
	{
		outfile << "** Static solution: starting from cycle 1 **" << endl;
		cout << "** Static solution: starting from cycle 1 **" << endl;
	}

	// iteration starts
	for (int iter=1;iter!=total_step+1;++iter)
	{
		// get current time
		double dt = dat.CalcStableTimeStep();
		time += dt;
		if (time<=tend)
		{
			// find the last step if total run time reaches
			if (time+dt>tend)
			{
				last_step = iter;
			}
			// *** beginning of static solution ***
			vector<vector<double>> k_matrix;
			vector<double> loa_vector, displ_vector;
			AssembleEquation(dat, time, k_matrix, loa_vector);
			Gauss(dat, k_matrix, loa_vector, displ_vector);
			// *** end of static solution ***

			// write output file for first step, last step and each dump frequency
			if (iter%dump==0 || iter==1 || iter==last_step)
			{
				WriteDetail(dat, iter, time);
				WriteTecplot(dat, iter, time);
			}

			for (vector<vector<double>>::iterator itor=k_matrix.begin();itor!=k_matrix.end();++itor)
			{
				(*itor).clear();
			}
			k_matrix.clear();
			loa_vector.clear();
			displ_vector.clear();

		}
		else
		{
			break;
		}
	}
	outfile.close();
}

void pdSolve::AssembleEquation(pdDataManager& dat, const double time, vector<vector<double>>& matrix, vector<double>& load)
{
	// This function assemble stiffness matrix (matrix) and load vector (load) for equilibrium equation.
	// Assume |eta|<<1 so the peridynamic function is linearized

	vector<int> row_of_node; // store the row number for 1st dof of each node
	vector<int> bc_row_flag; // store the boundary condition flag for each node (1 for BC and 0 for w/o BC)
	vector<double> bc_row_value; // store the value of boundary condition

	// save the row number, start from 0
	int nrow = 0;
	for (int i=0;i!=dat.GetNumNode();++i)
	{
		row_of_node.push_back(nrow);
		nrow += 3;
	}
	nrow = int(row_of_node.size()); // dimension of matrix (nrow by nrow square)

	// save the boundary condition flag and value for each node
	bc_row_flag.assign(nrow, 0);
	bc_row_value.assign(nrow, 0.);
	for (int i=0;i!=dat.GetNumNode();++i)
	{
		// Get the row number of node i in three dof
		int irowx = row_of_node[i];
		int irowy = irowx + 1;
		int irowz = irowx + 2;
		// Get boundary condition, 0 means no boundary condition is assigned to this node
		if (dat.GetNode(i)->GetBC() != 0)
		{
			pdBdryCondition* bc = dat.GetNode(i)->GetBC();
			pdBdryCondition::BCType bc_type = bc->GetType();
			if (bc_type==pdBdryCondition::DISPL) // prescribed displacement boundary condition
			{
				if (bc->GetDirFlagX1()==1)
				{
					bc_row_flag[irowx] = 1;
					bc_row_value[irowx] = bc->GetValue(1);
				}
				if (bc->GetDirFlagX2()==1)
				{
					bc_row_flag[irowy] = 1;
					bc_row_value[irowy] = bc->GetValue(2);
				}
				if (bc->GetDirFlagX3()==1)
				{
					bc_row_flag[irowz] = 1;
					bc_row_value[irowz] = bc->GetValue(3);
				}
			}
			else if (bc_type==pdBdryCondition::VELOC) // prescribed velocity boundary condition
			{
				if (bc->GetDirFlagX1()==1)
				{
					bc_row_flag[irowx] = 1;
					bc_row_value[irowx] = bc->GetValue(1)*time;
				}
				if (bc->GetDirFlagX2()==1)
				{
					bc_row_flag[irowy] = 1;
					bc_row_value[irowy] = bc->GetValue(2)*time;
				}
				if (bc->GetDirFlagX3()==1)
				{
					bc_row_flag[irowz] = 1;
					bc_row_value[irowz] = bc->GetValue(3)*time;
				}
			}
		}
	}

	// zero k matrix and load vector
	for (int n=0;n!=nrow;++n)
	{
		vector<double> column(nrow, 0.);
		matrix.push_back(column);
	}
	load.assign(nrow, 0.);

	// loop over nodes to set up K matrix
	for (int i=0;i!=dat.GetNumNode();++i)
	{
		pdNode* node_i = dat.GetNode(i);
		// Get the row number of node i in three dof
		int irowx = row_of_node[i];
		int irowy = irowx + 1;
		int irowz = irowx + 2;
		// Get the node coordinate
		double x1i = node_i->GetX1();
		double x2i = node_i->GetX2();
		double x3i = node_i->GetX3();
		pdMaterial* mat_i = node_i->GetMaterial();

		// *** degree of freedom in x1 direction starts ***
		if (bc_row_flag[irowx])
		{
			matrix[irowx][irowx] = 1.;
			load[irowx] = bc_row_value[irowx];
		}
		else
		{
			// loop over its family nodes
			int numfam = node_i->GetNumFamily();
			for (int n=0;n!=numfam;++n)
			{
				pdNode* node_j = node_i->GetFamilyNode(n);
				int j = node_j->GetID();
				// Get the row number of node j in three dof
				int jrowx = row_of_node[j];
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
				double volj = node_j->GetNodeVolume();
				// fill the k matrix and load vector
				matrix[irowx][irowx] += c11*volj;
				if (!bc_row_flag[jrowx])
				{
					matrix[irowx][jrowx] -= c11*volj;
				}
				else
				{
					load[irowx] += c11*volj*bc_row_value[jrowx];
				}
                if (!bc_row_flag[irowy])
				{
					matrix[irowx][irowy] += c12*volj;
				}
				else
				{
					load[irowx] -= c12*volj*bc_row_value[irowy];
				}
				if (!bc_row_flag[jrowy])
				{
					matrix[irowx][jrowy] -= c12*volj;
				}
				else
				{
					load[irowx] += c12*volj*bc_row_value[jrowy];
				}
				if (!bc_row_flag[irowz])
				{
					matrix[irowx][irowz] += c13*volj;
				}
				else
				{
					load[irowx] -= c13*volj*bc_row_value[irowz];
				}
				if (!bc_row_flag[jrowz])
				{
					matrix[irowx][jrowz] -= c13*volj;
				}
				else
				{
					load[irowx] += c13*volj*bc_row_value[jrowz];
				}
			}
		}
		// *** degree of freedom in x1 direction ends ***

		// *** degree of freedom in x2 direction starts ***
		if (bc_row_flag[irowy])
		{
			matrix[irowy][irowy] = 1.;
			load[irowy] = bc_row_value[irowy];
		}
		else
		{
			// loop over its family nodes
			int numfam = node_i->GetNumFamily();
			for (int n=0;n!=numfam;++n)
			{
				pdNode* node_j = node_i->GetFamilyNode(n);
				int j = node_j->GetID();
				// Get the row number of node j in three dof
				int jrowx = row_of_node[j];
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
				double c21 = spring*xi2*xi1/pow(r, 3.0);
				double c22 = spring*xi2*xi2/pow(r, 3.0);
				double c23 = spring*xi2*xi3/pow(r, 3.0);
				double volj = node_j->GetNodeVolume();
				// fill the k matrix
				if (!bc_row_flag[irowx])
				{
					matrix[irowy][irowx] += c21*volj;
				}
				else
				{
					load[irowy] -= c21*volj*bc_row_value[irowx];
				}
                if (!bc_row_flag[jrowx])
				{
					matrix[irowy][jrowx] -= c21*volj;
				}
				else
				{
					load[irowy] += c21*volj*bc_row_value[jrowx];
				}
				matrix[irowy][irowy] += c22*volj;
				if (!bc_row_flag[jrowy])
				{
					matrix[irowy][jrowy] -= c22*volj;
				}
				else
				{
					load[irowy] += c22*volj*bc_row_value[jrowy];
				}
				if (!bc_row_flag[irowz])
				{
					matrix[irowy][irowz] += c23*volj;
				}
				else
				{
					load[irowy] -= c23*volj*bc_row_value[irowz];
				}
				if (!bc_row_flag[jrowz])
				{
					matrix[irowy][jrowz] -= c23*volj;
				}
				else
				{
					load[irowy] += c23*volj*bc_row_value[jrowz];
				}
			}
		}
		// *** degree of freedom in x2 direction ends ***

		// *** degree of freedom in x3 direction starts ***
		if (bc_row_flag[irowz])
		{
			matrix[irowz][irowz] = 1.;
			load[irowz] = bc_row_value[irowz];
		}
		else
		{
			// loop over its family nodes
			int numfam = node_i->GetNumFamily();
			for (int n=0;n!=numfam;++n)
			{
				pdNode* node_j = node_i->GetFamilyNode(n);
				int j = node_j->GetID();
				// Get the row number of node j in three dof
				int jrowx = row_of_node[j];
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
				double c31 = spring*xi3*xi1/pow(r, 3.0);
				double c32 = spring*xi3*xi2/pow(r, 3.0);
				double c33 = spring*xi3*xi3/pow(r, 3.0);
				double volj = node_j->GetNodeVolume();
				// fill the k matrix
				if (!bc_row_flag[irowx])
				{
					matrix[irowz][irowx] += c31*volj;
				}
				else
				{
					load[irowz] -= c31*volj*bc_row_value[irowx];
				}
                if (!bc_row_flag[jrowx])
				{
					matrix[irowz][jrowx] -= c31*volj;
				}
				else
				{
					load[irowz] += c31*volj*bc_row_value[jrowx];
				}
				if (!bc_row_flag[irowy])
				{
					matrix[irowz][irowy] += c32*volj;
				}
				else
				{
					load[irowz] -= c32*volj*bc_row_value[irowy];
				}
				if (!bc_row_flag[jrowy])
				{
					matrix[irowz][jrowy] -= c32*volj;
				}
				else
				{
					load[irowz] += c32*volj*bc_row_value[jrowy];
				}
				matrix[irowz][irowz] += c33*volj;
				if (!bc_row_flag[jrowz])
				{
					matrix[irowz][jrowz] -= c33*volj;
				}
				else
				{
					load[irowz] += c33*volj*bc_row_value[jrowz];
				}
			}
		}
		// *** degree of freedom in x3 direction ends ***
	}

	if (dat.GetRank()==0)
	{
		for(int i=0;i!=nrow;++i)
		{
			for (int j=0;j!=nrow;++j)
			{
				cout.width(12);
				cout.flags(ios::left);
				cout.precision(4);
				cout << matrix[i][j] << "\t";
			}
			cout << endl;
		}
	}
}

void pdSolve::Gauss(pdDataManager& dat, vector<vector<double>> m, vector<double> rhs, vector<double>& sol)
{
	//	Here is the routine for Gauss-Jordan elimination with full pivoting:
	// this function solves equation [matrix]{sol}={rhs} where matrix [m] is n by n
	// the solution is saved in the sol container
	
	// get the dimension of matrix m (which is must be same as the length of rhs vector)
	if (m.size()!=rhs.size())
	{
		if (dat.GetRank()==0)
		{
			cout << "Gauss: the dimension of matrix is different from RHS vector " << endl;
		}
		exit(0);
	}
	else
	{
		int n = int(m.size()); // number of row (same as column)
		int irow, icol;
		//	The integer containers ipiv, indxr, and indxc are used for bookkeeping on the pivoting. 
		vector<int> indxc(n, 0); 
		vector<int> indxr(n, 0);
		vector<int> ipiv(n, 0);
		for (int i=0;i!=n;++i) 
		{ 
			// This is the main loop over the columns to be reduced. 
			double big=0.;
			for (int j=0;j!=n;++j) 
			{
				// This is the outer loop of the search for a pivot element. 
				if (ipiv[j] != 1)
				{	
					for (int k=0;k!=n;++k) 
					{
						if (ipiv[k] == 0) 
						{	
							if (fabs(m[j][k]) >= big) 
							{	
								big = fabs(m[j][k]);
								irow = j;
								icol = k;
							}
						}
						else if (ipiv[k] > 1) 
						{
							if (dat.GetRank()==0)
							{
								cout << "Gauss: Singular Matrix" << endl;
							}
							exit(0);
						}
					}
				}
			}
			++(ipiv[icol]);

		//	We now have the pivot element, so we interchange rows, if needed, to put the pivot
		//	element on the diagonal. The columns are not physically interchanged, only relabeled:
		//	indxc[i], the column of the ith pivot element, is the ith column that is reduced, while
		//	indxr[i] is the row in which that pivot element was originally located. If indxr[i]
		//	6 = indxc[i] there is an implied column interchange. With this form of bookkeeping, the
		//	solution b's will end up in the correct order, and the inverse matrix will be scrambled
		//	by columns. 
			if (irow != icol)
			{
				for (int l=0;l!=n;++l)
				{
					swap(m[irow][l], m[icol][l]);
				}
				double hold = rhs[irow];
				rhs[irow] = rhs[icol];
				rhs[icol] = hold;
			}

		// We are now ready to divide the pivot row by the pivot element, located at irow and icol.	
			indxr[i]=irow; 
			indxc[i]=icol;
			if (m[icol][icol] == 0.0) 
			{
				if (dat.GetRank()==0)
				{
					cout << "gaussj: Singular Matrix-2" << endl;
				}
				exit(0);
			}
			double pivinv = 1.0/m[icol][icol];
			m[icol][icol]=1.0;
			for ( int l=0;l!=n;++l) 
			{
				m[icol][l] *= pivinv;
			}
			rhs[icol] *= pivinv;
			for (int ll=0;ll!=n;++ll) //Next, we reduce the rows... except for the pivot one, of course.
			{
				if (ll != icol) 
				{	
					double dum = m[ll][icol];
					m[ll][icol]=0.0;
					for (int l=0;l!=n;++l) 
					{
						m[ll][l] -= m[icol][l]*dum;
					}
					rhs[ll] = -dum*rhs[icol];
				}
			}
		}

	//	This is the end of the main loop over columns of the reduction. It only remains to unscram-
	//	ble the solution in view of the column interchanges. We do this by interchanging pairs of
	//	columns in the reverse order that the permutation was built up.
		for (int l=n-1;l!=-1;--l)
		{	
			if (indxr[l] != indxc[l])
			{	
				for (int k=0;k!=n;++k)
				{	
					swap(m[k][indxr[l]],m[k][indxc[l]]);
				}
			}
		} // And we are done.
		indxc.clear();
		indxr.clear();
		ipiv.clear();

		for(int i=0; i!=n; ++i){	//copy the solution to sol 
			sol[i] = rhs[i];
		}
	}
}

void pdSolve::WriteDetail(const pdDataManager& dat, const int step, double const time)
{
	// This function outputs nodal data for each dump frequency, first and last time step.

	if (dat.GetRank()==0)
	{
		string keywd = dat.GetWkDir() + "pds_output_";
		keywd = dat.KeywordAppend(keywd, step);
		keywd += ".dat";
		ofstream outfile(keywd.c_str(), ios::out);
		outfile << "Detail data output mode\n";
		outfile.setf(ios::scientific);
		outfile.precision(5);
		// problem Setup data
		outfile << "Problem: " << dat.GetProblemTitle() << endl
			<< "Total nodes = " << dat.GetNumNode() << endl
			<< "Current step = " << step << endl
			<< "Current time = " << time << endl
			<< "Total work done through boundary = " << dat.GetWorkBdryStep() << endl
			<< "Total elastic energy = " << dat.GetElastStep() << endl
			<< "Total kinetic energy = " << dat.GetKinetTotal() << endl
			<< "Total force through plane x1 = " << dat.GetX1force() << " is " << dat.GetForce1() << endl
			<< "Total force through plane x2 = " << dat.GetX2force() << " is " << dat.GetForce2() << endl
			<< "Total force through plane x3 = " << dat.GetX3force() << " is " << dat.GetForce3() << endl;
		// nodal data
		for (int i=0;i!=dat.GetNumNode();++i)
		{
			outfile << *(dat.GetNode(i)) << endl;
		}
	outfile.close();
	}
}

void pdSolve::WriteTecplot(const pdDataManager& dat, const int step, const double time)
{
	// This function outputs nodal data with Tecplot formatting.

	if (dat.GetRank()==0)
	{
		string keywd = dat.GetWkDir() + "pds_tecplot_";
		keywd = dat.KeywordAppend(keywd, step);
		keywd += ".dat";
		ofstream outfile(keywd.c_str(), ios::out);
		outfile << "Tecplot formatting output mode\n";
		outfile.setf(ios::scientific);
		outfile.precision(5);
		// write header
		outfile << "VARIABLES = \n ID X1 X2 X3 U1 U2 U3 V1 V2 V3 MatID BCID bdryEff DmgFrac \n"
			<< "ZONE I=" << dat.GetGridDimensionX1() << ", J=" << dat.GetGridDimensionX2() 
			<< ", K=" << dat.GetGridDimensionX3() << ", F=POINT\n\n";
		// write nodal data
		for (int i=0;i!=dat.GetNumNode();++i)
		{
			outfile << *(dat.GetNode(i)) << endl;
		}
		outfile.close();
	}
}

void pdSolve::WriteEnergy(const pdDataManager& dat, const int step, const double time)
{
	// This function outputs total energy data.

	if (dat.GetRank()==0)
	{
		string keywd = dat.GetWkDir() + "pds_output_e_";
		keywd = dat.KeywordAppend(keywd, step);
		keywd += ".dat";
		ofstream outfile(keywd.c_str(), ios::out);
		outfile << "Total energy output mode\n";
		outfile.setf(ios::scientific);
		outfile.precision(5);
		outfile << step << "\t" << time << "\t" << dat.GetKinetTotal() << "\t" << dat.GetElastStep() << endl;
		outfile.close();
	}
}

void pdSolve::WriteEachStep(const pdDataManager& dat, const int step, const double time, ofstream& outfile)
{
	// This function writes out data at given noda for each time step.

	if (dat.GetRank()==0)
	{
		string stepFileKeyword = dat.GetWkDir() + "pds_steps.dat";
		ofstream stepFile(stepFileKeyword.c_str(), ios::app);
		stepFile.setf(ios::scientific);
		stepFile.precision(5);
		// write nodal data
		double gridSpacing = dat.GetGridSpacing();
		for (int i=0;i!=dat.GetNumNode();++i)
		{
			pdNode* node = dat.GetNode(i);
			if (node->GetX1() == 0.5 * gridSpacing && node->GetX2() == 0.5 * gridSpacing
			&& node->GetX3() == 0.5 * gridSpacing)
			{
				outfile << step << "\t" << time << "\t" << *node << endl;
			}
		}
	}
}

void pdSolve::Debug(pdDataManager& dat)
{
	// This function writes out useful data for debugging
	dat.Debug();
}

void pdSolve::PostFormatData(const pdDataManager& dat)
{
	// This function converts output files into desired formatted data file.

	int step;
	cout << " Please give the step number you want to read: ";
	cin >> step;
	string ifPath = dat.GetWkDir() + "pds_tecplot_";
	ifPath = dat.KeywordAppend(ifPath, step);
	ifPath += ".dat";
	ifstream inFile(ifPath.c_str(), ios::in);
	if (!inFile)
	{
		cerr << "ReadTecplotFile: pdMTEL tecplot output file not exist!" << endl;
		exit(0);
	}
	else
	{
		string ofPath = dat.GetWkDir() + "pds_format_";
		ofPath = dat.KeywordAppend(ofPath, step);
		ofPath += ".dat";
		ofstream outFile(ofPath.c_str(), ios::out);
		dat.PostFormatData(inFile, outFile);
		inFile.close();
		outFile.close();
	}
}

void pdSolve::PostAFD(pdDataManager & dat)
{
	int step;
	cout << " Please give the step number you want to read: ";
	cin >> step;
	string ifPath = dat.GetWkDir() + "pds_tecplot_";
	ifPath = dat.KeywordAppend(ifPath, step);
	ifPath += ".dat";
	ifstream inFile(ifPath.c_str(), ios::in);
	if (!inFile)
	{
		cerr << "ReadTecplotFile: pdMTEL tecplot output file not exist!" << endl;
		exit(0);
	}
	else
	{
		string ofPath = dat.GetWkDir() + "pds_AFD_";
		ofPath = dat.KeywordAppend(ofPath, step);
		ofPath += ".dat";
		ofstream outFile(ofPath.c_str(), ios::out);
		dat.PostAFD(inFile, outFile);
		inFile.close();
		outFile.close();
	}
}

//void pdSolve::PostFormatData(const pdDataManager& dat)
//{
//	// This function converts output files into desired formatted data file.
//
//	if (dat.GetRank()==0)
//	{
//		ifstream finput;
//		string ofpath = dat.GetWkDir() + "pds_format_data.dat";
//		ofstream foutput(ofpath.c_str(), ios::out);
//
//		string keywd = dat.GetWkDir() + "emu.plt.";
//		string end = ".0";
//		//input the time step interval
//		int step_start, step_final, step, freq;
//		cout << "start time step = " << endl;
//		cin >> step_start;
//		cout << "dump frequency = " << endl;
//		cin >> freq;
//		cout << "end time step = " << endl;
//		cin >> step_final;
//		step = (step_final-step_start)/freq+1;
//
//		for(int i=0;i<=step;++i)
//		{
//			if (step_start==0)
//			{
//				int f = i*freq;
//				if (f>step_final)
//				{
//					f = step_final;
//					++i;
//				}
//				cout << "current time step in reading is " << f << endl;
//				string emu = dat.KeywordAppend(keywd, f);
//				emu += end;
//				keywd = dat.GetWkDir() + "emu.plt.";
//				finput.open(emu.c_str(), ios::in);
//				dat.ReadEMUPlotFile(finput, foutput);
//				finput.close();
//			}
//			else if (step_start==1)
//			{
//				int f = 1;
//				if (f>step_final)
//				{
//					f = step_final;
//					++i;
//				}
//				cout << "current time step in reading is " << f << endl;
//				string emu = dat.KeywordAppend(keywd, f);
//				emu += end;
//				keywd = dat.GetWkDir() + "emu.plt.";
//				finput.open(emu.c_str(), ios::in);
//				dat.ReadEMUPlotFile(finput, foutput);
//				finput.close();
//				step_start = 0;
//			}
//			else
//			{
//				int f = step_start + i*freq;
//				if (f>step_final)
//				{
//					f = step_final;
//					++i;
//				}
//				cout << "current time step in reading is " << f << endl;
//				string emu = dat.KeywordAppend(keywd, f);
//				emu += end;
//				keywd = dat.GetWkDir() + "emu.plt.";
//				finput.open(emu.c_str(), ios::in);
//				dat.ReadEMUPlotFile(finput, foutput);
//				finput.close();
//			}
//		}
//		foutput.close();
//	}
//}
