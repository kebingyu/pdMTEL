#include "pdSolve.h"
#include "pdDataManager.h"
#include "windows.h"
#include <iostream>
using namespace std;

int SelectMenu();
string GetWkDir();
enum MenuOption { DYNAMIC_SOLVE=1, STATIC_SOLVE, DEBUG, POSTAFD, POSTFORMAT, EXIT};
#include "debug_new.h"
// *** memory leak detect ***
#define _CRTDBG_MAP_ALLOC  // this is used to indicate where the memory leak happens
#include<stdlib.h>
#include<crtdbg.h>

int main()
{
	_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF ); // memory leak detector will be called once exe file exits

	// initial mpi environment and working folder
	pdDataManager dat;
	dat.InitMPI();
	// choose calculation types
	int option = SelectMenu();

	// start tick counting
	DWORD tick_o = GetTickCount();
	switch(option)
	{
	case DYNAMIC_SOLVE:
		{
			pdSolve::InitialProblem(dat);
			pdSolve::SolveDynamic(dat);
			break;
		}
	case STATIC_SOLVE:
		{
			pdSolve::InitialProblem(dat);
			pdSolve::SolveStatic(dat);			
			break;
		}
	case DEBUG:
		{
			pdSolve::InitialProblem(dat);
			pdSolve::Debug(dat);
			break;
		}
	case POSTAFD:
		{
			pdSolve::InitialProblem(dat);
			pdSolve::PostAFD(dat);
			break;
		}
	case POSTFORMAT:
		{
			pdSolve::InitialProblem(dat);
			pdSolve::PostFormatData(dat);
			break;
		}
	default:
		{
			if (dat.GetRank()==0)
			{
				string ofpath = dat.GetWkDir() + "pds_output.dat";
				ofstream outfile(ofpath.c_str(), ios::app);
				outfile << endl << "user quits without doing any calculation!" << endl;
				outfile.close();
			}
			//MPI::Finalize();
			return 0;
		}
	}

	// find out the total run time
	int total_tick = GetTickCount() - tick_o;
	int hour = total_tick/3600000;
	int minute = (total_tick/1000 - hour*3600)/60;
	int second = (total_tick/1000 - hour*3600 - minute*60);
	if (dat.GetRank()==0)
	{
		string ofpath = dat.GetWkDir() + "pds_output.dat";
		ofstream outfile(ofpath.c_str(), ios::app);
		outfile << "The total CPU tick = " << total_tick << " (millisecond)" << endl;
		outfile << "The total CPU time is about " << hour << " hours, " << minute << " minutes, and " << second << " seconds" << endl;
		outfile.close();
	}
	//MPI::Finalize();

//	_CrtDumpMemoryLeaks(); // memory leak detector
	int c;
	c = cin.get();
	return 0;
}

int SelectMenu()
{
	int i;
	//if (MPI::COMM_WORLD.Get_rank() == 0)
	{
		cout << " Please choose an option to run the program: \n"
			<< "   - Solve the problem using dynamic solution (1) \n"
			//<< "   - Solve the problem using static solution (2) \n"
			<< "   - Debug (3) \n"
			<< "   - Post: Read output file and calculate AFD on given node (4) \n"
			<< "   - Post: Convert output file to formatted data file (5) \n"
			<< "   - Any other numbers to exit " << endl
			<< " Your input: ";
		cin >> i;
	}
	return i;
}

