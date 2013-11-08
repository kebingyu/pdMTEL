#include <ctime>
#include <cstdio>   
#include "math.h"   
#include <cstdlib>  
#pragma once;

const unsigned long maxshort=65536L; 
const unsigned long multiplier=1194211693L; 
const unsigned long adder=12345L; 

class RandomNumber 
{ 
private: 
	// data
	unsigned long _randSeed; 

public: 
	// constructor
	RandomNumber(unsigned long s); // seed generator, default means system gives seed  

	// function
	unsigned short Random(unsigned long n);
	double RandBetween(double r1, double r2);
};

RandomNumber::RandomNumber(unsigned long s) 
{ 
	if (s==0)
	{
		_randSeed = (int)time(0); // seed generates according to system time
	}
	else
	{
		_randSeed = s; // user-defined seed
	}
} 

unsigned short RandomNumber::Random(unsigned long n) 
{
	// generate a random unsigned short from 0 to n-1
	_randSeed = multiplier*_randSeed + adder; 
	return (unsigned short)((_randSeed>>16)%n); 
}

double RandomNumber::RandBetween(double r1, double r2) 
{ 
	// generate a random double from r1 to r2
	double rand = Random(maxshort)/double(maxshort); // random float from zero to one
	return rand*(r2-r1)+r1;
} 
 


