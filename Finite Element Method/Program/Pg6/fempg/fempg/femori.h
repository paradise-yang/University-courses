#pragma once
#include<fstream>
#include<string>
#include"ConjugateGradient.h"

using namespace std;

#define epsilon1 1e-1
#define epsilon2 1e-2
#define epsilon3 1e-3
#define epsilon7 1e-7

//真实解函数u(x)；
double u(double epsilon, double x) 
{ 
	if (epsilon >= 1e-2)
	{
		double c1 = epsilon * x + x * x / 2;
		double c2 = (exp(x / epsilon) - 1) / (exp(1.0 / epsilon) - 1);
		return c1 - (epsilon + 0.5) * c2;
	}
	else
	{
		double c1 = epsilon * x;
		double c0 = x * x / 2;
		double c2 = (x - 1.0) / epsilon;
		return c0 + c1 - (epsilon + 0.5) * exp(c2);
	}
	return (0.5 + epsilon1) * (1.0 - 1.0 * exp(x / epsilon1)) / (exp(1.0 / epsilon1) - 1.0) + 0.5 * x * x + epsilon1 * x; 
}
double u_2(double x){ 
	double c1 = epsilon7 * x;
	double c0 = x * x / 2;
	double c2 = (x - 1.0) / epsilon7;
	return c0 + c1 - (epsilon7 + 0.5) * exp(c2);
}
//函数f(x)；
double f(double x){	return x;}
