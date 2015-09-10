#include "headfile.h"
#include "extern.h"
#include "myfun.h"

float f6(int a)
{
/*
	This is the f6 function as described in the Handbook of
	Genetic Algorithms, p.8
*/
	double num, denom, f6;
	double errorf6;

	num=(sin(sqrt((xx[0][a]*xx[0][a])+(xx[1][a]*xx[1][a]))))  *
		 (sin(sqrt((xx[0][a]*xx[0][a])+(xx[1][a]*xx[1][a])))) - 0.5;
	denom=(1.0 + 0.001 * ((xx[0][a] * xx[0][a]) + (xx[1][a]*xx[1][a]))) *
		(1.0 + 0.001 * ((xx[0][a] * xx[0][a]) + (xx[1][a]*xx[1][a])));
	f6=(double) 0.5 - (num/denom);
	errorf6=1 - f6;
	return errorf6;
}

float sphere(int a, int b)
{
	/* This is the familiar sphere model
		int a: index of particles   b:dimension */

	double result;
	int i;

	result=0.0;

	for (i=0;i<b;i++)
	{
		result += xx[i][a]*xx[i][a];
	}

	return result;
}

float rosenbrock(int a, int b)
{
	/* this is the Rosenbrock function
		a: index of the particles; b:dimension */

	int i;
	double result;

	result=0.0;

	for (i=1;i<b;i++)
	{
		result +=100.0*(xx[i][a]-xx[i-1][a]*xx[i-1][a])*(xx[i][a]-xx[i-1][a]*xx[i-1][a]) + (xx[i-1][a]-1)*(xx[i-1][a]-1);
	}

	return fabs(result);
}

float rastrigrin(int a, int b)
{
	/* This is the generalized Rastrigrin function
		a:index of the particles; b:dimension */

	int i;
	double result;

	result=0.0;

	for (i=0;i<b;i++)
	{
		result +=xx[i][a]*xx[i][a] - 10.0*cos(2.0*3.141591*xx[i][a])+10.0;
	}

	return result;
}

float griewank(int a,int b)
{
	/* This is the generalized Griewank function
		a:index of the particles; b:dimension */

	int i;
	double result_s,result_p;

	result_s=0.0;
	result_p=1.0;

	for (i=0;i<b;i++)
	{
		result_s +=xx[i][a]*xx[i][a];
		result_p *=cos(xx[i][a]/sqrt(i+1));
	}
	result_s =result_s/4000.0 - result_p +1;

	return result_s;
}


