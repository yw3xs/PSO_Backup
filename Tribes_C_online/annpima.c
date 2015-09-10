// By Rui Mendes, 2002
// Modified by Maurice Clerc (MC), 2004
//#include <math.h>

//#define sigmoid(x_)	1.0 / (1.0 + exp(-(x_)))

//#define INPUT		7
//#define HIDDEN		7
//#define OUTPUT		1
//#define DIMENSIONS	INPUT * HIDDEN + HIDDEN + HIDDEN + 1 = 64

static double netPIMA(double *inputs, struct vector weights) {
	double h[7];
	int i, j;
	int n = 0;
	double acc;

	for(i = 0; i < HIDDEN; i++) {
		for(acc = 0.0, j = 0; j < INPUT; j++)
			acc += inputs[j] * weights.x[n++];
		h[i] = sigmoid(acc + weights.x[n++]);
	}

	for(acc = 0.0, i = 0; i < HIDDEN; i++)
		acc += h[i] * weights.x[n++];

	acc=acc + weights.x[n++];
	acc=sigmoid(acc);
//	if (acc>0) acc=1; else acc=0; // Threshold (MC)

	return acc;
}

//#define CASES	200

double ANNPIMA(struct vector weights) {
	int CASES=200;
	int i; //, j;
	//int par;
	double acc = 0.0, ret;

	static double input[200][7] =
#include "pima.input.dat"
;

	static double output[200] = 
#include "pima.output.dat"
;

	INPUT=7; HIDDEN=7; OUTPUT=1;

	for(i = 0; i < CASES; i++) {

		ret = netPIMA(input[i], weights);
		ret -= output[i];

		acc += ret * ret;
	}

	return sqrt(acc/(double)CASES);
}
