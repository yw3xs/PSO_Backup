// By Rui Mendes
//#include <math.h>

//#define sigmoid(x_)	1.0 / (1.0 + exp(-(x_)))

//#define INPUT		4
//#define HIDDEN		4
//#define OUTPUT		1
//#define DIMENSIONS	INPUT * HIDDEN + HIDDEN + HIDDEN + INPUT + 1 = 29

static double netSERVO(double *inputs, struct vector weights) {
	double h[4];
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

	for(i = 0; i < INPUT; i++)
		acc += inputs[i] * weights.x[n++];

	return acc + weights.x[n++];
}

//#define CASES	167

double ANNSERVO(struct vector weights) {
	int CASES=167;
	int i;
	double acc = 0.0, ret;

	static double input[167][4] =
#include "servo.input.dat"

	static double output[167] = 
#include "servo.output.dat"

	INPUT=4; HIDDEN=4; OUTPUT=1;
	for(i = 0; i < CASES; i++) {

		ret = netSERVO(input[i], weights);
		ret -= output[i];

		acc += ret * ret;
	}

	return sqrt(acc/(double)CASES);
}
