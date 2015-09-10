// By Rui Mendes
//#include <math.h>

//#define sigmoid(x_)	1.0 / (1.0 + exp(-(x_)))

//#define INPUT		1
//#define HIDDEN		8
//#define OUTPUT		1
//#define DIMENSIONS	INPUT * HIDDEN + HIDDEN + HIDDEN + INPUT + 1  = 26

static double netSTS(double *inputs, struct vector weights) {
	double h[8];
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

//#define CASES	80

double ANNSINSIMP(struct vector weights) {
   int CASES=80;
  int i;
	double acc = 0.0, ret;
   
	static double input[80][1] =
#include "sinsimp.input.dat"
;

	static double output[80] = 
#include "sinsimp.output.dat"
;

  INPUT=1;HIDDEN=8; OUTPUT=1;
  
	for(i = 0; i < CASES; i++) {

		ret = netSTS(input[i], weights);
		ret -= output[i];

		acc += ret * ret;
	}

	return sqrt(acc/(double)CASES);
}
