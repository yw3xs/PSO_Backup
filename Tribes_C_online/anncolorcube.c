// By Rui Mendes, 2002
/* Modified by Maurice Clerc (MC), 2004
Note: The network is too complicated
Just one output would be enough, by using another colour coding
for example -1=black, 0=grey,1=white

 */

//#include <math.h>

//#define sigmoid(x_)	1.0 / (1.0 + exp(-(x_)))

//#define INPUT		3
//#define HIDDEN		8
//#define OUTPUT		3
//#define DIMENSIONS	INPUT * HIDDEN + HIDDEN + HIDDEN * OUTPUT + OUTPUT =46

 void netCC(double *inputs, struct vector weights) {
	double h[8];
	int i, j;
	int n = 0;
	double acc;
	int max;

	for(i = 0; i < HIDDEN; i++) {
		for(acc = 0.0, j = 0; j < INPUT; j++)
			acc += inputs[j] * weights.x[n++];
		h[i] = sigmoid(acc + weights.x[n++]);
	}

	for(j = 0; j < OUTPUT; j++) {
		for(acc = 0.0, i = 0; i < HIDDEN; i++)
			acc += h[i] * weights.x[n++];
		retp[j] = sigmoid(acc + weights.x[n++]);
	}
/*
	max=0; // (Added MC) Just one 1 and the others to 0
	for (j=1;j<OUTPUT;j++) if(retp[j]>retp[max]) max=j;
	for (j=0;j<OUTPUT;j++) retp[j]=0;
	retp[max]=1;
*/

}

//#define CASES	27

double ANNCOLORCUBE(struct vector weights) {
	int CASES=27;
	int i, j;
	double acc = 0.0;
	//double ret[3];

	static double input[27][3] =
#include "colorcube.input.dat"
;

	static double output[27][3] =
#include "colorcube.output.dat"
;

	INPUT=3; HIDDEN=8;OUTPUT=3;
	for(i = 0; i < CASES; i++) {

		netCC(input[i], weights);

		for(j = 0; j < OUTPUT; j++) {
			retp[j] -= output[i][j];
			acc += retp[j] * retp[j];
		}
	}

	return sqrt(acc/(double)(CASES * OUTPUT));
}
