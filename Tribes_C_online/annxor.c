// By Rui Mendes, 2002
//#include <math.h>

//#define sigmoid(x_)	1.0 / (1.0 + exp(-(x_)))

//#define INPUT		2
//#define HIDDEN		2
//#define OUTPUT		1
//#define DIMENSIONS	INPUT * HIDDEN + HIDDEN + HIDDEN + 1 = 9

static double netXOR(double *inputs, struct vector weights) {
	double h[2];
	int i, j;
	int n = 0;
	double acc;

	for(i = 0; i < 2; i++) {
		for(acc = 0.0, j = 0; j < 2; j++)
			acc += inputs[j] * weights.x[n++];
		h[i] = sigmoid(acc + weights.x[n++]);
	}

	for(acc = 0.0, i = 0; i < 2; i++)
		acc += h[i] * weights.x[n++];

	return sigmoid(acc + weights.x[n++]);
}

//#define CASES	4

double ANNXor(struct vector weights) {
	int i, j;
	int par;
	double acc = 0.0, ret;
	double input[2];

	for(i = 0; i < 4; i++) {
		for(j = 0; j < 2; j++)
			input[j] = ( i & (1 << j) ) >> j;

		ret = netXOR(input, weights);

		for(j = 0, par = 0; j < 2; j++)

			par += input[j];

		ret -= par % 2;

		acc += ret * ret;
	}

	return sqrt(acc/(double)4);
}
