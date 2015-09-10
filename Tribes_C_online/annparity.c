// By Rui Mendes, 2002
//#include <math.h>

//#define sigmoid(x_)	1.0 / (1.0 + exp(-(x_)))

//#define INPUT		4
//#define HIDDEN		4
//#define OUTPUT		1
//#define DIMENSIONS	INPUT * HIDDEN + HIDDEN + HIDDEN + 1  = 25

static double netNBP(double *inputs, struct vector weights)
{
	double h[4];
	int i, j;
	int n = 0; // Arc rank
	double acc;

	for(i = 0; i < HIDDEN; i++)
  {
    acc=0;
    for(j = 0; j < INPUT; j++)
			acc += inputs[j] * weights.x[n++];
		h[i] = sigmoid(acc + weights.x[n++]);
	}
  acc=0;
	for(i = 0; i < HIDDEN; i++)
		acc += h[i] * weights.x[n++];
acc= sigmoid(acc + weights.x[n++]);

//printf("\n weights: ");for (i=0;i<4;i++) printf(" %f",weights.x[i]);printf("=> %f",acc);

	return acc;
}

//#define CASES	16

double ANNParity4(struct vector weights)
{
 int CASES;
  int i, j;
	int par;
	double acc = 0.0, ret;
	double input[4];

 CASES=16;
  //CASES=nprogr;
  INPUT=4; HIDDEN=4; OUTPUT=1;

	for(i = 0; i < CASES; i++) {
		for(j = 0; j < INPUT; j++)
			input[j] = ( i & (1 << j) ) >> j;

		ret = netNBP(input, weights);
		for(j = 0, par = 0; j < INPUT; j++)
			par += input[j];

		ret -= par % 2;

		acc += ret * ret;

	}
  
  acc=  sqrt(acc/(double)CASES);
  //printf("\n weights: ");for (i=0;i<4;i++) printf(" %f",weights.x[i]);printf("=> %f",acc);
	return acc;
}
