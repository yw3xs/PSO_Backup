int 			alea(int min, int max);
int 			alea_diff(int min,int max, int num);
double 			alea_float(double min, double max);
double			alea_normal(double mean,double std_dev);
struct vector_i	bits(double x, int D);
int				K_optim(int N, double eps);
double 			MAX(double a,double b);
double 			MIN(double a,double b);
double			number(int d,struct vector_i x);
struct vector 	rand_in_hypersphere(int dimen, double radius,double non_unif);
struct   vector   random_permut(struct vector p);
double			regranul(double x,double granul);
struct roots solve_polyn_3(double a0,double a1,double a2,double a3);
#define sigmoid(x_)	1.0 / (1.0 + exp(-(x_)))

