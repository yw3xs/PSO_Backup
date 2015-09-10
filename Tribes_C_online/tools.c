// ------------------------------------------------------------------- ALEA
int alea(int min, int max) /* Random integer number between min and max */
{
int		ir;
double 	r;

r=alea_float(0,1);

ir=(int)(min+r*(max+1-min));
if (ir>max) ir=max;
return ir;
}


/*------------------------------------------------------------------- ALEA_DIFF
*/
int alea_diff(int min,int max, int num) /* Generate randomly an integer number #

num */
{
int	temp1,temp2;

if (num==min)
	{
	temp1=alea(min+1,max);
	return temp1;
	}

if (num==max)
	{
	temp1=alea(min,max-1);
	return temp1;
	}

temp1=alea(min,num-1);
temp2=alea(num+1,max);

if (alea(0,100)<50)

	{
	return temp1;
	}
else

	{
	return temp2;
	}
}

// ---------------------------------------------------------------------- ALEA_FLOAT
double alea_float(double a, double b) /* Random  number between a and b */
{

double 	r;


// Normally, RAND_MAX = 32767 = 2^15-1
// Not good. Replace it by KISS (for example)

//r=(double)rand()/RAND_MAX; // r in [0,1]
//return a+r*(b-a);

	r=(double)rand_kiss();
	r=r/RAND_MAX_KISS;
//printf("\n %f",r);
	r=a+r*(b-a);
	return r;
}

 // -------------------------------------------------------------------------- ALEA_NORMAL
  double alea_normal(double mean, double std_dev)
 {
 /*
     Use the polar form of the Box-Muller transformation to obtain
     a pseudo random number from a Gaussian distribution
 */
        double x1, x2, w, y1;
        //double  y2;

         do {
                 x1 = 2.0 * alea_float(0,1) - 1.0;
                 x2 = 2.0 * alea_float(0,1) - 1.0;
                 w = x1 * x1 + x2 * x2;
         } while ( w >= 1.0);

         w = sqrt( -2.0 * log( w ) / w );
         y1 = x1 * w;
        // y2 = x2 * w;
          y1=y1*std_dev+mean;
         return y1;
  }

//============================================================= BITS
struct vector_i	bits(double x, int D)
{
// Gives a binary representation of x
	struct vector_i bit;
	int				d;
	int			y;

	for (d=0;d<D;d++) bit.x[d]=0;

	y=(int)(x+0.5);

//printf("\n %f %i = ",x,y);
	for(d=0;d<D;d++)
	{
		bit.x[d]=y%2;
		if (y==0 || y==1) {goto end;}
		y=(y-bit.x[d])/2;

	}

end:
bit.size=D;


//for (d=0;d<D;d++) printf("%i ",bit.x[d]);

	return bit;
}

//============================================================= COEFF_SC
double coeff_SC(int D)
{

	//  The D-cube whose edge is a
	// and  the D-sphere whose radius is r=a*coeff_S_C
	// have the same volume
	double coeff;
	double	c1,c2;
	int		d;
	double	d3;
	double	x;

	int option=0; // 0 => exact value
				// 1 => mean value

if (D==1) return 0.5;

	x=(double)D;
	if ((2*(int)(D/2)==D) || option==1) // D even
	{
		d3=1; for (d=2;d<=D/2;d++) d3=d3*d; // (D/2)!
		c1=pow(d3,1/x)/sqrt(pi);
		if (option==0) return c1;
	}

	 // D odd
	{
		d3=1; for (d=2;d<=D;d++) d3=d3*d; // D!
		c2=0.5*pow(d3,1/x)/pow(pi,0.5-0.5/x);
		if (option==0) return c2;
	}

	coeff=(c1+c2)/2;

return coeff;
}
 // ----------------------------------------------------------------------------- K_OPTIM
int	K_optim(int N,double eps)
{
// Compute the optimal number of random information links
int K;
//printf("\n K_optim. N %i, eps %f", N, eps); printf("\n");
if(N<=1) return 1;

if (eps<1/infinite) eps=1/infinite;

K=0.5+(2/(double)N)*log(eps)/log(1-1/(double)N);
//printf("\n K_optim %i",K);
return K;
}


// ------------------------------------------------------------------------- MAX
double MAX(double a,double b)
{
if (a>b) return a; return b;
}

// ------------------------------------------------------------------------- MIN
double MIN(double a,double b)
{
if (a<b) return a; return b;
}

//===========================================================
double number(int d,struct vector_i x)
{
// Compute the value (base 10) of the bit string

  int D;
  double z;
  D=x.size;
  if (d==D-1) return x.x[D-1];
  z=x.x[d]+2*number(d+1,x);
  return z;
}

// ------------------------------------------------------------------------- RAND_IN_HYPERSHERE
struct vector rand_in_hypersphere(int D, double radius,double non_unif)
{
/*  ******* Random point in a hypersphere ********
 Maurice Clerc 2003-07-11

Put  a random point inside the hypersphere S(0,1) (center 0, radius 1).

Uniform distribution if non_unif=1.
Normal (Gaussian) distribution if non_unif<0 (the abs value is then the standard

deviation)
*/

int 	j;
double   length;
double      pw;
double      r;
 double     sigma;
struct	vector	x;

x.size=D;
pw=1/(double)D;

// ----------------------------------- Step 1.  Direction
    length=0;
   for (j=0;j<D;j++)
   {
          x.x[j]=alea_normal(0,1);
          length=length+  x.x[j]*x.x[j];
   }

   length=sqrt(length);

 //----------------------------------- Step 2.   Random radius
   if (non_unif>0)  // Pure hyperspherical distribution. Uniform is non_unif=1/
 {
    r=alea_float(0,1);
    r=pow(r,pw*non_unif);
  }
  else
  {
        sigma=-non_unif;
        r=fabs(alea_normal(0,sigma));
  }

       for (j=0;j<D;j++)
   {
          x.x[j]=radius*r*x.x[j]/length;
   }
return x;
}
//============================================================= RANDOM_PERMUT
 struct   vector   random_permut(struct vector p)
 {
    int                j, j0;
    int                 jmax;
    int                 k;
    struct vector p0,p1;

    p1.size=p.size;
    p0=p;
// Random permutation of the coordinates

	k=0;
	jmax=p.size-1;

	next_k:

	j0=alea(0,jmax); // Random integer in [0,...,jmax]
	p1.x[k]=p0.x[j0];

	// Compact p0
	if (j0<jmax)
		{
		for (j=j0;j<jmax;j++)
			{
			p0.x[j]=p0.x[j+1];
			}
		}
	jmax=jmax-1;
	if (jmax>=0)
		{
		k=k+1;
		goto next_k;
		}
      return p1;
 }

 // ------------------------------------------------------------------------- REGRANUL
double	regranul(double x,double granul)
{ // Modify a value according to the granularity of the search space
	double xp;
if (granul<almostzero) return x;// Pseudo-continuous
								// (depending on the machine)

if (x>=0) xp=granul*floor (x/granul+0.5); //"1/granul" is the minimum
												// distance between two points
else xp= (-granul*floor(-x/granul+0.5));

return xp;

}

//================================================== KISS
/*

 the idea is to use simple, fast, individually promising
 generators to get a composite that will be fast, easy to code
 have a very long period and pass all the tests put to it.
 The three components of KISS are
        x(n)=a*x(n-1)+1 mod 2^32
        y(n)=y(n-1)(I+L^13)(I+R^17)(I+L^5),
        z(n)=2*z(n-1)+z(n-2) +carry mod 2^32
 The y's are a shift register sequence on 32bit binary vectors
 period 2^32-1;
 The z's are a simple multiply-with-carry sequence with period
 2^63+2^32-1.  The period of KISS is thus
      2^32*(2^32-1)*(2^63+2^32-1) > 2^127
*/



static ulong kiss_x = 1;
static ulong kiss_y = 2;
static ulong kiss_z = 4;
static ulong kiss_w = 8;
static ulong kiss_carry = 0;
static ulong kiss_k;
static ulong kiss_m;



void seed_rand_kiss(ulong seed)
{
    kiss_x = seed | 1;
    kiss_y = seed | 2;
    kiss_z = seed | 4;
    kiss_w = seed | 8;
    kiss_carry = 0;
}

ulong rand_kiss()
{
    kiss_x = kiss_x * 69069 + 1;
    kiss_y ^= kiss_y << 13;
    kiss_y ^= kiss_y >> 17;
    kiss_y ^= kiss_y << 5;
    kiss_k = (kiss_z >> 2) + (kiss_w >> 3) + (kiss_carry >> 2);
    kiss_m = kiss_w + kiss_w + kiss_z + kiss_carry;
    kiss_z = kiss_w;
    kiss_w = kiss_m;
    kiss_carry = kiss_k >> 30;
    return kiss_x + kiss_y + kiss_w;
}

//================================================== SOLVE_POLYN_3
struct roots solve_polyn_3(double a0,double a1,double a2,double a3)
{
/*
input: a0 + a1*x + a2*x*x + a3*x*x*x
output: real root(s). One or three.
  */
struct roots z;

double	alpha;
double	b0,b1,b2;
double	c1,c2;
double	discr;
double	f1,f2;
double	p,q,r,s;
double	rho,theta;
double	z1;

b0=a0/a3; b1=a1/a3; b2=a2/a3;
alpha=b2/3;
p=b1-3*alpha*alpha;
q=b0-alpha*(p+alpha*alpha);

r=q*q/4+p*p*p/27;
if (r<0)
{
	s=-r;
	rho=sqrt(q*q/4+s);
	theta=acos(q/(2*rho));
	z1=-pow(rho,1.0/3);
	z1=z1*2*cos(theta/3);

	c1=-p/z1-q/(z1*1);
	c2=-q/z1;
	discr=-c2+c1*c1/4;
	z.z[1]=-c1/2+sqrt(discr)-alpha;
	z.z[2]=-c1/2-sqrt(discr)-alpha;
	z.z[0]=z1-alpha;
	z.status=3; // Three real roots
	}
	else
	{
	r=sqrt(r);
	f1=-q/2 + r;
	if (f1>=0) f1=pow(f1,1.0/3);
	else f1=-pow(-f1,1.0/3);

	f2=-q/2 - r;
	if (f2>=0) f2=pow(f2,1.0/3);
	else f2=-pow(-f2,1.0/3);
	z.z[0]=f1+f2-alpha;
	z.status=1; // Just one real root
	}
printf("\n status %i",z.status);
return z;
}
