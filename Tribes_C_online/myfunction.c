// ----------------------------------------------------------------------------- MYFUNCTION
/*
   You can add any function you want.
   Don't forget to also modify functions.txt
*/
struct f MyFunction(struct position pos,int funct,double target,int level)
{

// Position evaluation

double	a_1,a_2,a_3;
double	alpha,beta,gamma,delta,p,q;
double	better;
double	b;
struct vector_i bit,bit_x;
int		bit_sum;
double   c,c1,c2;
int			choice;
int		d,dmax;
double	discr;
int		D;
int 			DD;
double	f_model;
double   f;
double      f1,f2;
double	gen[Max_DD]; // For Moving Peaks
int		i;
int		ix,ix2;
int     j;
int     k;
int		m;
int		n1,n2;
int		num;
//struct position post={0};
struct position post;
double      pr;
double			product;
double      r;
double			rho;
double			s;
double  sigma;
double			sum1,sum2;
double			theta;
struct f 			total;
int		used[Max_DD]={0};
double	x1,x2,x3,x4,x5,x6,x7,x8,x9;
double			x;
double      xd;
double			xid;
double	y,y1,y2,y3,y4;
double	z1,z2,z3;
struct roots z;

// Data for model tuning example (see case 22)

static float 	data[5][2] =
{
{	10	,	1480	},
{	20	,	2200	},
{	30	,	5000	},
{	40	,	5880	},
{	50	,	11590	}
};


static double	Fmax=1000.0;
static double	S=189000.0;
static double	lmax=14.0;
static double	dmin=0.2;
static double	Dmax=3.0;
static double	Fp=300;
static double	spm=6.0;
static double	sw=1.25;
static double	G=11500000;

// Variables specific to Coil compressing spring
// See also constrain.c
double Cf,K,sp,lf;

//----------------
int d2=5;

// Data for Master Mind  (hard coded solution)
static int data_M[4] =
{ 1,2,3,4};

/*
// Specific Informational table for 4 different colors
static int r_w[5][5]=
{

{16	,192,	120,	20,	1},
{152,	192,	24,	0,	0},
{312,	108,	6,	0,	0},
{136,	8,	0,	0	,0},
{9	,0,	0	,0,	0}
};
 */
// For Foxholes problem
static int a[2][25] =
{
   {-32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32,
      -32, -16, 0, 16, 32, -32, -16, 0, 16, 32  },
   {-32, -32, -32, -32, -32, -16, -16, -16, -16, -16,
      16, 16, 16, 16, 16, 32, 32, 32, 32, 32 }
};

// For Jeannet-Messine problem (cf. ROADEF 2003, p 273)
static double a1[6]=
{ 0.5, 0.3, 0.8, 0.1, 0.9, 0.12};
static double a2[6]=
{-0.5, 0.6, 0.1, 1.5, -1, 0.8};

// For polynom fittint
   int const M=60;
   double py, dx=(double)M;
 //-----------------------------------------
 
eval_f[level]=eval_f[level]+1;
eval_f_tot[level]=eval_f_tot[level]+1;

if (problem[level].printlevel>2)
	printf("\n my_function. tot. eval_f_tot %6.0f",eval_f_tot[level]);

DD=pos.p.size;
total.size=pos.f.size;

switch (funct)
{

case 1:  //DD-sum
	total.f[0] = 0;
	for( d=0;d<DD;d++)                       //for each dimension
	{
		total.f[0] = total.f[0] + pos.p.x[ d];               //add DD if i's coordinates
	}

	break;

case 2:  // DD-product.
	total.f[0] = 1;
	for( d=0;d<DD;d++)
		total.f[0] = total.f[0]*pos.p.x[ d];

	break;

case 3: //Parabole () Sphere
	//      Global min f(x)=0, x(i)=0
	total.f[0] = 0;
	for( d=0;d<DD;d++)
	{
		total.f[0] = total.f[0] + pos.p.x[ d]*pos.p.x[ d];
	}
	break;

case 4: // Rosenbrock's valley, Banana function
//  Global min f(x)=0, x(i)=1
	total.f[0]=0;

	for (d=0;d<DD-1;d++)
		{
		xid=1-pos.p.x[d];
		total.f[0]=total.f[0]+xid*xid;
		xid= pos.p.x[d]*pos.p.x[d]-pos.p.x[d+1];
		total.f[0]=total.f[0]+100*xid*xid;
		}
	break;


case 5: // Clerc's f1, Alpine function, min 0 in (0,...,0)
	total.f[0]=0;
	/*
	for( d=0;d<DD;d++)
		{
		xid=sin(pos.p.x[d]);
		total.f[0]=total.f[0]+sqrt(fabs(pos.p.x[d]*xid));
		}
		*/
	for( d=0;d<DD;d++)
	{
		xid=pos.p.x[d];
		total.f[0]=total.f[0]+fabs(xid*sin(xid)+0.1*xid);
	}
	break;

case 6: // Griewank's function. Min =0 at (100,100 ... 100)
	if (DD==2) goto Gr2;
	total.f[0]=0;
	product=1;
	for (d=0;d<DD;d++)
		{
		xid=pos.p.x[ d]-100;
		total.f[0]=total.f[0]+xid*xid;
		product=product*cos (xid/sqrt(d+1));
		}
	total.f[0]=total.f[0]/4000 -product +1;
	break;

Gr2:	// -------- 2D only 
	x1=pos.p.x[0];
	x2=pos.p.x[1];
	total.f[0]=(x1*x1+x2*x2)/2 - cos(two_pi*x1)*cos(two_pi*x2)+1;
	break;


case 7:  // Rastrigin. Minimum value 0. Solution (0,0 ...0)
	k=10;

	total.f[0]=0;
	for (d=0;d<DD;d++)
		{
		xid=pos.p.x[ d];
		total.f[0]=total.f[0]+xid*xid - k*cos(two_pi*xid);
		}
	total.f[0]=total.f[0]+DD*k;

	break;

case 8:  // Fifty-fifty problem
	total.f[0]=tot_fifty_fifty(pos);
	break;


case 9:  // Ackley
	sum1=0;
	sum2=0;
	x=DD;
	for (d=0;d<DD;d++)
		{
		xid=pos.p.x[ d];
		sum1=sum1+xid*xid;
		sum2=sum2+cos(two_pi*xid);
		}
	total.f[0]=(-20*exp(-0.2*sqrt(sum1/x))-exp(sum2/x)+20+E);
	break;


case 10: //Foxholes 2D
	total.f[0]=0;
	for (j=0; j<25; j++)
   		{
		sum1=0;
      	for (d=0; d<2; d++)
      		{
         	sum1 =sum1+ pow (pos.p.x[d] - a[d][j],6);
      		}
      	total.f[0]=total.f[0]+1/(j+1+sum1);
   		}
	total.f[0]=1/(0.002+total.f[0]);
	break;

case 11: // ======================= Apple trees
	total.f[0]=apple_trees(pos);
	post= homogen_to_carte(pos); // For Leonardo visualization (for example)
	break;


case 12: // ======================== El Farol
/*
1 dimension = 1 Irishman
Just two values/dim.: 0 <=> stay at home
                      1 <=> go to the pub

Note: the current swarm sw is a global variable
*/
	sum1=60; // Maximum "tolerance": maximum number of people each Irishman does accept in the pub
	total.f[0]=0;
	for (d=0;d<DD;d++)
	{
		total.f[0]=total.f[0]+pos.p.x[d];
	}

	total.f[0]=total.f[0]*(1-total.f[0]/sum1);
	break;

 case 13: //  ======================Fermat

/* Suggested parameters:
- dimension 3
- xmin 2
- xmax 100
- target 0
- granul 1
- eps 0





You will then find three integer numbers x,y,z
so that x^2 + y^2 = z^2
*/

total.f[0] = 0;
dmax=pos.p.size-1;
for( d=0;d<dmax;d++)  //for each dimension

	{
	total.f[0]=total.f[0]+pos.p.x[d]*pos.p.x[d];

	}
total.f[0]=fabs(total.f[0]-pos.p.x[dmax]*pos.p.x[dmax]);
break;

case 14: //  ===========================================Knapsack
/* Suggested parameters:
- dimension 10
- xmin 1

- xmax 100
- target 100
- granul 1
- eps 0
*/
total.f[0] = 0;
	for( d=0;d<pos.p.size;d++)
		total.f[0] = total.f[0] + pos.p.x[ d];
break;


case 15: //  ===================================== (x*x+y*y-1)*(x*x+y*y-1)+(sin(10*x)-y)*(sin(10*x)-y)
// dimension =2
x=pos.p.x[0];

y=pos.p.x[1];
total.f[0]=(x*x+y*y-1)*(x*x+y*y-1)+(sin(10*x)-y)*(sin(10*x)-y);
break;

case 16: //  =====================================  Intersection of two circles
//dimension =2
x=pos.p.x[0];
y=pos.p.x[1];
total.f[0]=fabs(x*x+y*y-1)+fabs((x-2)*(x-2)+y*y-4);
break;

case 17: //  =====================================  fabs(sin(x)-y-sqrt(3)/2+1)+fabs(x*x-log(y)-(pi/3)*(pi/3))
// dimension =2 , in [0,10]
x=pos.p.x[0];
y=pos.p.x[1];
total.f[0]=fabs(sin(x)-y-sqrt(3)/2+1)+fabs(x*x-log(y)-(pi/3)*(pi/3));
break;

case 18: //  ===================================== 2D Linear system
// dimension =2, solution (5,2)
x=pos.p.x[0];
y=pos.p.x[1];
total.f[0]=fabs(3*x+2*y-19)+fabs(2*x-y-8);
break;

case 19: //  ======================================== sin wave
// dimension 1
x=pos.p.x[0];
total.f[0]=sin(x);
break;


case 20: //  =================== 3D linear system. Chinese problem of the 100 fowls
// Several integer solutions
//printf("\n\n");
	total.f[0]=fabs (pos.p.x[0]+pos.p.x[1]+pos.p.x[2]-100);
	total.f[0]=total.f[0]+fabs (5*pos.p.x[0]+3*pos.p.x[1]+(1/3)*pos.p.x[2]-100);
break;


case 21: // ===================================== Magic square
d2=(int)(sqrt(pos.p.size)); // Nb of lines = Nb of columns=sqrt(dimension)
total.f[0]=0;

for (i=0;i<d2-1;i++)  // Rows
					// For each pair of rows, compute the difference x
					// of the sums. "distance" to minimize =x*x
{

	for (j=i+1;j<d2;j++) 
	{
		x=0;	
		for (d=0;d<d2;d++)
			x=x+pos.p.x[i*d2+d]-pos.p.x[j*d2+d];
	}
	total.f[0]=total.f[0]+x*x;
}

for (i=0;i<d2-1;i++)  // Columns
{
	for (j=i+1;j<d2;j++)
	{
		x=0;
		for (d=0;d<d2;d++)
		x=x+pos.p.x[i+d*d2]-pos.p.x[j+d*d2];
		total.f[0]=total.f[0]+x*x;
	}
}

break;

case 22: // ======================== Model tuning ================
/* We have some data and a 2 parameters model
Find the "best" set of parameters
*/

/*
Model lambda*(D^mu)
column 1: D
column 2: function value
*/

total.f[0]=0;
for (d=0;d<d2;d++)
	{
	f_model=pos.p.x[0]*pow(data[d][0],pos.p.x[1]);
	total.f[0]=total.f[0]+(f_model-data[d][1])*(f_model-data[d][1]);
	}
 total.f[0]=sqrt(total.f[0]);  
break;

case 23: // =================== 4 positions 6 colors Master Mind

/*
Just a test. PSO alone is quite bad for this problem.
If you are really interested in Master Mind, you should have a look at my "ultimate" program,
on my math web site: you can't beat it.
*/

// "Normal" estimation
n1=0; // right color, right position
n2=0; // right color, wrong position


for (d=0;d<4;d++) //Right position
{
	ix=(int)pos.p.x[d];
	if (data_M[d]!=ix) continue;
		 n1=n1+1;// Right position
		used[ix]=1;
}


for (d=0;d<4;d++) // Wrong position
{
	ix=(int)pos.p.x[d];
	if (used[ix]>0) continue;

	for (d2=0;d2<4;d2++)
	{
		if (d2==d) continue;
			ix2=data_M[d2];
			if (ix2!=ix) continue;
				if (used[ix2]>0) continue;
					n2=n2+1;
					used[ix]=1;
	}
}

total.f[0]=400-100*n1 - n2;
break;


case 24: //======================= Catalan's conjecture  x^m - y^n = 1 has just one integer solution (3^2 - 2^3 = 1)
/* dimension = 2x2 =4
  a) the solution (3,2,2,3) is easily found
  b) you may try as many times as you want, with any initial position, you'll never find another one

    => the conjecture is _probably_ true
*/
total.f[0]=pow(pos.p.x[0],pos.p.x[1])-pow(pos.p.x[2],pos.p.x[3])-1;
if (pos.p.x[0]==pos.p.x[2]) total.f[0]=total.f[0]+pow(pos.p.x[0],pos.p.x[1]); // Just to avoid the local minimum x0=x2 and x1=x3
break;

case 25:
break;

case 26: //================================ Chinese-MIT problem

a_1=100;

a_2=200;
a_3=300;
total.f[0] = 0;
	for( d=0;d<pos.p.size;d++)
		total.f[0] = total.f[0]+ pos.p.x[d]; // Any positive function


n1=(int)(float)pos.p.size/3;
n2=2*n1;

x1=0;
for (d=0;d<n1;d++) x1=x1+pos.p.x[d];
x2=0;
for (d=n1-1;d<n2;d++) x2=x2+pos.p.x[d];


x3=0;
for (d=n2-1;d<pos.p.size;d++) x3=x3+pos.p.x[d];

x1=(a_1-x1);
x2=(a_2-x2);
x3=(a_3-x3);
total.f[0]=total.f[0]+x1*x1+x2*x2+x3*x3;
//total.f[0]=x1*x1+x2*x2+x3*x3;

break;

case 28: // ================================== Test surface 2D
	x1=pos.p.x[0];
	x2=pos.p.x[1];
	total.f[0]=20*x1*sin(3*x1*x2)+5*x2*sin(2*x1);

	break;


case 29: // =================          Cognitive harmony
/*
P  = Weight (square) matrix. Weights in [-1,1]
x = Activation vector. Activations in [0,1]
=> harmony = transp(x)Px, to maximize

In fact, we compute here a dissonance to minimize
*/

// Evaluate the position
MM=pos.p.size;
total.f[0]=0;
for (d=0;d<MM;d++) // Compute harmony
{
	for(m=0;m<MM;m++) {total.f[0]=total.f[0]+pos.p.x[d]*problem[level].P.val[d][m]*pos.p.x[m];}
}
total.f[0]=MM*(MM-1)-total.f[0]; // Dissonance
 //total.f[0]=MM-total.f[0];    //**** WARNING, THIS MAY GIVE NEGATIVE RESULT
break;

 case 30: // Non specific TSP
 /*
   Just a test, can't be very good. See specific PSO_for_TSP for better result.
   For this problem, you must use
   granularity=1 (normal) or 0 (continuous)
   all_different=1
 */
// Evaluate the position
total.f[0]=f_tsp(pos,1,1,-1,-1,level);
 break;

 case 31: // Sum of absolute values
   // Note: you have a classical integer problem by setting granularity to 1, in the problem description
   	total.f[0] = 0;
	for( d=0;d<DD;d++)                       //for each dimension
	{
		total.f[0] = total.f[0] + fabs(pos.p.x[ d]);               //add DD if i's coordinates
	}
 break;

  case 32: // Non specific QAP
 /*
   Just a test, can't be very good.
   For this problem, you must use
   dimension = MM
   granularity=1 (normal)
   all_different=1
 */

// Evaluate the position

total.f[0]=f_qap(pos,1,1,-1,-1,level);
//printf("\nmyfunction. total %f",total.f[0]);

 break;



 case 33: //  Multiobjective Lis-Eiben 1
 /*
    Each run with a different random initialization
    gives a point whose (f1,f2) is near of the Pareto front
*/
x1=pos.p.x[ 0];
x2=pos.p.x[ 1];
x1=x1*x1+x2*x2;
f1=pow(x1,1.0/8);

x1= pos.p.x[ 0]-0.5;
x2= pos.p.x[ 1]-0.5;
x2=x1*x1+x2*x2;

f2=pow(x2,1.0/4);

total.f[0]=f1;
total.f[1]=f2;
break;

 case 34:   // Multiobjective Schaffer
 // x in (-5,10]
    x=pos.p.x[ 0];
    f2=(x-5)*(x-5);
    if (x<=1) {f1=-x;  goto end_34;}
    if (x<=3) {f1=-2+x;  goto end_34; }
    if (x<=4) {f1=4-x;  goto end_34;}
    f1=-4+x;
    end_34:
 total.f[0]=f1+1;   // Just to have a sure positive value
total.f[1]=f2;
   break;

      case 35:   // Multiobjective F1
 // x1 in [0,1]
 //   x2 in [0,1]

    x1=pos.p.x[ 0];
    x2=pos.p.x[ 1];
    f1=(x1*x1+x2*x2)/2;
    x1=x1-2;
    x2=x2-2;
    f2=(x1*x1+x2*x2)/2;

 total.f[0]=f1;
total.f[1]=f2;
   break;

         case 36:   // Multiobjective F2
 // x1 in [0,1]
 //   x2 in [0,1]


    x1=pos.p.x[ 0];
    x2=pos.p.x[ 1];
    f1=x1;
    x=1+9*x2;
    f2= x-sqrt(f1*x);

 total.f[0]=f1;
total.f[1]=f2;
   break;
   
    case 37:   // Multiobjective F3
 // x1 in (0,1]
 //   x2 in (0,1]

    f1=pos.p.x[ 0];
    x2=pos.p.x[ 1];
    x=1+9*x2;
    f2=x*(1-f1*f1/(x*x));

    
 total.f[0]=f1;   
total.f[1]=f2;
   break;

     case 38:   // Multiobjective F4
 // x1 in (0,1]
 //   x2 in (0,1]

    f1=pos.p.x[ 0];
    x2=pos.p.x[ 1];
    x=1+9*x2;
    x3=f1/x;
    f2=x*(1-pow(x3,0.25)-pow(x3,4));

 total.f[0]=f1;
total.f[1]=f2;
   break;

         case 39:   // Multiobjective F5
 // x1 in (0,1]
 //   x2 in (0,1]

    f1=pos.p.x[ 0];
    x2=pos.p.x[ 1];
    x=1+9*x2;
    x3=f1/x;
    f2=x*(1-pow(x3,0.5)-x3*sin(5*two_pi*f1));

 total.f[0]=f1;
total.f[1]=f2; 
   break;

   case 40:   // Multiobjective F6 (Deb)
 // x1 in (0,1]
 //   x2 in (0,1]

    f1=pos.p.x[ 0];
    x2=pos.p.x[ 1];  // For 2D
 //x2=0;  For 1D
    x=1+10*x2;
    x3=f1/x;

    f2=x*(1-x3*x3-x3*sin(4*two_pi*f1));

 total.f[0]=f1;
total.f[1]=f2;
   break;
   
case 41:
// Multiobjective   Coello F3
    f1=pos.p.x[ 0];
    x2=pos.p.x[ 1];

    x=11+x2*x2-10*cos(two_pi*x2);
    if (x>=f1)
    {
      f2=x*(1-sqrt(f1/x));
    }
    else
       f2=0;


 total.f[0]=f1;
total.f[1]=f2;
 break;
      
 case 42: // For implicit function
    problem[1].target=0;
    problem[1].eps=0.1;
    problem[1].printlevel=0;
    problem[1].funct=42;
	problem[1].Max_Eval=problem[0].Max_Eval;
    coeff=pos;
   post=PSO(1,problem[1].Max_Eval);
   total.f[0]=post.f.f[0];


 break;

 case 43: // Part 2 for implicit function . Here a spherical one
// coeff is a global variable

    DD= coeff.p.size;
    total.f[0]=0;
  for (d=0;d<DD;d++)
  {
     x=coeff.p.x[d];
     total.f[0]=total.f[0]+x*x;
  }
  total.f[0]=total.f[0]-2+ pos.p.x[ 0];
   break;

case 44:   // Jeannet_Messine
/*
  Ref. ROADEF 2003, p 273
  i ={0,1,2,3,4,5} (see static data a1 and a2)
  x2 in [-15,25]
  x3 in [3,10]
  min -112.5   in position (3,-7.5, 10)
  The authors define 5 different methods, and
their best result is: 3271 evaluations
TRIBES is _far_ better
*/

 i=pos.p.x[ 0];
x2=pos.p.x[ 1];

x3= pos.p.x[ 2];

total.f[0]=20*a1[i]*x2*x2 + 2*a2[i]*x2*x3;
break;

case 45: // MINLP X. Yan
total.f[0]=MINLP(pos,1); // Constraints
total.f[1]=MINLP(pos,0);// Function to minimize
break;

case 46: // Tripod function (Louis Gacogne)
// on [-100, 100], min 0
x1=pos.p.x[0];
x2= pos.p.x[1];

if(x2<0)
{
	total.f[0]=fabs(x1)+fabs(x2+50);	
}
else
{
	if(x1<0)
		total.f[0]=1+fabs(x1+50)+fabs(x2-50);
	else
		total.f[0]=2+fabs(x1-50)+fabs(x2-50);
}
 break;
 
 case 47: // Tripod function (Louis Gacogne)
// on [-100, 100], min 0
// using hierarchical multiobjective optimization
// constraint= in circle radius 50, center ( 0,-50)
n1=1;
n2=0;
x1=pos.p.x[0];
x2= pos.p.x[1];


 r=    x1*x1+(x2+50)*(x2+50);
 r=sqrt(r);

total.f[n1]=constrain_positive(50-r);
   //printf("\n%f, %f, %f,  %f",x1,x2, r,total.f[0]);
if(x2<0)
{                           
	total.f[n2]=fabs(x1)+fabs(x2+50);
}
else
{
	if(x1<0)
		total.f[n2]=1+fabs(x1+50)+fabs(x2-50);
	else
		total.f[n2]=2+fabs(x1-50)+fabs(x2-50);
}
break;

case 48: // DeJong f4
	//      Global min f(x)=0, x(i)=0

	total.f[0] = 0;
	for( d=0;d<DD;d++)
	{
		total.f[0] = total.f[0] + pos.p.x[ d]*pos.p.x[ d]*pos.p.x[ d]*pos.p.x[ d];
	}
	break;

 case 49: // Pressure vessel 
	 	  // Ref New Optim. Tech. in Eng. p 638
/* D=4
  1.1 <= x1 <= 12.5        granularity 0.0625
  0.6 <= x2 <= 12.5         granularity 0.0625
  0 .0 <= x3 <= 240
  0.0 <= x4 <= 240
  constraints
  g1:= 0.0193*x3-x1 <=0
  g2 := 0.00954*x3-x2<=0
  g3:= 750*1728-pi*x3*x3*(x4-(4/3)*x3)  <=0
*/           
	x1=pos.p.x[0];
	x2= pos.p.x[1];
	x3=pos.p.x[2];
	x4= pos.p.x[3];
 
   f=0.6224*x1*x3*x4 + 1.7781*x2*x3*x3 + 3.1611*x1*x1*x4 + 19.84*x1*x1*x3;
 
  // Constraints
  y=0.0193*x3-x1; if ( y>0) {c= 1+pow(10,10)*y;  f=f*c*c; }
  y=  0.00954*x3-x2; if (y>0) {c=1+y; f=f*c*c;  }
  y = 750*1728-pi*x3*x3*(x4+(4.0/3)*x3);  if (y>0) {c=1+y; f=f*c*c;  }
  
  total.f[0]=f;
 break;  
 
 case 50: // Coil compression spring 
	 	  // Ref New Optim. Tech. in Eng. p 644

	x1=pos.p.x[0];

	x2= pos.p.x[1];
	x3= pos.p.x[2];

	f=pi*pi*x2*x3*x3*(x1+2)*0.25;
	
	// Constraints
	Cf=1+0.75*x3/(x2-x3) + 0.615*x3/x2;
	K=0.125*G*pow(x3,4)/(x1*x2*x2*x2);
	sp=Fp/K;
	lf=Fmax/K + 1.05*(x1+2)*x3;

	y=8*Cf*Fmax*x2/(pi*x3*x3*x3) -S;
	if (y>0) {c=1+y; f=f*c*c*c;}

	y=lf-lmax;
	if (y>0) {c=1+y; f=f*c*c*c;}

	y=sp-spm;
	if (y>0) {c=1+y; f=f*c*c*c;}

	//y=sp+(Fmax-Fp)/K + 1.05*(x1+2)*x3 - lf;
	y=sp-Fp/K;

	if (y>0) {c=1+pow(10,10)*y; f=f*c*c*c;}

	y=sw- (Fmax-Fp)/K;
	if (y>0) {c=1+pow(10,10)*y; f=f*c*c*c;}

  total.f[0]=f;

	 break;


 case 51: // Gear train
	  // Ref New Optim. Tech. in Eng. p 634

	x1=pos.p.x[0];
	x2= pos.p.x[1];
	x3=pos.p.x[2];
	x4= pos.p.x[3];

	f=1/6.931 - x1*x2/(x3*x4); f=f*f;
	total.f[0]=f;
	break;

 case 52: // Pressure vessel as multiobjective (4 functions)
             
	x1=pos.p.x[0];
	x2= pos.p.x[1];
	x3=pos.p.x[2];
	x4= pos.p.x[3];
  
	 total.f[0]=0.6224*x1*x3*x4 + 1.7781*x2*x3*x3 + 3.1611*x1*x1*x4 + 19.84*x1*x1*x3;
 
  // Constraints

	y=0.0193*x3-x1; 
	total.f[1]=y+fabs(y);



	y=  0.00954*x3-x2;
	total.f[2]=y+fabs(y);

	y = 750*1728-pi*x3*x3*(x4+(4.0/3)*x3); 
	total.f[3]=y+fabs(y);  

 break;

 case 53: // Coil compression spring as multiobjective (6 functions)

	x1=pos.p.x[0];
	x2= pos.p.x[1];
	x3= pos.p.x[2];

	 total.f[0]=pi*pi*x2*x3*x3*(x1+2)*0.25;
	
	// Constraints
	Cf=1+0.75*x3/(x2-x3) + 0.615*x3/x2;
	K=0.125*G*pow(x3,4)/(x1*x2*x2*x2);
	sp=Fp/K;
	lf=Fmax/K + 1.05*(x1+2)*x3;

	y=8*Cf*Fmax*x2/(pi*x3*x3*x3) -S;
	total.f[1]=y+fabs(y);

	y=lf-lmax;
	total.f[2]=y+fabs(y);

	y=sp-spm;
	total.f[3]=y+fabs(y);

	y=sp-Fp/K;
	total.f[4]=y+fabs(y);

	y=sw- (Fmax-Fp)/K;
	total.f[5]=y+fabs(y);

	 break;

 case 54: // Pressure vessel by homeomorphism
// WARNING. TO COMPLETE
// May be not possible
	alpha=0.0193;
	beta=0.00954;
	gamma=750*1728;
	delta=240;

	y1=pos.p.x[0];
	y2=pos.p.x[1];
	y3=pos.p.x[2];
	y4=pos.p.x[3];

	x4=y4+delta;

	a_1=0.75*x4;
	b=0.75*(y3-gamma)/pi;

	z=solve_polyn_3(b,0,a_1,1);

	y=z.z[0]; choice=1;

root_choice:

	x3=y-a_1/3;
	x1=alpha*x3-y1;
	x2=beta*x3-y2;
printf("\nx1,x2,x3 %f %f %f",x1,x2,x3);


	if (x1<0 || x2<0)// || x3<0)
	{
		if(z.status>1)
		{
			if (choice==1) {choice=2; y=z.z[1]; goto root_choice;}
			if(choice==2) {choice=3;y=z.z[2]; goto root_choice;}
		}
		printf("\n error my_function case 54. No possible root");
		scanf("%i",&d);
	}


//	x=x3*x3*x3+a_1*x3*x3+b;
//	printf("\n x %f =? 0",x); // Should be equal to 0



   f=0.6224*x1*x3*x4 + 1.7781*x2*x3*x3 + 3.1611*x1*x1*x4 + 19.84*x1*x1*x3;

//	y=gamma-pi*x3*x3*(x4+(4.0/3)*x3); // Should be equal to y3;
//	printf("\nx3 %.12f, x4 %f,  %f =? %f",x3,x4,y,y3);

	printf("\n y %f %f %f %f",y1,y2,y3,y4);
	printf("\n x %f %f %f %f",x1,x2,x3,x4);
	printf("\n f %f",f);

	total.f[0]=f;
	 break;


case 55: // Neural Network Training. 2 Bit Parity (XOR)
	// Dim=9
	// min 0.0
      total.f[0]= ANNXor(pos.p);
break;
      
case 56: // Neural Network Training. 4 Bit Parity
	// Dim=25
	// A good result is 0.11
      total.f[0]= ANNParity4(pos.p);
   break;
case 57: // Neural Network Training. Three Color Cube
	// Dim=46
	// A good result is 0.2
	total.f[0]=ANNCOLORCUBE(pos.p);
	break;
case 58: // Neural Network Training. Diabetes in Pima Indians
	// Dim=64
	// A good result is 0.16
	total.f[0]=ANNPIMA(pos.p);
	break;

case 59: // Neural Network Training. Sin Times Sin
	// Dim=26
	// A good result is 0.23
	total.f[0]= ANNSINSIMP(pos.p);
	break;
case 60: // Neural Network Training. Rise Time Servomechanism
	// Dim=28
	// A good result is 0.45
	total.f[0]= ANNSERVO(pos.p);
	break;
case 61: // Moving Peaks 
	// Warning. Only for mono-objective non recursive optimization
	// (rank 0 in .f and level=0)
	for (d=0;d<DD;d++) gen[d]=pos.p.x[d];
	total.f[0]=eval_movpeaks (gen);
  //printf("\n %f %f",global_max,  total.f[0]);


	if(recent_change==1 && recurs==0)
	{
//printf("\n %f", best_result.f.f[0]);		
		recurs=1; // myfunction will be called again in reinit_swarm
				// so, it is will be then necessary to skip this part
		n_change=n_change+1;
	// Recompute the memorized f values

	reinit_swarm(level,0); // Re-initialise the swarm but keeping the best result
//	reinit_swarm(level,1); // Re-initialise the swarm without keeping the best result

	offline_error+=total_error(best_result.f);
	offline_error_cont+=total_error(best_result.f);
	recent_change=0;
	recurs=0;

	// The current particule has been evaluated twice
	eval_f_tot[level]=eval_f_tot[level]-1; 
	}	 
	break;

case 62: // Goldberg’s order 3 deceptive problem
	// Be sure dimension=3*k
	// Interval inside ]-0.5 1.5[
	// and granularity=1, so that the only possible x values
	// are 0 and 1
	// The min value is 0
f=0;
/* // Exact binary form
for (d=0;d<DD-2;d=d+3)

{
	bit_sum=0;
	for (i=d;i<d+3;i++) bit_sum+=pos.p.x[i];
	if (bit_sum==0) {f+=0.9; continue;}
	if (bit_sum==1) {f+=0.6; continue;}
	if (bit_sum==2) {f+=0.3; continue;}
	if (bit_sum==3) {f+=1; continue;}
}
*/
/* // Continuous form
for (d=0;d<DD-2;d=d+3)
{
	s=0;
	for (i=d;i<d+3;i++) s+=pos.p.x[i];
	s=MIN(s,3); s=MAX(s,0);
	if(s<=2) f+=0.9-0.3*s;
	else f+=0.3+0.7*(s-2);
}
*/

// Continuous integer form
for (d=0;d<DD-2;d=d+3)
{
	s=0;
	for (i=d;i<d+3;i++) s+=floor(pos.p.x[i]+0.5);
	s=MIN(s,3); s=MAX(s,0);
	if(s<=2) f+=0.9-0.3*s;
	else f+=0.3+0.7*(s-2);
}

total.f[0]=DD/3-f;
break;

case 63: // Mühlenbein's order 5 (MC variant), seen as a 1D function
        // Binary optim. Min 0 on 1057 = 100010001000
D=15; // Length of the bit string
f=0;
bit_x=bits(pos.p.x[0],D);
bit.size=5;
for (d=0;d<D-3;d=d+5)
{
	for (i=d;i<d+5;i++) bit.x[i-d]=bit_x.x[i];
   num=number(0,bit);

   if (num==0) {f+=3.0; continue;}
   if (num==1) {f+=4.0; continue;}
   if (num==3) {f+=2.0; continue;}
   if (num==7) {f+=1.0; continue;}
   if (num==31) {f+=3.5; continue;}   
}

total.f[0]=4*D/5-f;
	break;

 case 64: // Schaffer'f6  (2D)
 	x1=pos.p.x[0];
	x2= pos.p.x[1];
  x1=x1-1; x2=x2-1;    // Optional shift
  r = x1*x1 + x2*x2;
  y=  1 + 0.001*r ;
  y1 = y*y;
  y=sin(sqrt(r));
  y2 = y*y - 0.5;
    total.f[0] = 0.5 + y2/y1;
 break;

//-----------------------------------------------------------
case 99: //  Test functions
// Test
	x1=pos.p.x[0];
	x2= pos.p.x[1];
	f1=2*x1*x2*(1-x2)+ x2*x2;
	f2=2*(1-x1)*x2*(1-x2)+(1-x2)*(1-x2);
	total.f[0]=1-f1;
	break;
// polynomial fitting problem 
		// Cf. Differential Evolution Homepage, C code
		// on [-300 200]^9

	f=0; y=-1;
   dx = 2/dx;
   for (i=0;i<=M;i++)
   {
      py = pos.p.x[0];
      for (j=1;j<DD;j++)
      {
		py = y*py + pos.p.x[j];
      }
      if (py<-1 || py>1) f+=(1-py)*(1-py);
      y+=dx;
   }
   py = pos.p.x[0];
   for (j=1;j<DD;j++) py=1.2*py+pos.p.x[j];
   py = py-72.661;
   if (py<0) f+=py*py;
   py = pos.p.x[0];
   for (j=1;j<DD;j++) py=-1.2*py+pos.p.x[j];
   py =py-72.661;
   if (py<0) f+=py*py;
		total.f[0]=f;
   break;
// Stochastic
	total.f[0] = 0;
	for( d=0;d<DD;d++)
	{
		xd=pos.p.x[d];
		total.f[0]=	total.f[0]+(d+1)*pow(xd,4)+alea_normal(0,1);
	}
	break;	

// Step
	total.f[0] = 0;
	for( d=0;d<DD;d++)
	{
		total.f[0]=	total.f[0]+(int)pos.p.x[d];
	}
	break;

// Beatle (Luc Heylen)
 // Min 0 on (1,1)
 // search space [-3 5]^2 => succes rate 76%
	x1=pos.p.x[0];
	x2= pos.p.x[1];
f = 100*(x2-x1*x1)*(x2-x1*x1)+(1-x1)*(1-x1);
f2=100*(-x2+2-x1*x1)*(-x2+2-x1*x1)+(1-x1)*(1-x1);
//f=log(0.2+f+10*(f*(5.1-mod(sqrt(f+f2),5))))+1.61;
y1=sqrt(f+f2);
y2=y1-5*((int)y1/5);
y3= 0.2+f+10*(f*(5.1-y2));
//f=log(y3)-log(0.2); // Modified function (just one solution)
 f=log(y3)+1.61; // Original function
  total.f[0]=f;
break;
 // Ripple (Luc Heylen)
// Min 0 on (1,1)
// search space [-3 5]^2  => success rate 64%
	x1=pos.p.x[0];
	x2= pos.p.x[1];

y1=1-exp(-3*(x1-1)*(x1-1)-0.3*(x2-1)*(x2-1));
y2 = 100*(x2-x1*x1)*(x2-x1*x1)+(1-x1)*(1-x1);
f=(12+sqrt(y2)/2-10*cos(20*y1)-abs((cos(3*(x1-1))+cos(10*(x2-1)))))/10;
total.f[0]=f;
break;
 

// Dynamic Parabola (Sphere)
  sigma=0.01;
  pr=0.5;
  if (alea(0,1)<pr) r=0; else r=1;
	total.f[0] = 0;
	for( d=0;d<DD;d++)
	{
     //xd= pos.p.x[ d]*(1+r*alea_normal(0,sigma));
     xd= pos.p.x[ d] + r*alea_normal(0,sigma); 
    total.f[0] = total.f[0] +xd*xd ;
	}
  DYN=1;
	break;
// Leon 2D
// on [-1.2 1.2]^2 min 0 on (1 1)
	x1=pos.p.x[0];
	x2= pos.p.x[1];
	f1=x2-x1*x1*x1;
	 total.f[0]=100*f1*f1 +(1-x1)*(1-x1);
	 break;

// Schaffer 2D
// on [-100 100]^2  min 0 on (0 0)
	x1=pos.p.x[0];
	x2= pos.p.x[1];
	f1=1+0.001*(x1*x1+x2*x2);
	f1=f1*f1;
	 total.f[0]=0.5+(sin(sqrt(x1*x1+x2*x2))-0.5)/f1;
	 break;

// Bukin6 2D
// on [-15 -5]x[-3 3] min 0 on (-10 1)
	x1=pos.p.x[0];
	x2= pos.p.x[1];
	 total.f[0]=100*sqrt(abs(x2-0.01*x1*x1)) +0.01*abs(x1+10);
	 break;

// Schwefel 2D
// on [-500 500]^2, min=-837.9658 on (420.9687, 420.9687)
	x1=pos.p.x[0];
	x2= pos.p.x[1];
	 total.f[0]=-x1*sin(sqrt(abs(x1)))-x2*sin(sqrt(abs(x2)));
break;


  // somme i*x_i
  total.f[0] = 0;
	for( d=0;d<DD;d++)
	{
		total.f[0] = total.f[0] + (d+1)*pos.p.x[ d];
	 }

   break;

	 // Livre sur l'OEP. Chapitre Contraintes
	// homéomorphisme disque/4=carré
    y1=pos.p.x[ 0]; // in [0,1]
    y2=pos.p.x[ 1]; // in [0,1]

	x1=y1*cos(pi*y2/2); 
	x2=y1*sin(pi*y2/2);

	f1=(x1-1)*(x1-1)+(x2-1)*(x2-1); 

	total.f[0]=f1;

 break;
	

	 // Livre sur l'OEP. Chapitre Contraintes
	// Méthode multicritère
  // Valeur min 0.17157
    x1=pos.p.x[ 0]; // in [0,2]
    x2=pos.p.x[ 1]; // in [0,2]

	f1=(x1-1)*(x1-1)+(x2-1)*(x2-1); // Parabola
	f2=x1*x1+x2*x2 -1; // disc
	f2=f2+fabs(f2); // inside the disc

 total.f[0]=f1;
 total.f[1]=f2;
 break;

	 // Livre sur l'OEP. Chapitre Contraintes
	// Méthode monocritère
    x1=pos.p.x[ 0]; // in [0,2]
    x2=pos.p.x[ 1]; // in [0,2]

	f1=(x1-1)*(x1-1)+(x2-1)*(x2-1); // Parabola
	f2=x1*x1+x2*x2 -1; // disc
	f2=f2+fabs(f2); // inside the disc

 total.f[0]=f1+ f2;

 break;


	// Schewefel's function
	total.f[0] = 42000;
	for( d=0;d<DD;d++)
	{
		xid=pos.p.x[d];
		total.f[0] = total.f[0] -xid*sin(sqrt(fabs(xid)));
	}
	break;


	// Ellipsoid     Global min f(x)=0, x(i)=0

	total.f[0] = 0;
	for( d=0;d<DD;d++)
	{
		total.f[0] = total.f[0] + (d+1)*pos.p.x[ d]*pos.p.x[ d];
	}
	break;

 //g03  on [0 1]^D.  min 0 on (sqrt(D)... sqrt(D))
 //A Constraint-Handling Mechanism for Particle Swarm Optimization
 //Gregorio Toscano Pulido and Carlos A. Coello Coello
 // CEC 2004
 f1=1; for( d=0;d<DD;d++) f1=f1* pos.p.x[ d];
    f1=f1*sqrt((double)DD);

 f2=0;  for( d=0;d<DD;d++) f2=f2+pos.p.x[ d]*pos.p.x[ d];

            total.f[0]  =1- f1;
            total.f[1]  = fabs(f2-1);      // Constraint

               total.f[0]=  total.f[0]    +           total.f[1] ;   
break;
 // on [-10 10]^2. min 80.70404  on (10, 10)
 	f1=0; f2=1;
	for( d=0;d<DD;d++)
	{
		f1 = f1+ pos.p.x[ d];
            f2=f2*cos(  pos.p.x[ d]);
	}
           total.f[0]  = 100-(f1-f2);
           break;
   
   x1=pos.p.x[ 0];
   f1=sin(5*pi*(pow(x1,0.75)-0.05));
      total.f[0]  =pow(f1,6);
      break;


  x1=pos.p.x[ 0];
  f1=sin(5*pi*x1);
  total.f[0]  =pow(f1,6);

//printf("\n %f %f %f   ", x1,f1,total.f[0]); 
break;

x1=pos.p.x[ 0];
f1=fabs(sin(x1));
f2=fabs(x1);
total.f[0]=f1;
total.f[1]=f2;
break;
// Levy n° 3 on [-10,10]^2. Min=-176.542
// 760 local minima, 18 global ones.
// About 2400 evaluations with strategies (12 12 12 1) and fuzziness option on
// but only 850 with (21 21 21 21) and fuzziness off

// with success rate= 100%
// Best result by Parsopoulos (2002) with Stretching PSO: 7000

x1=pos.p.x[ 0];
x2=pos.p.x[ 1];

f1=0;
f2=0;
for (i=1;i<=5;i++)
{
	x=(double)i;
	f1=f1+x*cos((x-1)*x1+x);
	f2=f2+x*cos((x+1)*x2+x);
}

total.f[0]=f1*f2;
break;

// Integer programming. Rudolph F1
// [-100, 100]^D, min 0
// Use granularity=1
// D=5 strategy (12 12 12 1)=> 360 eval. Best result Parsopoulos (2002) 8800
// D=25........................7500 ...............................160000)

f1=0;

for (d=0;d<DD;d++)
{
	f1=f1+fabs(pos.p.x[d]);
}

total.f[0]=f1;
break;

// XOR function (Blum 1989)(Vrahatis 2000)(Parsopoulos 2002)
// on [-1,1]^9
x1=pos.p.x[ 0];
x2=pos.p.x[ 1];

x3=pos.p.x[ 2];
x4=pos.p.x[ 3];
x5=pos.p.x[ 4];
x6=pos.p.x[ 5];
x7=pos.p.x[ 6];
x8=pos.p.x[ 7];
x9=pos.p.x[ 8];


f1=pow(1+exp(-x7/(1+exp(-x1-x2-x5))-x8/(1+exp(-x3-x4-x6))-x9),-2);
f1=f1+pow(1+exp(-x9-x7/(1+exp(-x5))-x8/(1+exp(-x6))),-2);
f1=f1+pow(1-1/(1+exp(-x9-x7/(1+exp(-x1-x5))-x8/(1+exp(-x3-x6)))),2);
f1=f1+pow(1-1/(1+exp(-x9-x7/(1+exp(-x2-x5))-x8/(1+exp(-x4-x6)))),2);

total.f[0]=f1;
break;


// Freudenstein-Roth  on [-5,5]^2. Min=0. About 3000 evaluations
x1=pos.p.x[ 0];
x2=pos.p.x[ 1];

total.f[0]=pow(-13+x1+((5-x2)*x2-2)*x2,2) +pow(-29+x1+((x2+1)*x2-14)*x2,2);
break;


x1=pos.p.x[ 0];
x2=pos.p.x[ 1];

 total.f[0]=cos(x1)*cos(x1)+sin(x2)*sin(x2);
 break;


  // Circle

x1=pos.p.x[ 0];
x2=pos.p.x[ 1];

 total.f[0]=fabs(x1*x1+x2*x2-1);
break;

 
 // Multiobjective
x1=pos.p.x[ 0];
x2=pos.p.x[ 1];
f1=x1;
f2=1-x1*x1-x2*x2;

total.f[0]=f1;
total.f[1]=f2;
break;

 /*
     min distance(x-A)
     A=(2,0)
     with x1*x1 + x2*x2 -1=0
*/
  	total.f[0] = -1;
	for( d=0;d<DD;d++)
	{
		total.f[0] = total.f[0] + pos.p.x[ d]*pos.p.x[ d];
	 }

      total.f[0]=fabs(total.f[0]);
      
      x=sqrt((2-pos.p.x[ 0])* (2-pos.p.x[ 0])+(0-pos.p.x[ 1])* (0-pos.p.x[ 1]) ); 
    total.f[0]=total.f[0]+ x;
    break;



/*

     min xD
     with x1*x1 + x2*x2 + x3*x3 + ...xD*xD-1=0
*/
  	total.f[0] = -1;
	for( d=0;d<DD;d++)
	{
		total.f[0] = total.f[0] + pos.p.x[ d]*pos.p.x[ d];
	 }
 
      total.f[0]=fabs(total.f[0]);
    x=problem[level].H.min[DD-1]- pos.p.x[ DD-1];

    total.f[0]=total.f[0]+ fabs(x);
  break;


  
default:
	printf("\n ERROR, unknown objective function %i",funct);
	total.f[0]=NA;
break;
}

total.f[0]=fabs(target-total.f[0]); // Warning. For  a multiobjective problem
									// "target" is only for the first function

if(problem[level].printlevel>2)
   printf("\n my_function => %f",total.f[0]);

return total;
}

// ----------------------------------------------------------------------------- TOT_FIFTY_FIFTY
double tot_fifty_fifty(struct position pos)

// Sum of first numbers = sum of last numbers

{
int		a;
int		d;
int		d2;
double	tax_twice;
double	tot_f;
double	xa,xb;

a=(int)(0.5*(pos.p.size+1));

xa=0;

for (d=0;d<a;d++) xa=xa+pos.p.x[d];

xb=0;
for (d=a;d<pos.p.size;d++) xb=xb+pos.p.x[d];

tot_f=(float)fabs(xa-xb);
return tot_f;

}

