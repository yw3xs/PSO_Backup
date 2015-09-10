double MINLP(struct position pos, int option)
{

/*
Problem submitted by X. Yan (2003-09)
  pos dimension is 12*3*2=72
  The first 1 to 36 are called X, for X[3][12]. X is binary
  the last  37 to 72 are called P, for P[3][12]

  Here the data are hard coded, just for test, but, of course,
  in real problems like this onen, they should be read on a file.

  NOTE. The code is not optimized. Some values like Xit, Pit etc. are computed
  several times. This on purpose for the moment, in order to have something clean
  and easy to understand.
  

*/
//static int  N={3};
//static int T={12};
int N=3;
int T=12;
static int shift={36}; // N*T

static double D[12] =
{
170,250,400,520,700,1050,1100,800,650,330,400,550
} ;

static  double SR[12]=
{
20,25,40,55,70,95,100,80,65,35,40,55
};

static double Pmax[3]={600,400,200};
static double Pmin[3]={100,100,50};
static double a[3]={500,300,100};
static double b[3]={10,8,6};
static double c[3]={0.002,0.0025,0.005};
static int Tup[3]={3,3,3};
static int Tdown[3]={3,3,3};
static double ST[3]={450,400,300};
static double coef[3]={0.8,0.8,0.8};
static int X_init[3]={0,1,1};
//static int Yc_init[3]={0,3,3};
//static int Yd_init[3]={3,0,0};
 
double   constr_1[12];
double   constr_3[12];
double   constr_4a[3][12], constr_4b[3][12];
double   constr_6[3][12];
double   constr_7[3][12];
double f;
double   Fit;
int         h;
int      i;
double   Pit,Pit_1;
int         rank;
int      t;
double   TC;
double   Xih,Xit,Xit_1;
//double   Yc[3][12],Yd[3][12];
double   Z;

int   constr[8]={0};
// For test purpose, you can activate or not any part of the composite objective function
constr[0]=1;

constr[1]=1;
constr[3]=1;
constr[4]=1;
constr[6]=1;
constr[7]=0;

 TC=0;
//printf("\nMINLP. N=%i, T=%i, shift=%i",N,T,shift);
//display_position(pos);
 if(option==1) goto constr1; 
 
 // Main objective function


 if(constr[0]==0) goto constr1;
 for (i=0;i<N;i++)
 {
      for (t=1;t<T;t++)
      {
                rank=  i*T+t;
                Xit=pos.p.x[ rank];
                Xit_1=pos.p.x[rank-1];
                rank=rank+shift;
                Pit=pos.p.x[ rank];
                
                Fit=a[i]+b[i]* Pit+c[i]*Pit*Pit;
                TC=TC+(Fit+ST[i]*(1- Xit_1))*Xit;        
      }
       // Special case t=0
          rank=i*T;
          Xit=pos.p.x[ rank];
          Xit_1=X_init[i];
          Pit=pos.p.x[ rank];

          Fit=a[i]+b[i]* Pit+c[i]*Pit*Pit;
          TC=TC+(Fit+ST[i]*(1- Xit_1))*Xit;          
 }

 if (option==0) return TC;

 constr1:

 if (constr[1]==0) goto constr3;
 //printf("\n constraints 1");
 // Constraints 1
 for (t=0;t<T;t++)
 {
    Z=0;
    for (i=0;i<N;i++)
    {
       rank=  i*T+t;
       Xit=pos.p.x[ rank];
       rank=rank+shift;
       Pit=pos.p.x[ rank];
       Z=Z+Pit*Xit;
    }
    Z=Z-D[t];
    constr_1[t]=constrain_null(Z);
//printf("\n constr_1[t] %f",constr_1[t]);
 }

constr3:
if (constr[3]==0) goto constr4;
 //printf("\n Constraints 3");
 // Constraints 3
 
 for (t=0;t<T;t++)
 {
    Z=0;
    for (i=0;i<N;i++)
    {
       rank=  i*T+t;
       Xit=pos.p.x[ rank];
       rank=rank+shift;
       Pit=pos.p.x[ rank];
       Z=Z+Pmax[i]*Xit;
    }
    Z=Z-D[t]-SR[t];
    constr_3[t]=constrain_positive(Z);
 //printf("\n constr_3[t] %f",constr_3[t]);
 }

 constr4:
 if(constr[4]==0) goto constr6;
 //printf("\n Constraints 4");
 // Constraints 4
 
  for (t=0;t<T;t++)
 {
    for (i=0;i<N;i++)                                                                                                     
    {
    rank=  i*T+t;
       Xit=pos.p.x[ rank];
       rank=  rank+shift;
       Pit=pos.p.x[rank];
       Z=Pit-Xit*Pmin[i];
       constr_4a[i][t]=constrain_positive(Z);
 //if (constr_4a[i][t]>0) printf("\n rank %i Xit %f Pit %f Pmin %f c4a %f",   rank-shift,Xit, Pit,Pmin[i],constr_4a[i][t]);

       Z=Xit*Pmax[i]-Pit;
       constr_4b[i][t]=constrain_positive(Z);
  //if (constr_4b[i][t]>0) printf("\n rank %i Xit %f Pit %f Pmax %f c4b %f",   rank-shift,Xit,Pit,Pmax[i],constr_4b[i][t]);
    }     
 }

 constr6:
 if(constr[6]==0) goto constr7;
 //printf("\n Constraints 6");
 // Constraints 6
 //for (t=0;t<T;t++)
   for (t=1;t<T;t++)
 {
    for (i=0;i<N;i++)
    {
       rank=  i*T+t+shift;
       Pit=pos.p.x[ rank];
       Pit_1= pos.p.x[ rank-1];
       Z=coef[i]*Pmax[i]-fabs(Pit-Pit_1);
       constr_6[i][t]=constrain_positive(Z);
    }
  //printf("\n %f",   constr_6[t]);
 }

 constr7:
 if(constr[7]==0) goto combine;
 //printf("\n Constraints 7");
 // Constraints 7

 //for (t=0;t<T;t++)
  for (t=1;t<T;t++)
 {
    for (i=0;i<N;i++)
    {
 //printf("\n i,t %i %i",i,t);
       rank=  i*T+t;
       Xit=pos.p.x[ rank];
       Z=0;
        if (Xit==1)
        {
           if (t>=Tup[i]-1)
           {
              for(h=t-Tup[i]+1;h<=t-1;h++)
              {
                 rank=  i*T+h;
                 Xih=pos.p.x[ rank];
                 Z=Z+Xih;
              }
              Z=Tup[i]-Z;
           }
        }
        else
        {
           if(t>=Tdown[i]-1) 
           {
              for(h=t-Tdown[i]+1;h<=t-1;h++)
              {
                 rank=  i*T+h;
                 Xih=pos.p.x[ rank];
                 Z=Z+1-Xih;
              }
              Z=Tdown[i]-Z;
           }
        }
        constr_7[i][t]=constrain_positive(Z);
 //printf("\n i %i, t %i, %f",   i,t,constr_7[i][t]);
    }
 }

 combine:
 //printf("\n Combine constraints");
 // Combine constraints
 f=0;
if(constr[1]>0)
{
     for (t=0;t<T;t++)
      {
         f=f+constr_1[t];
      }
}
if(constr[3]>0)
{
 for (t=0;t<T;t++)
      {
         f=f+constr_3[t];
      }
}

if(constr[4]>0)
{
      for (i=0;i<N;i++)
      {
         for (t=0;t<T;t++)
         {
           f=f+constr_4a[i][t]+constr_4b[i][t];
         }
      }
 }

 if(constr[6]>0)
{
      for (i=0;i<N;i++)
      {
         for (t=1;t<T;t++)
         {
            f=f+constr_6[i][t];
         }
      }
 }

  if(constr[7]>0)
{
      for (i=0;i<N;i++)
      {
         for (t=1;t<T;t++)
         {
            f=f+constr_7[i][t];
         }
      }
 }
         
f=f+TC;
 
 return f;
}
