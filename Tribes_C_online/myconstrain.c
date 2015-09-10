
// ----------------------------------------------------------------------------- ALL_DIFFERENT
struct position	all_different(struct position pos, float eps,int level )
{
/*
 Modify the position so that all components are different
 Should be used only when granularity>O
 Useful for some combinatorial problems, like Knapsack or Fifty-fifty
*/

int				d,d1;
float			dx;
int				k;
float			max_loop;
float			nloop;
struct position	post;


post=pos;

for (d=1;d<pos.p.size;d++) // For each component ...
{
	if (problem[level].H.dx[d]>0) {dx=problem[level].H.dx[d];} else {dx=eps;}
	max_loop=1+(problem[level].H.max[d]-problem[level].H.min[d])/dx;
	k=1;
	nloop=0;
	loop1:
	//... if component # of all previous, OK, ...
	for (d1=0;d1<d;d1++)
	{
		if (post.p.x[d]==post.p.x[d1]) // ... else, modify slightly and retry
		{
			loop2:
			post.p.x[d]=pos.p.x[d]-k*dx;
			if (k<0) {k=-k+1;} else {k=-k;}
			if (post.p.x[d]<problem[level].H.min[d]) goto loop2;
			if (post.p.x[d]>problem[level].H.max[d]) goto loop2;
				nloop=nloop+1;
				if (nloop==max_loop)
				{
					printf("\n WARNING. Some components are equal (%i and %i) Impossible to modify this particule:",d1,d);
					display_position(post);
					return post;
				}
				goto loop1;
		}
	}
}
return post;
}


// ----------------------------------------------------------------------------- CONSTRAIN
struct particle constrain(struct particle part,int level)
{ /*
 Keep the particle in the definition domain
*/
int				cmax;
int				d;
double			dd,dd1;
int				discr,discr_v;
int				n_discrete;
struct particle partt;
double			pr;
double			x_discr;			
double			x;


cmax=(int)problem[level].max_fr;
partt=part;

if (problem[level].printlevel>2)
	printf("\n constrain.  label %i",partt.label);

   // Keep it in the box (hyperparallelepid)
if(confin_interv<0) goto granul;
//printf("\n ici %f",partt.x.p.x[0]);

for (d=0;d<partt.x.p.size;d++)
{
	if( partt.x.p.x[d]>problem[level].H.max[d])
	{
         partt.x.p.x[d]=problem[level].H.max[d];
//printf(" => %f",problem[level].H.max[d]);
         continue;
	}

	if( partt.x.p.x[d]<problem[level].H.min[d])
	{
            partt.x.p.x[d]=problem[level].H.min[d];
	}
}

granul:
  if (problem[level].printlevel>2)
	  printf("\n Granularity taken into account");
   for (d=0;d<partt.x.p.size;d++)
   {
          if(problem[level].H.dx[d]>0)

             partt.x.p.x[d]=regranul( partt.x.p.x[d],abs(problem[level].H.dx[d]));
   }


all_diff:

if (problem[level].all_diff>0) // Eventually, modify so that all components are different
		partt.x=all_different(partt.x,problem[level].eps,level);


 //display_position(partt.x);

 // Discrete values for some special functions
if (problem[level].constrain==1)
{
	for (discr=0;discr<discrete_nb;discr++) // For each special discrete variable
	{
		discr_v=discrete[discr].d; // Rank of the discrete variable
		n_discrete=discrete[discr].size; // Number of acceptable values
		x=partt.x.p.x[discr_v];

 // Find the nearest discrete value
		x_discr=discrete[discr].v[n_discrete-1];

		if (x>=x_discr) { partt.x.p.x[discr_v]=x_discr; break;}

		for (d=0;d<n_discrete-1;d++)
		{
			dd=discrete[discr].v[d]; dd1=discrete[discr].v[d+1];
			if(x<dd || x>=dd1) continue;
	 // dd<= x < dd1
			if (x-dd<dd1-x) 
			{ partt.x.p.x[discr_v]=dd; goto end1;}
			else
			{ partt.x.p.x[discr_v]=dd1; goto end1;}
		}
printf("\n ERROR constrain.c. I find no discrete value for %f",x); 
end1:;
	}
}

return partt;
}



//---------------------------------------------------------------------------------------------------------- CONSTRAIN_NULL
double   constrain_null(double A)
{
  // Check "how much" A is not null

return fabs(A);
}
 //---------------------------------------------------------------------------------------------------------- CONSTRAIN_POSITIVE
double   constrain_positive(double A)
{
  // Check "how much" A is not positive

   if( A<0) return -A;
   return 0;
}

//---------------------------------------------------------------------------------------------------------- DICHOTOM
  struct position dichotom(struct position p0, struct position p1, int level)
  {
   /*
    The constraint is defined by g(position)>=0
    We suppose we have funct(p0)>=0 and funct(p1)<0

  */
    struct position a0,a1,b;
    int  d,D;
    double  eps;
    int  funct;
    double g0,g1,g;

    a0=p0;
    a1=p1;
    D=p0.p.size;

     eps= problem[level].target;
	 funct=problem[level].funct;

 dicho:
    g0=Myconstrain(a0,level);
    g1=Myconstrain(a1,level);
    if(g0-g1<=eps) {b=a0;goto end;}

    //b=(a0+a1)/2
     for (d=0;d<D;d++)
     {
        b.p.x[d]=(a0.p.x[d]+a1.p.x[d])/2;
     }

      g=Myconstrain(b,level);
       if (g==0) goto end;
          if (g<0) a1=b;  else a0=b;
          goto dicho;

  end:
     // Evaluate position
      b.f=MyFunction(b,funct,eps,level);       //evaluation of the position

     return b;
  }

  // ----------------------------------------------------------------------------- DOMAIN_POS
struct position domain_pos(struct position pos,struct search_space H)
{
// Modify the position according to the "granularity" of the search space
int				d;
struct position post;

post=pos;
for (d=0;d<pos.p.size;d++)
{
	post.p.x[d]=regranul(pos.p.x[d],H.dx[d]);
}

return post;
}


  //---------------------------------------------------------------------------------------------------------- MY_CONSTRAIN 
double   Myconstrain(struct position pos, int level)
//
//WORK IN PROGRESS
//A domain is defined by a function dom(x)<=0
 //    Return 1 if the position is inside the domain

{
int                     constraint;
int                     d;
int                     DD;
double               dom;
int                     type;
 double              r;

 DD=  pos.p.size;
 
 constraint=0;     // TEST
 
 switch (constraint)
 {
   case 0:  // Sphere    center 0, radius 1
      r=1; // Radius

     dom=0;
   for (d=0;d<DD;d++)
   {
         dom=dom+pos.p.x[d]*pos.p.x[d] ;
   }
    dom=dom-r*r;
   return dom;

    
     case 99: // test
     
// TO COMPLETE
    
   
    break;
 }
return dom;
}


                  
