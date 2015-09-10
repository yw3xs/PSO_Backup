         // ----------------------------------------------------------------------------- DISPLAY_I_GROUP
void display_i_group(FILE *f_trace,struct particle part, int option)
{
int	i;

if (option==0)
{
printf("\n i_group of particle %i:   ",part.label);
for (i=0;i<part.ig.size;i++)
	printf("%i ",part.ig.label[i]);
}

if (option==1)
{
fprintf(f_trace,"\n i_group of particle %i   ",part.label);
for (i=0;i<part.ig.size;i++)
	fprintf(f_trace,"%i ",part.ig.label[i]);
}
}


 // ----------------------------------------------------------------------------- DISPLAY_MEMO
void   display_memo(int level)
{
  // Display memorized positions
  int d;
  int i;


  printf ("\n\n Memo level %i, size %i",level,memo[level].size);
  printf("\n Rank  Error      Position");
  for (i=0;i<memo[level].size;i++)
  {
    printf("\n %i.  ",i+1);

    printf("  %f",memo[level].error[i]);
    for (d=0;d<memo[level].x[i].p.size;d++)
    {
       printf("%10.3f ",memo[level].x[i].p.x[d]);
    }
  }

}

// ----------------------------------------------------------------------------- DISPLAY_problem
void display_problem(int level)
{
int 	d;
float 	p1,p2,p3;
int 	times;

printf("\n\n ***********************************************************************");
printf("\n\n PROBLEM");
printf("\n funct...........: %2i",problem[level].funct);
printf("\n Dimension.......: %2i",problem[level].DD);
printf("\n Search space:");
printf("\n xmin    xmax    dx (granul)");
p1=0;p2=0;
p3=-1;// Arbitrary unused value for dx

for (d=0;d<problem[level].H.dim;d++)
	{
	if (problem[level].H.min[d]!=p1 || problem[level].H.max[d]!=p2 || problem[level].H.dx[d]!=p3)
		{
		if (d>0) printf(" (%i times)",times);
		p1=problem[level].H.min[d];
		p2=problem[level].H.max[d];
		p3=problem[level].H.dx[d];
		printf("\n %f  %f  %f",p1 ,p2,p3);
		times=1;
		}
	else
		{
		times=times+1;
		}
	}
printf(" (%i times)",times);
if(problem[level].constrain>0) printf("\n constrain function %i",problem[level].constrain);
printf("\n target..........: %4.3f",problem[level].target);
printf("\n eps.............: %4.10f",problem[level].eps);
printf("\n printlevel......: %2i",problem[level].printlevel);
printf("\n times...........: %2i",problem[level].times);
if (confin_interv<0) printf("\nWarning, no intervall constraint");
if (problem[level].all_diff==1) printf("\nAll components of the solution must be different");
if (BIN==1) printf("\n This problem is seen as a special binary one");
if (no_best==1) printf("\n Best neighbour chosen at random");
if (no_best==2) printf("\n Best neighbour chosen using pseudo gradient");
if (problem[level].Max_Eval<0)
printf("\n Max eval. number: %i",-problem[level].Max_Eval);
else
printf("\n Max clock ticks.: %i",problem[level].Max_Eval);
if (problem[level].fuzzy>0) printf("\n I use fuzzification");
if (adapt==1) printf("\n Constant swarm size %i",problem[level].init_size);
printf("\n initial swarm size: %i",problem[level].init_size);
if (rand_hood>0) printf("\n I use random informer groups of size %i",rand_hood);
if (problem[level].save==-1) printf("\n Swarm size and  energies saved on energy.txt");
if(linkreorg==1) printf("\n I use link reorganization");
printf("\n Print level %i",problem[level].printlevel);


printf("\n---------\n");
}

// ----------------------------------------------------------------------------- DISPLAY_POSITION
void display_position(struct position pos)
{
int				d;
int				ncar=0;
struct position post;

post=pos;

printf("\n position size %i\n ",pos.p.size);
for (d=0;d<pos.p.size;d++)
{
	printf("%6.6f ",post.p.x[d]);
	ncar=ncar+10;
	if (ncar>75)
	{
		ncar=0;
		printf("\n");
	}
}

printf (" Error = ");
for (d=0;d<post.f.size;d++) printf(" %3.4f",post.f.f[d]);

}

// ----------------------------------------------------------------------------- DISPLAY_SWARM
void display_swarm(int type,int level)
{
	// TR is a global variable
int d;
int i;
int	j;
int	k;
int nb_f;

nb_f=problem[level].nb_f;
if (type<=1)printf("\n Current swarm");
if (type==2) printf("\n Best previous swarm");
if (type==3)
{
	printf("\n Particle f values");

goto end;
}

for (i=0;i<TR[level].size;i++)

{
	for (j=0;j<TR[level].tr[i].size;j++)
	{
		printf("\nLabel %i - f value(s): ",TR[level].tr[i].part[j].label);
		if (type==2)
         {
            for (d=0;d<nb_f;d++) printf(" %f / ",TR[level].tr[i].part[j].p_i.f.f[d]);   // best f value(s)
         }
		if (type==1)
         {
            for (d=0;d<nb_f;d++) printf(" %f  /",TR[level].tr[i].part[j].x.f.f[d]);   // current f value(s)
          }
		printf(". Position: ");
		for (k=0;k<TR[level].tr[i].part[j].x.p.size;k++)
		{
			if (type==2) printf("%f ",TR[level].tr[i].part[j].p_i.p.x[k]); // Best position
			if (type==1) printf("%f ",TR[level].tr[i].part[j].x.p.x[k]);     // Current position

		}
	}
}

end:;

}

// ----------------------------------------------------------------------------- DISPLAY_TRIBE
void display_tribe(FILE *f_trace,struct tribe tr, int option,int level)
{
int	i;

int	k;

if (option==0) // Display
{
	printf("\n label/status");
	for (i=0;i<tr.size;i++)
	{
		printf("%3i/%i ",tr.part[i].label,tr.part[i].status);
	}
}

if (option>=1)  // Save
{
	for (i=0;i<tr.size;i++)
	{
		fprintf(f_trace,"\n%3i/ ",tr.part[i].label);
		// Best position
				for (k=0;k<tr.part[i].p_i.p.size;k++)
               fprintf(f_trace,"%f ",tr.part[i].p_i.p.x[k]);
		// Links (labels)
				//fprintf(f_trace,"\n");
				display_i_group(f_trace,tr.part[i],1);

	}
}
}


// ----------------------------------------------------------------------------- PRINT_N_RUN
void print_N_run(int level,struct position best_result)
{
int               d;
double 		diff;
double		S;
double 		z,z1;
z=eval_f_tot[level]; // Evaluations here has not to be taken into account
z1=eval_f[level];

// Some stuff to print at the end

S=problem[level].target;
// Total error
diff=0;
  for(d=0;d<best_result.f.size;d++)
     diff=diff+ best_result.f.f[d]*best_result.f.f[d];
  diff=sqrt(diff);

printf("\n");

if(problem[level].printlevel>1)
{
	display_problem(level);
}

printf( "\nResult.................: ") ;  //show us how close it came
   for(d=0;d<best_result.f.size;d++) printf(" %4.10f",best_result.f.f[d]);
printf("  (wanted < %4.10f)",problem[level].eps);
if (fabs(S)>almostzero) printf( "\nRelative diff..........: %4.4f",fabs(diff/S));

if (problem[level].printlevel>0) display_position(best_result);

printf("\nTotal number of evaluations: %.0f",z);
printf("\nTotal number of clock ticks: %ld",clock_tick-1);

eval_f_tot[level]=z;
eval_f[level]=z1;

}


// ----------------------------------------------------------------------------- READ_PROBLEM
struct problem read_problem(FILE *f_problem, FILE *f_data,int level)
{
int		d,d_current,d_times;
int               m;
int               MM;
struct problem problemt;
float          t;
double z;

printf("\n-----------\n Read problem ");
strcpy(problemt.data,"None"); // Default. No special data file

problemt.init_file=0; // Default value: no data file for initial swarm. 1 else
strcpy(problemt.init_swarm_r,"None"); // Default value. Initial swarm NOT read from a file
problemt.used_constr=1; // Default value: use all constraints. <1 else (cf Frequency Assignment problem)

	//--------------------------------- Read problem file

   fscanf(f_problem,"%i",&problemt.nb_f);     // Number of objective functions (usually 1)
   if (problemt.nb_f<0)
   {
           problemt.nb_f=-problemt.nb_f;
           lexico=1;
   }
   else lexico=0;


fscanf (f_problem,"%i",&problemt.funct);

 if (problemt.funct==30) TSP=1;  else TSP=0;// Flag
  if (problemt.funct==32) QAP=1;  else QAP=0;// Flag


  if ( problemt.funct<0)
  {
     problemt.funct=-problemt.funct;
     problemt.fuzzy=1;
  }
  else
      problemt.fuzzy=0;
 /*
  PROGR=0;  // For progressive approach
  if (problemt.funct==55) PROGR=4;
  if (problemt.funct==56) PROGR=16;
  if (problemt.funct==57) PROGR=27;
  if (problemt.funct==58) PROGR=200;
  if (problemt.funct==59) PROGR=80;
  if (problemt.funct==60) PROGR=167;
 */


	fscanf (f_problem,"%i",&problemt.DD);

   	if (adapt==0)   // Complete adaptive method
    {
      // adapt=1;
       problemt.init_size=1; // Initial swarm size
    }
	else // Non adaptive, "adapt" is the swarm size
   {
      problemt.init_size=adapt; // Initial swarm size
      //adapt=0;
   }

	// Search space

	problemt.H.dim=problemt.DD;
	d_current=0;
	new_d:
	fscanf (f_problem,"%i",&d_times); // Number of times you define the same interval and granularity
	fscanf (f_problem,"%f",&problemt.H.min[d_current]);
	fscanf (f_problem,"%f",&problemt.H.max[d_current]);
	fscanf (f_problem,"%f",&problemt.H.dx[d_current]);
	d_current=d_current+1;
	if (d_times>1)
	{
		for (d=d_current;d<d_current+d_times-1;d++)
		{
			problemt.H.min[d]=problemt.H.min[d-1];
			problemt.H.max[d]=problemt.H.max[d-1];
			problemt.H.dx[d]=problemt.H.dx[d-1];
		}
		d_current=d_current+d_times-1;
	}
  mincoordinate=problemt.H.min[0]; // For hypercube dynamic problems (see movpeaks)
  maxcoordinate=problemt.H.max[0];

  if(problemt.H.dx[0]<0) BIN=1; else BIN=0; // Special binary problem

	if (d_current<problemt.DD) goto new_d;

      fscanf (f_problem,"%i",&problemt.constrain); // Constraint function (if > 0)

	fscanf (f_problem,"%f",&t); problemt.target=t;
	fscanf (f_problem,"%f",&problemt.eps);
	fscanf (f_problem,"%i",&problemt.all_diff);

	fscanf (f_problem,"%i",&problemt.printlevel);
	fscanf (f_problem,"%i",&problemt.save);
	fscanf (f_problem,"%i",&problemt.times);
	fscanf (f_problem,"%i",&problemt.Max_Eval);
    fscanf (f_problem,"%i",&problemt.Max_Eval_2);
    fscanf (f_problem,"%i",&problemt.Max_Eval_delta);
printf("\n %i %i",problemt.Max_Eval_2,problemt.Max_Eval_delta);
    if (TSP==1) // Read graph the very first time
    {
       problemt.P=read_graph_tsp(problemt.printlevel);
       // More particles for big  TSP problems
        //  if (problemt.DD>20 && adapt==0)
			  problemt.init_size=problemt.DD;
        //else
        //problemt.init_size=1;
     }

     if (QAP==1) // Read the matrices
     {
              problemt.P=read_graph_qap(problemt.printlevel);
       // More particles for big  QAP problems
         // if (problem[level].DD>20)
         //if (adapt==0)
              problemt.init_size=problemt.DD;
        //else
        //problem[level].init_size==1;
     }

	 if (problemt.constrain==1) // There are some discrete variables
								// whose values have to be read on a file
	 {
		fscanf(f_discrete,"%i", &discrete_nb);

		for (d=0;d<discrete_nb;d++)
		{
			fscanf(f_discrete,"%i", &discrete[d].d); // Variable rank
			fscanf(f_discrete,"%i", &discrete[d].size); // Number of values

			for (m=0;m<discrete[d].size;m++)
			{
				fscanf(f_discrete,"%f",&t); // Values
				discrete[d].v[m]=t;
			}
		}
	 }

     if (problemt.funct==29) // Read matrix for Cognitive dissonance
{
	// Read matrix P  (weights)
	printf("\n Read matrix");

	f_matrix=fopen("matrix.txt","r");
	fscanf(f_matrix,"%i",&MM);
      problemt.P.size=MM;
	for (d=0;d<MM;d++)
	{
		for (m=0;m<MM;m++)
		{
			fscanf(f_matrix,"%f",&t);
			problemt.P.val[d][m]=t;
		}
	}
	printf("\n Weigths read. %i lines %i columns",MM,MM);
}
   DYN=0; // Flag for dynamic optimisation
	 if (problemt.funct==61) // Moving peaks
   {
		  geno_size=problemt.DD;
		 init_peaks();
     DYN=1;
     }


   problemt.N=problemt.init_size;
   if(problemt.times>Max_run)
   {
           printf("\nWARNING. I can't run more than %i times this problem",Max_run);
           printf(" \n (you should increase Max_run in def_struct.c)");
           problemt.times=Max_run;
   }
   return problemt;
}


// ----------------------------------------------------------------------------- SAVE_SWARM
void save_swarm(FILE *f_swarm,int level )
{
int d;
int i;
int	j;
int	k;
//printf("\nSave swarm on swarm.txt");

// label    f_values      position   best_f_values   best_position

for (i=0;i<TR[level].size;i++)
{
	for (j=0;j<TR[level].tr[i].size;j++)
	{
         // label
         if (i==0 && j==0) fprintf(f_swarm,"\n* %i ",TR[level].tr[i].part[j].label);
		else fprintf(f_swarm,"\n- %i ",TR[level].tr[i].part[j].label);
         // f_values
		for (d=0;d<TR[level].tr[i].part[j].x.f.size;d++)
               fprintf(f_swarm," %f ",TR[level].tr[i].part[j].x.f.f[d]);
           // position

		for (k=0;k<TR[level].tr[i].part[j].x.p.size;k++)
               fprintf(f_swarm,"%f ",TR[level].tr[i].part[j].x.p.x[k]);

           // best_f_values
		for (d=0;d<TR[level].tr[i].part[j].p_i.f.size;d++)
               fprintf(f_swarm," %f ",TR[level].tr[i].part[j].p_i.f.f[d]);
           // best_position
		for (k=0;k<TR[level].tr[i].part[j].p_i.p.size;k++)
               fprintf(f_swarm,"%f ",TR[level].tr[i].part[j].p_i.p.x[k]);
	}
}
}

