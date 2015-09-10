/*************************************************************************
     TRIBES, an fully adaptive particle swarm optimiser
                             -------------------
    last update           : 2005-01-12
    email                    : Maurice.Clerc@WriteMe.com
    Home page           :  http://www.mauriceclerc.net
 ***************************************************************************/
 /*
     Have fun, and keep me posted if you do something interesting with it,
     or if you find a bug (I'm pretty sure there are some!).
 */
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/* Recent updates
2005-01
2004-11 Specific binary strategies (and examples of binary problems).
  Some of them with a very simple "taboo" rule.
2004-07
- added a pseudo-gradient method to define the best informer
- rewritten the "remove worst particle" part. There was a minor bug.
  The algo is so robust that most of the time it was not important.
  But sometimes it did ... So, globally, performances are now a bit better.

- added link_reorg. Not very convincing, though.

- you can now use discrete variables whose values can't be simply
   computed by using min, max and granularity. The acceptable values are read
   on the file discrete.txt
*/

/* TO DO

- adding something so that you can cope with more complicated search spaces
   (not only hyperparallelepids) WORK STILL IN PROGRESS.
   The perfect way, for a (semi) continuous search space would be by homeomorphism.
   Quite difficult, though.

 */

/*
************** TRIBES: An adaptive parameter free Particle Swarm Optimiser  ******************************
This version does not really need any parameters at all, except, of course, the ones
defining the problem (file problems.txt):
- function
- dimension
- xmin,xmax, granularity (for the hyperparallepipedic search space)
- target
- admissible error
- if all components of the solution have to be different (useful for some combinatorial problems)

Granularity. If you _know_  that the solution, on some dimensions,  is on an "integer" point, or even a rational one,
you can give this information
Fo example, if for a given dimension, the solution is an integer, then you can set granularity to 1 for this dimension:
the process may be faster (example: functions 4 and 6, or 1 and 2 if S is an integer)

Also, you have to set granularity to 1 for some "combinatorial" problems like Knap-sack, Fifty-fifty or Magic Square

You can try to optimize some simple (and less simple) functions, for example:
1) find DD numbers that sum to S.
2) find DD numbers whose product is S
3) Sphere. Find DD numbers  whose sum of squares is S
4) Rosenbrock (Banana)
5) Clerc's f1 (Alpine)
6) Griewank
7) Rastrigin
8) Fifty-fifty Problem (in particular with granul=1)
9) Ackley
10) Foxholes  2D.
  etc.
(see also the file functions.txt).
--------
 TSP
There is now a Traveling Salesman Problem option (30). Not that good, but it works for small
graphs (typically <20 nodes). The best way seems to be also the simplest:
just setting granularity to 1 and activating the all_different option.
However, some features (in particular initialization and local search) are specific (see TSP.c)

Of course, you need a file describing the graph (full matrix TSPLIB format).
I also tried a continous relaxation of the problem in another version, but it does not give good results.

Note that PSO for TSP may work quite well, but only by adding some local search that are not in Tribes.
 ---------
QAP
As TSP is a particular case of QAP (Quadratic Assignment Problem), I added this option (32), which works
almost exactly like TSP.
---------
Multiobjective   (see, for example, 33 etc.)
Using Tribes for multiobjective problems is very easy, for each run gives a different solution,
thanks to the random initialization of the first swarm (usually just one particle, at the very beginning).

The trick is to be able to compare two positions, no matter one or more functions are to take into account.
This is done in better_than().

Some points:
- set accuracy to 0.
- it never really "converges", so the stop criterion is the max number of evaluations you give, Max_Eval.
- I have added a module that "cleans up" the result file (run.txt), if you want keeping only the
  non dominated  solutions (in file run_cleaned.txt).

In practice, it is better to run a lot of  times (say 500)  with a small Max_Eval (say 20),
and then to clean up  the result.


Note 1
Sometimes, you may prefer using only positive functions.
   Then just add big enough constants. See 34 (Schaffer) for example.

Note 2
Multiobjective approach is a nice way to take complicated constraints into account.
See, for example Pressure vessel example (code function 52)

--------
You can add your own function (see myfunction.c and functions.txt)

Just for fun, using the function 2 (DD-product), you can then factorize a (small!) integer number.

Note:
For test purpose (cf. file problems.txt), there _are_ some parameters, in particular:
 no_begin:
  you may choose the "best" neighbour at random, or using or not  pseudo-gradient method

 adapt:
If adapt=n, it means you want to use a swarm of constant size n


 (and the neighbourhood of each particle is the whole swarm)
For classical test functions, and with n= about 20, this version is often better
than the ones that don't use hyperspheres.

If adapt=0, it means you want to use complete adaptive method.
Not always better than the previous one, but you don't even has to "guess"
what the right swarm size is.

Strategies:
There is no reason to adopt the same strategy for a "good" particle and for a "bad" one.
  These problemeters define which one to use in which case.
Note that, though, strategy 1 is a  compromise for any case.

non_unif coefficients:
Using these coefficients, distributions in hyperspheres can be non uniform
(Gaussian, for example).
In practice, uniform distribution seems OK.

Fuzzification:
You may use fuzzy decision rules and/or fuzzy distributions.
Doesn't seems to be very useful.

----------------------------------------------

Equation
--------
For each particle, we have
.
x(t+1) = x(t) + w(x(t),x(t-1),p_i(t),p_g(t))
with
.
x(t) := position at time t
x(t-1) := position at time t-1
pi(t):= previous best position of the particle (at time t)
pg(t):= previous best position found so far in the i-group (informer group) of the particle (at time t),
        including the particle itself
w(x1,x2,x3) := a vectorial function of three positions

Unlike "classical" PS0s, there is no explicit velocity
and w is not at all depending on some more or less arbitrary "weights"
 cognitive or social coefficients (often called alpha, phi_1, phi_2).
 A promising area is defined using hyperspheres (cf. move_particle).


Download
--------
Well, you have it, havn't you ?
For more information about PSO
 http://www.mauriceclerc.net, Math stuff for PSO
And for even more information, see the Particle Swarm Central http://www.particleswarm.net

Compiling the program
-------------------
It is written in ANSI C, so you shouldn't have any problem.

Using the program
-----------------
See files problems.txt
"Problems" describe the problem(s) (function, dimension, search space etc.

*/
/* INFORMAL DESCRIPTION

Common part
-----------
At the very beginning, there are n particles {0,1,2,...,n},

n=3 (or even 1) for the complete adaptive version. So we have on the same time
- the swarm,
- a tribe T0,
- three i-groups I0, I1,I2 (informers group).
Each particle belongs to the tribe T0
For each particle, its i-group contains:
	i) the particle itself
	ii) the best particle of the tribe


At each time step each particle "moves" towards a more promising area,
using its knowledge (x,p_i,p_g).
- let H_i be the hypersphere radius=norm(p_i-p_g), center p_i
- let H_g be the hypersphere radius=norm(p_i-p_g), center p_g
- a point pp_i is randomly chosen in H_i
		(WARNING. It is not that easy, see rand_in_hypersphere)
- a point pp_g is randomly chosen in H_g
- x(t+1) is a weighted combination of pp_i and pp_g,

	according to the f values on p_i and p_g.

Note: if you don't use adaptive version, it means you have a global PSO
(the neighbourhood of each particle is the whole swarm).

Adaptive part
-------------
From time to time social adaptation is performed:
- for each tribe, one counts how many particles have improved their best performance (n)
	If n=0,the tribe is said "bad"
	if (n>=size(tribe)/2) the tribe is said "good"

- if there is at least one bad tribe, a new tribe is generated
- for each bad tribe:
	- the best particle is found and "generates" a new particle
		in the new tribe, for the moment purely at random,
	- this new particle is added to the i-group of this best particle

- for each particle of the new tribe, the i-group contains
	ii) the particle itself
	ii) the best particle of the new tribe,
	iii) the particle that has generated it (symmetry. Note that thanks to the i-group
			concept, you could build non symmetrical relationships)



- for each "good" tribe
	- the worst particle is removed (if it is bad)
	- in all i-groups it is replaced by the best one of the tribe (general case)
		or by the best particle of the i-group of the removed particle
		(mono-particle tribe case)
	- its i-group (except of course this removed particle itself) is merged to the one
		that replaces it.


Example
-------
Very beginning
T0 = {0,1,2}
i-group table
particle	i-group
	0		{0,1,2}
	1		{0,1,2}
	2		{0,1,2}
========
Adaptation 1, T0 is "bad", best particle 1.
A new tribe T1 is generated.
3 is generated by 1 in the new tribe T1

T0 = {0,1,2}
T1 = {3}

 i-group table
 --------------
	particle	i-group
		0		{0,1,2}
		1		{0,1,2,3}
		3		{1,3}
=========
Adaptation 2. T0 is good, best particle 0, worst 1, T1 is bad, best particle 3.
T2 is generated.
1 is removed (and replaced by 0 for the i-groups)
4 is generated by 0 into T2.


T0 = {0,2}
T1 = {3}
T2 = {4}
 i-group table
 --------------
	particle	i-group
		0		{0,2,3}

		2		{0,2}
		3		{0,3,4} (1 has been replaced by 0)
		4		{3,4}
==========
Adaptation 3. T0 is bad (best particle 2, worst particle 0),T1 is good, T2 is bad.
T3 is generated

0 is removed (replaced by 2 for i-groups)
6 is generated by 2 in T3
8 is generated by 4 in T3

T0 = {2,5}
T1 = {3}
T2 = {4}
T3 = {6,8}

 i-group table
 --------------
	particle	i-group
		2		{2,5,6}
		3		{3,7}
		4		{4,8}
		5		{2,5}
		6		{2,6,8}
		7		{3,7}
		8		{4,8}

etc.

 */

// DEPENDING on how you build the application you may or may not need to "include"
// the following files

 #include "def_struct.c"

 // External subroutines
 #include "movpeaks_mc.c"
 #include "read_display_save.c"
 #include "tools.c" // Some  mathematical tools like alea(), max(), min(), regranul() etc.
 #include "extra_tools.c" // Some "informative" tools, not strictly necessary, like distance(), energy() etc.
 #include "myfunction.c" // Write your own objective function in this file.
 #include "myconstrain.c"     // Write your own constrain function g(position)>=0
 #include "apple_trees.c" // For "apple trees" example

 #include "TSP.c"
 #include "QAP.c"
 #include "MINLP.c"

 #include "annxor.c"
 #include "annparity.c"
 #include "anncolorcube.c"
 #include "annpima.c"
 #include "annsinsimp.c"
 #include "annservo.c"

 //=========================================================================================
int main(int argc, char *argv[])
{
int           bidon;
int				i,j;
int            level;
float   Max_Eval;
int				nb_pb;
struct position result;
double t;

E=exp(1);
pi=(double)2*acos(0);
two_pi=2*pi;

 printf("\n t  %f",20.0); // Test of code generation. Should print 20.000000

// Initialize function names. Just for display
  if ((f_functs = fopen("functions.txt","r")) == NULL)
   {
   fprintf(stderr,"Can't open the file functions.txt\n");
   exit(1);
   }

i=1;
next_funct:
	fscanf(f_functs,"%i %s\n",&bidon,&functions[i]);
  if (bidon!=-1)
  {
	// printf("\n %i %s",bidon,functions[i]);
    i=i+1;
    goto next_funct;
  }
printf("\n-----------------------------------\n");

fclose (f_functs);


//===================  Other files  (keep them open during the whole process)


f_problem=fopen("problems.txt","r");
if (f_problem == NULL)
 {
 fprintf(stderr,"Can't open the file problems.txt\n");
 exit(1);
 }
 f_discrete=fopen("discrete.txt","r");
 if (f_discrete == NULL)
 {
 fprintf(stderr,"Can't open the file discrete.txt\n");
 exit(1);
 }
printf("\n Open files");
f_run=fopen("run.txt","w");
f_run_c=fopen("run_c.txt","w");
f_swarm=fopen("swarm.txt","w");
f_synth=fopen("synth.txt","w");
f_trace=fopen("trace.txt","w");
f_energy=fopen("energy.txt","w");
f_trace_run=fopen("trace_run.txt","w");
//-----------------

// Read strategies  (for tests. Will be hard coded later) . See move_particle()
for (i=0;i<2;i++)
{
      for (j=0;j<Max_status;j++)
         fscanf(f_problem, "%f",&strategies[i][j]);

}

fscanf(f_problem,"%i",&confin_interv); // If <0, no interval (boundary) confinement

fscanf(f_problem,"%i",&circular_hood); // Option. Usually equal to 0.
// If >0,  indicates the size of the circular neighbourhood
//    This option is just to rapidly compare with parametric PSO

if (circular_hood<0)
{
	if (circular_hood<-1)
	{
		circular_hood=-circular_hood;
		rand_hood=circular_hood;
	}
	else
	{
		//rand_hood=K_optim(problem[level].init_size,problem[level].eps);
		// Just for info. : init_size and eps are not yet known
	}
}
else
{
	rand_hood=0;
}

// Check whether positions have to be memorized or not
MEMO=0;
for (j=0;j<Max_status;j++)
{
   if (strategies[0][j]==18 || strategies[0][j] ==19) MEMO=1;
}

fscanf(f_problem, "%i",&no_best);
fscanf(f_problem,"%i",&adapt); // Usually 0. If >0, it is a constant swarm size
fscanf(f_problem,"%i",&linkreorg);

fscanf(f_problem,"%i",&nb_pb_max); // Number of different problems

nb_pb=0;
problem:
printf("\n Problem number %i",nb_pb);
// Set level to the main one (normal case)
 level=0;

  // Description of the problem(s)
  problem[level]=read_problem(f_problem,f_data,level); // READ THE PROBLEM DESCRIPTION

   if (TSP==1 || QAP==1) // For TSP or QAP problem, the graph has been read in read_problem
                                       // Check it if is consistent with the data
   {
          if(problem[level].DD!=problem[level].P.size)
          {
            printf("\n ERROR. Graph size %i not consistent with problem data  %i",problem[level].P.size, problem[level].DD);
           return EXIT_SUCCESS;
          }
   }

  display_problem(level);
   times=0; // you can run several times the same problem (problemeter problem[level].times)

   // Print titles on the f_run file
  fprintf(f_run,"Function Run Target"  );
  for (i=0;i<problem[level].nb_f;i++) fprintf (f_run," f%i",i+1);
  fprintf(f_run," Eval_nb Duration Added Removed"  );
   for (i=0;i<problem[level].DD;i++) fprintf( f_run," x%i",i+1);

  if (problem[level].funct==11) // Apples trees (more generally, use of homogeneous coordinates
  {
    fprintf(f_run," Real position");
    for (i=1;i<problem[level].DD-1;i++) fprintf( f_run," _");
   }

  // Print titles on the f_synth file
  fprintf(f_synth,"Function Dim. Target Wanted_error Nb_of_runs Mean_eval._nb Std_dev") ;
  fprintf(f_synth," Success_rate");
  fprintf(f_synth," Min_eval Max_eval Mean_swarm Mean_error Std_dev Mean_error_cont Std_dev Min_error Max_error Duration") ;

  problem[1]=problem[0];
  problem[1].DD=1;


  // Compute coefficients for some mouve options
  	coeff_S_C=coeff_SC(problem[level].DD);
	phi=2/0.97725;
	khi=1/(phi-1+sqrt(phi*phi-2*phi));
	cmax=khi*phi;

  // Display function and dimension
  printf("\n%s %iD",functions[problem[level].funct],problem[level].DD);

  seed_rand_kiss(1); // Initialize pseudo random number generator

nprogr=1; // For progressive problems like Neural Network Training

//core:
printf("\n %i %i",problem[level].Max_Eval, abs(problem[level].Max_Eval_2));
for (Max_Eval=problem[level].Max_Eval;fabs(Max_Eval)<=fabs(problem[level].Max_Eval_2);Max_Eval=Max_Eval+problem[level].Max_Eval_delta)
{
   // Arbitrary very high values for best of the best, to begin
  for (i=0;i<problem[level].nb_f;i++)  BEST.f.f[i]=infinite;
  printf("\n**********\n Run PSO with Max_Eval= %.0f",Max_Eval);
  times=0;
  result=PSO(level, Max_Eval); // ******* CORE OF THE PROGRAM *******
}
//nprogr=nprogr+1; if (nprogr<=PROGR) goto core;
nb_pb=nb_pb+1; if(nb_pb<nb_pb_max) goto problem;

  return EXIT_SUCCESS;
}

/*============================ S U B R O U T I N E S =====================================*/
/* Except the first one, they are in alphabetical order                                   */

struct position PSO(int level, float Max_Eval)  // Main routine.
{
int				add_option=1; // 0 => generate just random particles
							  //  1 => generate also a hopefully good one
int            already;
int             bad;
//struct position best_result;

int         better;
int            chance;

int				cycle,cycle_size;
int				d;
//int				dim;
double         duration;
double         duration_tot=0;
double         Ek; // Kinetic energy
double			Ep; // Potential energy
double         eps_run[Max_run][2];
double         error,error_cont;
double         eval_run[Max_run];
int				failure_tot;
struct f         f0,f1;
int             good;
int				i,i0;
struct i_group    ig,igt;
int 			Iter;
double            improv;
int                  k;
int				j;
int               l;
int         label_best,label_worst;
int         labelt;
int         m;
int               max_add;
double			max_eps;

double               max_f;
double			mean_eps,mean_eps_cont;
double          mean_eval;
double			mean_swarm_size;
double			min_eps;
int				n_add;
int             n_connect;
int             n_local_remove;
int				n_remove;

double			nb_eval_max;
double			nb_eval_min;
double			nb_no_progress;
int                     new_label;
struct	particle	part_best;
struct particle	partg;
struct	particle	partt;
struct position   post;
double               pr;
struct f     previous_best;

int                     rank;
int                     rank_best,rank_worst;
double               sigma_eps,sigma_eps_cont;
double               sigma_eval;

int                     size;
double      strateg[nb_strateg]; // Just for information (see move_particle)
int			swarm_size;
clock_t 		ticks0,ticks1, ticks2;
int                     TRsize0;
//double			volume; // Just for info
//double			x_max,x_min;
//double      x1,x2;
//int				worst[Max_swarmsize]; // Label of worst particle in each tribe
double			z,zzz;

 zzz=0;
nb_eval_max=0;
nb_eval_min=infinite;
failure_tot=0; // Number of failures
min_eps=infinite;   // This is a _minimization_ problem
mean_eps=0;
max_eps=0;

for (i=0;i<nb_strateg;i++) strateg[i]=0; // Just for information about the strategies used
for (i=0;i<Max_status;i++) status_count[i]=0; // Just for information about particle status

//======================================== You may solve several times the same problem

times_loop:
 memo[level].size=0;
 /* Useful only when using pivot method
     Note: for implementation of a future ReHope method (Restart)
     it would be interesting to put this instruction BEFORE the restart loop
     in order to keep memory of best results. In such a case, you also have to put the
      other "... =0" before the loop

 */
eval_f_tot[level] =0; // Number of objective function evaluations
eval_f[level]=0; // Number of evaluations counter
offline_error=0; // Useful only for dynamic optimisation
offline_error_cont=0;

 ticks1=clock(); // Just for information. To evaluate processor time needed

chance=0; // Equal to 1 if a solution is found

nb_no_progress=0;// If this number becomes > parapet, then => Failure
Iter=1;		// Number of iterations. Just for information
n_add=0;
n_change=1;

n_remove=0;
mean_swarm_size=0;
label[level]=0; // Label of the first particle. Just for later visualisation
 ticks0=clock();

Xmin.size=problem[level].DD;
Xmax.size=Xmin.size;
 for(d=0;d<Xmin.size;d++) // Initialize Xmin and Xmax
 {
	 Xmin.x[d]=problem[level].H.min[d];
	 Xmax.x[d]=problem[level].H.max[d];
 }

//srand((unsigned)time(NULL)); // Re initialize random number generation, according to the time
//srand(1);  // Re initialize random number generation with the same seed

	// Check improvement every cycle_size time steps

	cycle_size=problem[level].init_size; // Initial value


							// ******** SWARM INITIALIZATION *********
if (problem[level].printlevel>=0)
   printf("\n------------------------------------------------------------------------------------------------------------");
if (problem[level].printlevel>0)
	printf(" \n Run %i/%i",times+1,problem[level].times);

// Initialize the first tribe
best_result.p.size= problem[level].DD;
best_result.f.size= problem[level].nb_f;
previous_best.size=problem[level].nb_f;
TR[level].size=1; // Number of tribes
TR[level].tr[0].size=problem[level].N; // Tribe size
 max_f=0;

 // Generate particles and keep the best result found so far
 //printf("\n Run %i",times+1);
 if (TSP==1)  // Special initialization for TSP
 {
    TR[level].tr[0]=init_swarm_tsp(problem[level].N,problem[level].target,problem[level]. printlevel,level);
    goto end_init;
 }

  if (QAP==1)  // Special initialization for QAP
 {
    TR[level].tr[0]=init_swarm_qap(problem[level].N,problem[level].target,problem[level].printlevel,level);
    goto end_init;
 }

 // Normal initialization (random)
   for (i=0;i<TR[level].tr[0].size;i++)
   {
	TR[level].tr[0].part[i]=init_particle(0,level,dummy_part);
   }

 end_init:    // End of initialisation

//   display_swarm(1,level);

best_result=TR[level].tr[0].part[0].x;
previous_best=best_result.f;

for (i=1;i<TR[level].tr[0].size;i++)
{
//display_position( TR[level].tr[0].part[i].x);
   // Look for the best result
    better=better_than(TR[level].tr[0].part[i].x.f, best_result.f,level);

 	if (better==1)
	{
      previous_best= best_result.f;
    best_result=TR[level].tr[0].part[i].x;
	}
	offline_error_cont+=total_error(best_result.f);

}  // next i

if (BIN==1)
	for (d=0;d<best_result.p.size;d++)
		best_result.p.x[d]=floor(best_result.p.x[d]+0.5);

    // Memorize the best value found so far
		error=total_error(best_result.f);
		offline_error=error;
		offline_error_cont=error;


    if (error<min_eps) min_eps=error;

// In case of progressive approach, keep the best of the best
if (nprogr>1)    TR[level].tr[0].part[0].x=BEST;

 // Display best result
//if (problem[level].printlevel>0)
 {
    printf("\n Best result after initialization");
    display_position(best_result);
    printf("\n");
 }

 // Check if solution found by chance
if (DYN==0)  if (error<=problem[level].eps) goto end_by_chance;

 // Memorize positions   ("Black board" optional method)
 if (MEMO==1)
 {
    for (i=0;i<TR[level].tr[0].size;i++)
    {
       add_memo(TR[level].tr[0].part[i].x,level);
    }
 }


// Initialize i-groups
if (circular_hood>0) goto circular; // Option to simulate classical PSO
if (rand_hood>0) goto rand_i_group;

	// First particle
TR[level].tr[0].part[0].ig.size=TR[level].tr[0].size; // i-group = all particles
for (j=0;j<TR[level].tr[0].part[0].ig.size;j++)
{
	TR[level].tr[0].part[0].ig.label[j]=TR[level].tr[0].part[j].label;
}

	// Same i-group for all the others
for (i=1;i<TR[level].tr[0].size;i++)
   TR[level].tr[0].part[i].ig=TR[level].tr[0].part[0].ig;
 goto before_cycles;

 //------------------------------
 circular: // W A R N I N G
 // Valid only for constant swarm size option (and, then, just one tribe, rank 0)
  size= TR[level].tr[0].size;

 for (i=0;i<size;i++)
 {
    TR[level].tr[0].part[i].ig.size=circular_hood;

    for (j=0;j<TR[level].tr[0].part[i].ig.size;j++)
    {
        rank=i+j-(int)(circular_hood/2.0);
        if (rank<0) rank=rank+size;
        else
           if (rank>=TR[level].tr[0].size)  rank=rank-size;

       TR[level].tr[0].part[i].ig.label[j]=TR[level].tr[0].part[rank].label;
    }
}

 goto before_cycles;
//--------------------------  Option random i-groups
rand_i_group:
swarm_size=problem[level].init_size;
for(k=0;k<TR[level].size;k++)
{
     size= TR[level].tr[k].size;
     for (i=0;i<size;i++)
     {
        TR[level].tr[k].part[i].ig.size=rand_hood;

        for (j=0;j<TR[level].tr[k].part[i].ig.size;j++)
        {
           rank=alea(0,swarm_size);
           TR[level].tr[k].part[i].ig.label[j]=TR[level].tr[k].part[rank].label;
        }
     }
}


//------------------------
before_cycles:
   if (problem[level].save==-1)  // Save swarm size and energies
										// Just for information
        {
           Ek=energy_k(level); Ep=energy_p(level);
		   swarm_size=problem[level].init_size+n_add-n_remove;
           fprintf(f_energy,"\n%i %i %i %i %f %f",(int)eval_f_tot[level],swarm_size,n_add,n_remove,Ek,Ep);

        }
	if (problem[level].save==-2) // Save info for convergence curve plotting
	{
		fprintf(f_trace_run,"Eval error error_cont n_add");
		fprintf(f_trace_run,"\n%.0f %f %f 0",eval_f_tot[level],error,error);
	}

recent_change=0;  // For Moving Peaks

// cycles:
cycle=0;  // Number of iterations since last adaptive update
                                                                                                                // ******** SWARM MOVES *********

do_Iter: //termination criterion is (currenteps<=eps or eval_f_tot>stop)

	if (problem[level].printlevel>0)
	{
		printf("\n  Best result after %.0f evaluations:  ",eval_f_tot[level]);
         for (i=0;i<best_result.f.size;i++) printf (" %f",best_result.f.f[i]);
		 if (DYN==1) printf("\n offline error %f",offline_error/n_change);

       }
     if (problem[level].printlevel>1)
     {
		display_position (best_result);
		printf("\n------ MOVING the %i tribe(s)",TR[level].size);
	}

	if (problem[level].printlevel>1)

	{

		printf("\n %i tribe(s)",TR[level].size);
	}
	if (problem[level].save>=2) fprintf(f_trace,"\n %i tribe(s)",TR[level].size);


	for (i=0;i<TR[level].size;i++) // For each tribe ...
	{
		if (problem[level].printlevel>1)
		{
			display_tribe(f_trace,TR[level].tr[i],0,level);
		}

		if (problem[level].printlevel>2)
		{
			for (j=0;j<TR[level].tr[i].size;j++)
				display_i_group(f_trace,TR[level].tr[i].part[j],0);
		}

		if (problem[level].save>0)
		{
			fprintf(f_trace,"\n  tribe %i",i);

			display_tribe(f_trace,TR[level].tr[i],1,level);
		}

		for (j=0;j<TR[level].tr[i].size;j++) // ... for each particle ...
		{
         // Find the "best informer" (the true one, or at random)
			partg=best_informer(TR[level].tr[i].part[j],no_best,level);

         // ******************* MOVE ************************;

			TR[level].tr[i].part[j]=move_particle(TR[level].tr[i].part[j],partg,level,0);
			if(MEMO==1)  add_memo(TR[level].tr[i].part[j].x,level);

			// Eventually update the best result
               better=better_than(TR[level].tr[i].part[j].p_i.f, best_result.f,level);


			   offline_error_cont+=total_error(best_result.f);

			if (better==1)
			{
			  // For dynamic optimisation, it is better
			  // to permanently evaluate the error

				offline_error+=-total_error(best_result.f);

				previous_best=best_result.f;

				best_result=TR[level].tr[i].part[j].p_i;

				if (BIN==1)
	for (d=0;d<best_result.p.size;d++)
		best_result.p.x[d]=floor(best_result.p.x[d]+0.5);

				error=total_error(best_result.f);
				offline_error+=error;
//   printf(" %f",  offline_error);

				if (DYN==1)
					error=offline_error/n_change; // For non recursive dynamic optim

                  	if (problem[level].printlevel>0)
                      {
                         printf("\n");
                         printf(" Eval: %.0f",eval_f_tot[level]);
                         printf(" Result: %.8f",error);
                      }

					if (problem[level].save==-2) // Save info for convergence curve plotting
					{
						fprintf(f_trace_run,"\n%.0f %f %f %i",eval_f_tot[level],error,offline_error_cont/(eval_f_tot[level]-n_add),n_add);
					}

                      if (problem[level].save>0)
                      {
                         fprintf(f_trace,"\n %f",best_result.f.f[0]);
                      }
                     if(DYN==0 && error<=problem[level].eps) chance=1;  else chance=0;
                      if(chance==1) goto end_by_chance; // A solution has been found

			}
   }  // next j

			// Arbitrary stop criterion
			ticks2=clock();
			clock_tick=ticks2-ticks1;

			if((Max_Eval<0 && eval_f_tot[level]>-Max_Eval) || Max_Eval>clock_tick)
			{
				failure_tot=failure_tot+1; // Just for information
      	       // if (problem[level].printlevel>0)
				printf("\n FAILURE eval_tot  %.0f > max= %.0f",eval_f_tot[level],fabs(Max_Eval));
				goto the_end;
			}

		} // next i


 	if (problem[level].printlevel>1) display_swarm(1,level);

   // Special local search for TSP
   if (TSP==1)
   {
      for (i=0;i<TR[level].size;i++) // For each tribe ...
	    {
         for (j=0;j<TR[level].tr[i].size;j++) // ... for each particle ...
		    {
			  TR[level].tr[i].part[j].x=local_search_tsp(TR[level].tr[i].part[j].x,problem[level].target,level);

			   // If improvement, memorize
               better=better_than(TR[level].tr[i].part[j].x.f, TR[level].tr[i].part[j].p_i.f,level);
               if (better==1)
               {
                 TR[level].tr[i].part[j].p_i= TR[level].tr[i].part[j].x;
               }
            }
      }
   } // end if(TSP==1)


      // Special local search for QAP
   if (QAP==1)
   {
      for (i=0;i<TR[level].size;i++) // For each tribe ...
		{
			for (j=0;j<TR[level].tr[i].size;j++) // ... for each particle ...
			{
               TR[level].tr[i].part[j].x=local_search_qap(TR[level].tr[i].part[j].x,problem[level].target,level);
               // If improvement, memorize
               better=better_than(TR[level].tr[i].part[j].x.f, TR[level].tr[i].part[j].p_i.f,level);
               if (better==1)
               {
                 TR[level].tr[i].part[j].p_i= TR[level].tr[i].part[j].x;
               }
            }

      }
   } // end if(QAP==1)


	// Total particle number as criterion
	cycle_size=problem[level].init_size+n_add-n_remove;
	swarm_size=cycle_size;

   if (cycle_size>Max_swarmsize)
   {
     printf("\nWARNING. You should increase Max_swarmsize (%i) in def_struct.c",Max_swarmsize);
    cycle=Max_swarmsize;

   }
	mean_swarm_size=mean_swarm_size+cycle_size;

  	Iter=Iter+1; // number of time steps. Just for information

	if (adapt>0) goto end_adapt; // Non adaptive option (constant swarm size)

    //                                ***********  ADAPTATION(S)  (begin) *************
  // Number of connections
  n_connect=0;
 for (i=0;i<TR[level].size;i++) // For each tribe
{
	for (j=0;j<TR[level].tr[i].size;j++)  // for each particle
	{
      n_connect=n_connect+TR[level].tr[i].part[j].ig.size;
	}
}

	cycle=cycle+1;
if (cycle<n_connect)
        goto end_adapt; // Adaptation is not made at each iteration
										 // just from time to time

		TRsize0=TR[level].size; // Memorize the current number of tribes

		// Look for Xmin and Xmax. Useful only if you use init 3 (see init_particle)
		if (swarm_size>1) // Needs several particles
		{
			for (d=0;d<problem[level].DD;d++) // loop on dimensions
			{
				Xmin.x[d]=problem[level].H.max[d];
				Xmax.x[d]=problem[level].H.min[d];

				for (i=0;i<TRsize0;i++) // Loop on tribes
				{
					for (j=0;j<TR[level].tr[i].size;j++) // loop on particles
					{
						z=TR[level].tr[i].part[j].p_i.p.x[d];
						if(z<Xmin.x[d]) Xmin.x[d]=z;
						if(z>Xmax.x[d]) Xmax.x[d]=z;

					}
				}
			}
		}

		// If improvement, "bad" particles are not useful => remove them
		//  If deterioration, it is worthwhile to generate a very different particle


		for (i=0;i<TRsize0;i++) // Loop on tribes, to check them
		{
			// Find the best particle in the current tribe
			label_best=TR[level].tr[i].part[0].label; // Label of the best
			f0=TR[level].tr[i].part[0].p_i.f;
			rank_best=0;

			if (TR[level].tr[i].size>1)
			{
				for (j=1;j<TR[level].tr[i].size;j++)
				{
					f1=TR[level].tr[i].part[j].p_i.f;
					if (better_than(f1,f0,level)==1)
					{
						label_best=TR[level].tr[i].part[j].label;
						f0=f1;
						rank_best=j;
					}
				}
			}

         		// Find the worst particle in the current tribe
           // (maybe for future version)
				label_worst=TR[level].tr[i].part[0].label; // Label of the worst
				f0=TR[level].tr[i].part[0].p_i.f;
				rank_worst=0;

				for (j=1;j<TR[level].tr[i].size;j++)
				{
                        f1=TR[level].tr[i].part[j].p_i.f;
                        if (better_than(f1,f0,level)!=1)
                        {
                           label_worst=TR[level].tr[i].part[j].label;
                           f0=f1;
                           rank_worst=j;
                        }
				}
				//worst[i]=TR[level].tr[i].part[rank_worst].label;

			// Count how many particles of the current tribe have improved their position
			improv=local_improv(TR[level].tr[i],level);

      if (problem[level].fuzzy==0)   // Crisp decision rules. Note a tribe may be both bad an good
      {
           bad=improv<TR[level].tr[i].size; // At least one particle has NOT improved its position
		   good=improv>0; // At least one particle has improved  its position
      }
      else   // Use fuzzy rules to decide whether the tribe is bad or good
      {
         z= TR[level].tr[i].size;
		 pr=alea_float(0,z);

		 bad=(improv<=pr);
		 good=1-bad;
      }


if (bad==0) bad=alea(0,1); // To be sure there are on the whole
							// more generations than deletions
	  if (bad==1)
        {
					// If it is the first bad tribe, generate a new empty one
					if (TR[level].size==TRsize0)
					{

						TR[level].size=TR[level].size+1;
						TR[level].tr[TRsize0].size=0;
					}

              // max_add= TR[level].tr[i].size+1;
             // max_add=alea(1,TR[level].tr[i].size); // TEST. May generate several particles
				max_add=1;
				if (add_option==1) max_add=2;

              if (good && (TR[level].tr[i].size>1 || TR[level].tr[i].part[0].status>0)) //The worst particle will be removed
                 max_add=max_add+1;

               for (m=0;m<max_add;m++)  // May generate several particles
               {

				   //Generate a particle in the new tribe (rank of the tribe=TRsize0)
                  if (TSP==1)
                  {
					TR[level].tr[TRsize0].part[TR[level].tr[TRsize0].size]=init_particle_tsp(problem[level].DD,problem[level].target,level);
					goto new_label;
                  }
                  if (QAP==1)
                  {
					TR[level].tr[TRsize0].part[TR[level].tr[TRsize0].size]=init_particle_qap(problem[level].DD,problem[level].target,level);
					goto new_label;
                  }

				  if (add_option==0 || (add_option==1 && m!=1)) // Random init
					TR[level].tr[TRsize0].part[TR[level].tr[TRsize0].size]=init_particle(0,level,dummy_part);

					if (add_option==1 && m==1) // Generate a hopefully good particle
					{
						// Find the best informer of the tribe
						part_best=TR[level].tr[i].part[0];
						for (j=0;j<TR[level].tr[i].size;j++) // ... for each particle ...
						{
						partg=best_informer(TR[level].tr[i].part[j],0,level);
					//	partg=best_informer(TR[level].tr[i].part[j],no_best,level);

						better=better_than(partg.p_i.f,part_best.p_i.f,level);
						if (better==1) part_best=partg;
						}

						// Generate a particle around the best position known by the best informer
						partg=TR[level].tr[i].part[label_best]; // Best particle of the tribe
						partt=move_particle(partg,part_best,level,1);
						partt=complete_part(partt,0,level);
						partt.ig.size=0;
						TR[level].tr[TRsize0].part[TR[level].tr[TRsize0].size]=partt;
					}
// ---

new_label:
					new_label=TR[level].tr[TRsize0].part[TR[level].tr[TRsize0].size].label;
					if(rand_hood>0) goto end_i_groups;

					TR[level].tr[TRsize0].size=TR[level].tr[TRsize0].size+1; // Increase the (new) tribe size
					n_add=n_add+1;

					// Begin to build its i-group:
					// with the label of the best particle that has generated the new particle

					// Add it to the i-group of the new particle
					rank=TR[level].tr[TRsize0].part[TR[level].tr[TRsize0].size].ig.size;
					TR[level].tr[TRsize0].part[TR[level].tr[TRsize0].size].ig.label[rank]=label_best;
					TR[level].tr[TRsize0].part[TR[level].tr[TRsize0].size].ig.size=rank+1; // Increase the i-group size

					// Update the i_group of the "generating" best particle ...

					rank=TR[level].tr[i].part[rank_best].ig.size;
					TR[level].tr[i].part[rank_best].ig.label[rank]=new_label;
					TR[level].tr[i].part[rank_best].ig.size=rank+1;
/*
					// ... or ... Update the i_groups of all particles of the "generating" tribe
					// (if you use this option, do not use the previous one)
					for (j=0;j<TR[level].tr[i].size;j++)
					{
						rank=TR[level].tr[i].part[j].ig.size;
						TR[level].tr[i].part[j].ig.label[rank]=new_label;
						TR[level].tr[i].part[j].ig.size=rank+1;
					}
*/

                  end_i_groups:;
               } // m<max_add
			}   // End "if (bad)"

			// If  enough improvement
		if (good==1)
			{
				n_local_remove=0;

				// If there is just one particle, check whether an informer is better
				if (TR[level].tr[i].size==1)
				{
					// Find the best informer
					partg=best_informer(TR[level].tr[i].part[0],no_best,level);

					if (better_than(partg.p_i.f,TR[level].tr[i].part[0].p_i.f,level)==1)
          { label_best=partg.label;
            goto remove;
          }

					goto exit_good;
				}
remove:
					// Remove the worst particle of the current tribe

					n_local_remove=n_local_remove+1;
					// Temporarily keep its i-group
					ig=TR[level].tr[i].part[rank_worst].ig;

					// Remove it from the tribe (compact the tribe)
					if (rank_worst<TR[level].tr[i].size-1)
					{
						for (j=rank_worst;j<TR[level].tr[i].size-1;j++)
						{
							TR[level].tr[i].part[j]=TR[level].tr[i].part[j+1];
						}
					}
					TR[level].tr[i].size=TR[level].tr[i].size-1;
					n_remove=n_remove+1;

					// Replace it by the best in all i-groups
					for (j=0;j<TR[level].size;j++) // For each i_group...
					{
						for (k=0;k<TR[level].tr[j].size;k++) // ... for each particle ...
						{
							// Check if the best is already in the i-group

							already=-1;
							for (l=0;l<TR[level].tr[j].part[k].ig.size;l++)
								if (TR[level].tr[j].part[k].ig.label[l]==label_best)
									{already=l;}

//check_inf:

							igt.size=0;
							for (l=0;l<TR[level].tr[j].part[k].ig.size;l++) //... for each informer...
							{
								labelt=TR[level].tr[j].part[k].ig.label[l];
								if (labelt!=label_worst)
								{igt.label[igt.size]=labelt; igt.size=igt.size+1;}
							}
             if(already==-1 && label_worst!=label_best)
               {igt.label[igt.size]=label_best; igt.size=igt.size+1;}
							TR[level].tr[j].part[k].ig=igt;

//next_part:;
						}
					}  // End "for each i-group"
                         if(rand_hood>0) goto end_i_groups_2;
					// Merge its i-group (except itself) into the one of the best
					for (j=0;j<ig.size;j++) // Loop on informers
					{
						if (ig.label[j]==label_worst) continue;
						size=TR[level].tr[i].part[rank_best].ig.size;
						for (k=0;k<size;k++)
						{
							if (ig.label[j]==TR[level].tr[i].part[rank_best].ig.label[k]) goto next_informer;
						}
						TR[level].tr[i].part[rank_best].ig.label[size]=ig.label[j];
						TR[level].tr[i].part[rank_best].ig.size=size+1;

					next_informer:;
					}
               end_i_groups_2:;
				} // End "if (good)"

exit_good:;
		} // End "for each tribe"

 if (rand_hood>0) goto empty_tribe;

       		// Complete the i-groups of each particle of the new tribe (if there is one)
		// with all particles of the tribe
      if (TR[level].size>TRsize0)
		{
			for (i=0;i<TR[level].tr[TRsize0].size;i++) // For each particle ...
			{						
				for (j=0;j<TR[level].tr[TRsize0].size;j++) // ... add the others as informers
				{
					rank=TR[level].tr[TRsize0].part[i].ig.size;
					labelt=TR[level].tr[TRsize0].part[j].label;
					TR[level].tr[TRsize0].part[i].ig.label[rank]=labelt;
					TR[level].tr[TRsize0].part[i].ig.size=rank+1;

				}
			}
		}

empty_tribe:

		// Eventually remove empty tribes
		for (i=0;i<TR[level].size;i++)
		{
			if (TR[level].tr[i].size>0) continue; // Not empty
                     i0=i;

				if (i<TR[level].size-1) // Compact the tribe list
				{
					for (j=i;j<TR[level].size-1;j++)
					{
						TR[level].tr[j]=TR[level].tr[j+1];
					}
					i=i-1;
				}
				TR[level].size=TR[level].size-1;
				if (problem[level].printlevel>0)

				{
					printf("\n %i tribes. No particle anymore in tribe %i!",TR[level].size+1,i0);
					printf(" => New tribe number %i",TR[level].size);
				}
				if (TR[level].size==0)
				{
					printf("\n NO PARTICLE ANYMORE!");
					goto the_end;
				}
				goto empty_tribe;
		}

if(linkreorg==1 && TR[level].size>1) // NEW 6.2 Optionnally reorganize links
{
	link_reorg(level);
}

		cycle=0;
		if(problem[level].save==1) save_swarm(f_swarm,level);

	end_adapt:;       //   ---------------------------  ADAPTATION(S)  (end)



 if(rand_hood>0)
 {
      swarm_size=problem[level].init_size+n_add-n_remove;
	  if(circular_hood==1) rand_hood=K_optim(swarm_size,problem[level].eps);
      for(k=0;k<TR[level].size;k++)
      {
         size= TR[level].tr[k].size;
         for (i=0;i<size;i++)
         {
            TR[level].tr[k].part[i].ig.size=rand_hood;

            for (j=0;j<TR[level].tr[k].part[i].ig.size;j++)
            {
               rank=alea(0,swarm_size);
               TR[level].tr[k].part[i].ig.label[j]=TR[level].tr[k].part[rank].label;
            }
         }
      }
 }

         if (problem[level].save==-1)  // Save swarm size and energies
										// Just for information
        {
           Ek=energy_k(level); Ep=energy_p(level);
		   swarm_size=problem[level].init_size+n_add-n_remove;
           fprintf(f_energy,"\n%i %i %i %i %f %f",(int)eval_f_tot[level],swarm_size,n_add,n_remove,Ek,Ep);
        }


		goto do_Iter; // Iterates till a criterion is met

end_by_chance:
		printf("\n CHANCE !");

the_end: // ------------------------------------------------------------------------------------------------------------------- END
	error_cont=offline_error_cont/(eval_f_tot[level]-n_add);
	if (DYN==0)
    error=total_error(best_result.f);
	else error=offline_error/n_change;// Mean error for dynamic optimisation

	if (problem[level].save==-2) // Save info for convergence curve plotting
	{
		fprintf(f_trace_run,"\n%.0f %f %f %i",eval_f_tot[level],error,error_cont,n_add);
	}

ticks2=clock();
clock_tick=ticks2-ticks1; // Processor time needed
duration = (double)(clock_tick/CLOCKS_PER_SEC);
duration_tot=duration_tot+duration;

eval_run[times]=eval_f_tot[level];
eps_run[times][0]=error;
eps_run[times][1]=error_cont;

			// Keep the min. and max. eps. Just for information
			if (error<min_eps) min_eps=error;
			if (error>max_eps) max_eps=error;

if (level>0)  goto end_PSO;

mean_swarm_size=mean_swarm_size/Iter;
Mean_swarm=Mean_swarm+mean_swarm_size;

if (problem[level].printlevel>=0)
{
   if(error<=problem[level].eps)
   {
      printf("\n SUCCESS");
      }
      else printf("\n FAILURE");

  printf ("\n Run %i => eps %f, eps_cont %f,eval: %.0f, duration: %.2f",times+1,error,error_cont,eval_f_tot[level],duration);
  printf(" + %i - %i",n_add,n_remove);
  printf(" (%.0f)",mean_swarm_size);
   display_position(best_result);



//fprintf(f_trace,"%f %f\n",error,eval_f_tot[level]);


// Save run result
// This print should be consistent with the titles written at the very beginning
fprintf (f_run,"\n%i %i %f",problem[level].funct,times+1,problem[level].target);
if(DYN==0)
for (i=0;i<best_result.f.size;i++) fprintf (f_run," %f",best_result.f.f[i]); // f value(s)
else fprintf(f_run," %.12f",error);// For non recursive dynamic optimisation

fprintf(f_run," %.0f %f %i %i",eval_f_tot[level],duration,n_add,n_remove);

for (i=0;i<best_result.p.size;i++) fprintf(f_run," %.12f",best_result.p.x[i]); // position

  if (problem[level].funct==11) // Apple trees
  {
    post= homogen_to_carte(best_result);
    printf("\n   Cartesian coordinates: ");
    for (i=0;i<post.p.size;i++)  printf(" %f",post.p.x[i]);
      fprintf(f_run,"   ");
  for (i=0;i<post.p.size;i++)  fprintf(f_run," %f",post.p.x[i]);
  }
}


if (eval_f_tot[level]>nb_eval_max) nb_eval_max=eval_f_tot[level];
if (eval_f_tot[level]<nb_eval_min) nb_eval_min=eval_f_tot[level];


if(problem[level].printlevel>0) print_N_run(level,best_result);

if(error<total_error(BEST.f)) BEST=best_result; // Best of the best  amongst runs

	times=times+1;
	if (times<problem[level].times) goto times_loop;

  //-----------------------------------------------------------------------------------------
  if (problem[level].printlevel<0) goto end_PSO;

		ticks2=clock();
		mean_eval=zzz/problem[level].times;
		mean_eps=mean_eps/problem[level].times;
		Mean_swarm=Mean_swarm/problem[level].times;
		duration=duration_tot/problem[level].times;

		mean_eval=0;
		mean_eps=0; mean_eps_cont=0;
//		times=problem[level].times;
		for (i=0;i<times;i++)
		{
			mean_eval=mean_eval+eval_run[i];
			mean_eps=mean_eps+eps_run[i][0];
			mean_eps_cont=mean_eps_cont+eps_run[i][1];
		}
		mean_eval=mean_eval/times;
		mean_eps=mean_eps/times;
		mean_eps_cont=mean_eps_cont/times;

		sigma_eval=0;
		for (i=0;i<problem[level].times;i++)
		{
			sigma_eval=sigma_eval+(eval_run[i]-mean_eval)*(eval_run[i]-mean_eval);
		}
		sigma_eval=sqrt (sigma_eval/problem[level].times);


		sigma_eps=0; sigma_eps_cont=0;
		for (i=0;i<problem[level].times;i++)

		{
			sigma_eps=sigma_eps+(eps_run[i][0]-mean_eps)*(eps_run[i][0]-mean_eps);
			sigma_eps_cont+=(eps_run[i][1]-mean_eps_cont)*(eps_run[i][1]-mean_eps_cont);
		}
		sigma_eps=sqrt (sigma_eps/problem[level].times);
		sigma_eps_cont=sqrt (sigma_eps_cont/problem[level].times);

		fprintf(f_synth,"\n");
		fprintf(f_synth,"%2i",problem[level].funct);
		fprintf(f_synth," %2i",problem[level].DD);
		fprintf(f_synth," %4.3f",problem[level].target);
		fprintf(f_synth," %4.10f",problem[level].eps);
		fprintf(f_synth," %i",problem[level].times);
		fprintf(f_synth," %.0f",mean_eval); // Mean number of evaluations
		fprintf(f_synth," %.0f",sigma_eval); // Standard deviation for eval
      z= 1- (double)failure_tot/(double)problem[level].times ;
		fprintf(f_synth," %.3f",z);  // Success rate

		fprintf(f_synth, " %.0f %.0f",nb_eval_min,nb_eval_max);
		fprintf(f_synth, " %4.0f ",Mean_swarm);

		fprintf(f_synth, " %4.10f ",mean_eps);
		fprintf(f_synth," %4.10f",sigma_eps); // Standard deviation for eps
		fprintf(f_synth, " %4.10f ",mean_eps_cont);
		fprintf(f_synth," %4.10f",sigma_eps_cont);

		fprintf(f_synth, " %4.10f ",min_eps);
		fprintf(f_synth, " %4.10f ",max_eps);

		fprintf(f_synth," / Time %4.4f",duration);

		printf("\nMEAN NUMBER OF EVALUATIONS: %.0f",mean_eval); // Mean number of evaluations
		printf(" sigma %.0f",sigma_eval); // Standard deviation
		printf(" [%.0f %.0f ]",nb_eval_min,nb_eval_max);

		printf("\nMean eps: %4.10f",mean_eps);
		printf(" sigma %4.10f",sigma_eps);
		printf(" [%4.10f, %4.10f]",min_eps,max_eps);
		printf("\nMean eps cont.: %4.10f",mean_eps_cont);
		printf(" sigma %4.10f",sigma_eps_cont);


		printf("\nMean swarm:%4.0f ",Mean_swarm);

		printf("\nTOTAL NUMBER OF FAILURES: %i/%i", failure_tot,problem[level].times);
//		if (eval_f_sup>0) printf(" (after %.0f evaluations)",eval_f_sup);
  printf("\n Success rate %.3f", z );


		duration = duration_tot/(float)problem[level].times;
		printf(" / MEAN TIME %4.4f  ",duration);

/*
   if (level==0)    // Display strategies used
   {
      printf("\n strategies used: ");
      for (i=0;i<nb_strateg;i++)
      if (strateg[i]>0)

         printf("\n %i: %.0f times  ",i,strateg[i]) ;
     }
 */

 if (level==0) // Display particle status
 {
   printf("\n\n Status  Count    Rate ");
    z=0; for (i=0;i<Max_status;i++) z=z+status_count[i];
   for (i=0;i<Max_status;i++) printf("\n %i      %.0f    %.0f",i,status_count[i],100.0*status_count[i]/z);
 }

 if (problem[level].printlevel>1) display_swarm(1,level);

 save_swarm(f_swarm,level);
fprintf(f_swarm,"\n%i\n",-1); // End of save swarm file
      fclose(f_swarm);
         fprintf(f_run,"\n%i\n",-1); // End of save run file
         fclose(f_run);

  if (problem[level].nb_f==1) goto end_PSO;

  // Multiobjective problem. Clean up the result file
printf("\n Multiobjective cleaning");
  f_run=fopen("run.txt","r");
  f_run_clean=fopen("run_cleaned.txt","w");

  clean_run(f_run,f_run_clean,problem[level].nb_f,problem[level].DD,level);

end_PSO:
if (problem[level].printlevel>0)
{
   if(MEMO==1) display_memo(level);
}
return best_result;

}


 // ----------------------------------------------------------------------------- ADD_MEMO
void add_memo(struct position x, int level)
{
 /*
   Put a new position in memo, so that memo is sorted by increasing error.
   Note: memo[2] is a global variable
   Useful only if MEMO option is used.
 */

 double error,error1;
  int  i;
  int rank_p;

 error=total_error(x.f);
 if (memo[level].size==0)  // The very first time
 {
       memo[level].x[memo[level].size]=x;
       memo[level].error[memo[level].size]=error;
       memo[level].size=memo[level].size+1;
       rank_p=0;
       goto end;
 }


 // Find the first worst already stored
 for (i=0;i<memo[level].size;i++)
 {
   error1=memo[level].error[i];
   if (error1>=error) {rank_p=i;goto shift;}
 }

 rank_p=memo[level].size;
 goto insert;

 shift:  // Make place for the new position to store

  for (i=memo[level].size;i>rank_p;i=i-1)
  {
     memo[level].x[i]=  memo[level].x[i-1] ;
     memo[level].error[i]=memo[level].error[i-1];
  }



 insert:  // Put the new position
   memo[level].size=memo[level].size+1;
  if (memo[level].size>=Max_Memo-1) memo[level].size=Max_Memo-1;

   memo[level].x[rank_p]=x;
   memo[level].error[rank_p]=error;

   end:
 //printf("\n add_memo (size %i) rank %i, error %f",memo[level].size,rank_p,error);
 //display_memo(level)
;

}

// ----------------------------------------------------------------------------- BEST_INFORMER
struct particle	best_informer(struct particle part,int no_best,int level)
{
/*
   Look for the best informer of the particle.
   If no_best=1, choose it at random.
    2 => pseudo gradient
	0 (default) => really the best
	3 => 0 or 2 at random


*/
int               better;
char			bid;
double			delta1,delta2;
double			dist;
double			error;
double			grad;
int				i;
int				j;

int				j_best;
int				k;
int				k_best;

struct f			f0;
int				labelt;

switch (no_best)
{
case 1:// Choose it at random

	i=alea(0,part.ig.size-1);
    labelt=part.ig.label[i];
	// Find the particle with this label
	for (j=0;j<TR[level].size;j++) // For each tribe ...
		for (k=0;k<TR[level].tr[j].size;k++) // ... for each particle in the tribe
		{
			if (TR[level].tr[j].part[k].label==labelt)
               {
                j_best=j;
                k_best=k;
                break;
               }
		}
	printf("\nERROR. best_informer 1. Can't find particle with label %i",labelt);
	display_swarm(1,level);
	scanf("%c",&bid);
    return TR[level].tr[j].part[k];

case 2: // Pseudo gradient
case2:
	j_best=-1;
	k_best=0;
	grad=0;

	for (i=0;i<part.ig.size;i++) // For each informer
	{
		labelt=part.ig.label[i];

		// Find the particle with this label
		for (j=0;j<TR[level].size;j++) // For each tribe ...
		{
			for (k=0;k<TR[level].tr[j].size;k++) // ... for each particle in the tribe
			{
				if (TR[level].tr[j].part[k].label==labelt) goto compare2;
			}
		}
		printf("\nERROR. best_informer (2). Can't find particle with label %i",labelt);
		printf("\n Particle label %i. ig list:",part.label);
		for (k=0;k<part.ig.size;k++) printf(" %i",part.ig.label[k]);
		display_swarm(1,level);
		scanf("%c",&bid);

		break;


	compare2:
		error=total_error(part.p_i.f);
		delta1=error-total_error(TR[level].tr[j].part[k].p_i.f);
		dist=distance(part.p_i,TR[level].tr[j].part[k].p_i);
		if (dist<=0) continue;


		delta2=delta1/dist;
		if (delta2>grad)
		{
			j_best=j;k_best=k; grad=delta2;
		}
	}
  if(j_best>=0) break;
  // If no better than the particle itself
   //if(j_best<0) return part;  // return the part
   if(j_best<0) goto case0; // Use the classical method
   if(j_best<0) {j_best=0;k_best=0;break;} // Choose particle 0


case 3: // Random 0 or 2

	if (alea(0,1)<0.5) goto case0; else goto case2;

default: // Find really the best
case0:
	j_best=0;
	k_best=0;
	//f0=TR[level].tr[j_best].part[k_best].p_i.f;
	f0.size=TR[level].tr[j_best].part[k_best].p_i.f.size;
	for (i=0;i<f0.size;i++) f0.f[i]=infinite;

/*
printf("\n best_informer. Particle %i. Informers: ",part.label);
for (i=0;i<part.ig.size;i++)
{
	labelt=part.ig.label[i];
	printf(" %i",labelt);
}
*/
	for (i=0;i<part.ig.size;i++) // For each informer
	{
		labelt=part.ig.label[i];
		// Find the particle (tribe j, particle k) with this label
		for (j=0;j<TR[level].size;j++) // For each tribe ...
		{
			for (k=0;k<TR[level].tr[j].size;k++) // ... for each particle in the tribe
			{
				if (TR[level].tr[j].part[k].label==labelt) goto compare;
			}
		}
		printf("\nERROR. best_informer (default). Can't find particle with label %i",labelt);
		printf("\n Particle label %i. ig list:",part.label);
		for (k=0;k<part.ig.size;k++) printf(" %i",part.ig.label[k]);
		display_swarm(1,level);
		scanf("%c",&bid);

	compare:
        better=better_than( TR[level].tr[j].part[k].p_i.f,f0,level);
		if (better!=1) continue;
			j_best=j;
			k_best=k;
			f0=TR[level].tr[j_best].part[k_best].p_i.f;
	}
	break;
} // end switch no_best

// end:
 //printf("\n label %i, k_best %i",part.label,k_best);
return TR[level].tr[j_best].part[k_best];
}


// ----------------------------------------------------------------------------- BETTER_THAN
int	better_than(struct f f1,struct f f2, int level)
{


/* return 1 if f1 is better than f2
          0 if f2 is worse than f2
          2 if they are equivalent
Useful mainly for multiobjective problem
f1 and f2 are supposed to be the same size
*/

 int d;
 double eps;

  if(lexico>0) goto lexico;

 for( d=0;d<f1.size;d++)
 {
    if (f1.f[d]>f2.f[d]) return 0;
 }
 // All f1's f_values are <=
  for( d=0;d<f1.size;d++)
 {
      if (f1.f[d]<f2.f[d]) return 1;  // At least one is <
 }
return 2;    // They are all equal

 lexico: // Lexicographic sort

 eps=problem[level].eps/ f1.size;

   for( d=0;d<f1.size;d++)
 {
	  if (fabs(f1.f[d]-f2.f[d])+eps<MIN(f1.f[d],f2.f[d])) continue;
	  if (f1.f[d]<eps && f2.f[d]<eps) continue;
	 if (f1.f[d]<f2.f[d]-eps)   return 1;
     if (f1.f[d]>=f2.f[d]-eps)   return 0;
 }

 return 2;
}


//------------------------------------------------------------------------------ CLEAN_RUN
void clean_run(FILE *f_run,FILE *f_run_clean, int nb_f,int DD,int level)

{

 /*
   For multiobjective problem, keep only the non dominated points
  f_run: 3 numbers, nb_f values, 4 numbers, position (DD values)
  EOF = -1
 */
 int     better;
 float bid;
 int     d;

  int eof;
  int    i,j,k;
  int    ibid;
  char line[400];
  int    n_clean;
  int    n_pos;
  struct position pos[Max_run];
int      clean[Max_run];


  printf("\nMultiobjective. A moment, please, I clean up the result file");

if (nb_pb_max>1)
{
	printf("\n WARNING. I can clean up the run file for you are running");
	printf("\n           the program on several problems");

}

// Shunt the title line
fgets( line,400,f_run);

// Read the positions
   n_pos=0;
read_pos:
        fscanf(f_run,"%i",&eof);

        if(eof==-1) goto clean_up;

        fscanf(f_run,"%i %f",&ibid,&bid);

        for (d=0;d<nb_f;d++)
        {
           fscanf(f_run," %f",&bid);
           pos[n_pos].f.f[d]=bid;
        }

        fscanf(f_run,"%f %f %f %f",&bid,&bid, &bid,&bid);

        for (d=0;d<DD;d++)
        {

           fscanf(f_run," %f",&bid);
            pos[n_pos].p.x[d]=bid;
        }
        pos[n_pos].f.size=nb_f;
        pos[n_pos].p.size=DD;

         n_pos=n_pos+1;
         goto read_pos;

clean_up:
// Write the titles
for(d=0;d<nb_f;d++)    fprintf(f_run_clean," f%i", d+1);
for(d=0;d<DD;d++)    fprintf(f_run_clean," x%i", d+1);

 n_clean=0;
 // For each result, look if it is dominated by another one
 // If not, write it in the f_run_clean file


 for (i=0;i<n_pos;i++)
 {
       if (i==n_pos-1) goto already;

       // Compare to not already checked positions
       for (j=i+1;j<n_pos;j++)
       {
         better=better_than(pos[j].f,pos[i].f,level);
         if (better==1) goto next_i; // pos(j) is better. Don't keep pos(i)
       }

 already:
         // Compare to already saved positions
        for (j=0;j<n_clean;j++)
        {
           k=clean[j];
           better=better_than(pos[k].f,pos[i].f,level);

           if (better!=0) goto next_i;  // pos(k) is better or equal. Don't keep pos(i)
         }

         clean[n_clean]=i;
         n_clean=n_clean+1;
 next_i:;

 }

 // Write the cleaned result

 for (i=0;i<n_clean;i++)
 {
    j=clean[i];
    fprintf(f_run_clean,"\n");
      for(d=0;d<nb_f;d++)    fprintf(f_run_clean," %f", pos[j].f.f[d]);
      for(d=0;d<DD;d++)    fprintf(f_run_clean," %.12f", pos[j].p.x[d]);
 }
 printf("\n OK, finished. I keep %i non dominated positions. See file run_cleaned.txt",n_clean);
}



// ----------------------------------------------------------------------------- COMPLETE_PART
struct particle complete_part(struct particle par, int option,int level)
 // WARNING: position must be already initialized. Complete the particle
{
struct particle	part;

part=par;

if (option==0) // Completely new particle
{
	part.label=label[level];
	label[level]=label[level]+1;
}

//part.ig.size=1;
//part.ig.label[0]=part.label;

part.x.p.size=problem[level].DD;

part=constrain(part,level); // Keep in the search space

part.x.f.size=problem[level].nb_f;
part.x.f=MyFunction(part.x,problem[level].funct,problem[level].target,level); //evaluation of the position
part.p_i=part.x; // previous best = current position
part.status=0; // no improvement
part.prev_x=part.x;
return part;

}

// ----------------------------------------------------------------------------- HOMOGEN_TO_CARTE
struct position	homogen_to_carte(struct position pos)
{
/* Convert a position given by its homogeneous coordinates into
the same position given by its cartesian coordinates. The homogeneous coordinates
are supposed to be given relatively to the unit points:
(0, ... 0), (1,0,...0) (0,1,0...,0) ... (0,...,0,1)

See, for example, Apple Trees problem

*/
int				d;

struct position post;
double			s;

post.p.size=pos.p.size-1;

s=0;
for (d=0;d<pos.p.size;d++)
{
		s=s+pos.p.x[d];
	}
if (fabs (s)<almostzero)
{
	printf("\n WARNING. Incorrect homogeneous coordinates. Return null point");
	for (d=0;d<post.p.size;d++)

		post.p.x[d]=0;
	return post;
}

for (d=0;d<post.p.size;d++)
	post.p.x[d]=pos.p.x[d+1]/s;

return post;
}




// ----------------------------------------------------------------------------- INIT_PARTICLE
struct particle	init_particle(int option,int level,struct particle part0)
{ // Initialisation of a particle
int				d;
double			dist_min;
int				loop;
int				loop_max;

int				random; // See below the different kinds of random generation

int				place;
int				rank;
struct particle part={0};

//random=alea(0,3);  // 0,1,2,3
random=alea(0,2);
//if (problem[level].nb_f>1)
//   random=2*alea(0,1); // 0 or 2    (mainly useful for multiobjective problem)
//else
 //  random=0;
//random=alea(0,1); if(random==1) random=2;

part.x.p.size=problem[level].DD;

if(option>0) part=part0; // Re-init position, but no the links
if (option==2) goto complete; // Just recompute the position

dist_min=0;
loop=0;
loop_max=problem[level].DD;

switch (random)
{
case 0:	// Random generation anywhere in the search space
for( d = 0;d<problem[level].DD;d++)  //for each dimension
{
	part.x.p.x[d] = alea_float(problem[level].H.min[d],problem[level].H.max[d]);
}


break;

case 1: // Random generation on an edge of the search space

for( d = 0;d<problem[level].DD;d++)
{

	place=alea(0,2);

	switch (place)
	{
	case 0:
		part.x.p.x[d] = problem[level].H.min[d];
		break;
	case 1:
		part.x.p.x[d] = problem[level].H.max[d];
		break;
	case 2:
		part.x.p.x[d] = (problem[level].H.min[d]+problem[level].H.max[d])/2;
		break;
	}
}
break;

case 2: // On the frontier
	for( d = 0;d<problem[level].DD;d++)
	{
		part.x.p.x[d] = alea_float(problem[level].H.min[d],problem[level].H.max[d]);
	}
	rank=alea(0,part.x.p.size-1);
	place=alea(0,1);
	if (place==0) 	part.x.p.x[rank]=problem[level].H.min[d];
	else part.x.p.x[rank]=problem[level].H.max[d];
	break;

case 3: // Inside the memory swarm
//printf("\n init_particle %f %f",Xmin.x[0],Xmax.x[0]);

for( d = 0;d<problem[level].DD;d++)  //for each dimension
{
	part.x.p.x[d] = alea_float(Xmin.x[d],Xmax.x[d]);
}

}

complete:
part=complete_part(part,option,level);
//display_position(part.x);

if (problem[level].printlevel>1)
   printf("\n init => %f",total_error(part.x.f));

// For some strategies, more features are needed
for( d = 0;d<problem[level].DD;d++)  //for each dimension
{
     part.mod[d]=0; // Has not been modified (for Taboo option)
     part.v[d]=0; // Or random, if you want
  }
return part;
}

// ----------------------------------------------------------------------------- LINK_REORG
void link_reorg(int level)
{

/* Reorganise information links between tribes
The tribe list TR[level] is a global variable

Note 1: symmetry of the information graph is not garanteed
Note 2: connectivity of the information graph is not garanteed
Note 3:  this is just a "research" option. Does not seem to be very efficient
*/
int					better;
struct	particle	best;
double				dist;
double				dist_far;
double				dist_near;
struct	particle	farest;
int					far_ig_rank;
int					i,i2,i3,i4;
int					i_best;
int					i_near;
int					inform;
int					label;
struct	particle	nearest;
int					not_inform;
struct	particle	part;
int					t,t2,t3;
int					t_near;


for(t=0;t<TR[level].size;t++) // For each tribe
{
// ... find the best particle B
	best=TR[level].tr[t].part[0]; i_best=0;

	for (i=1;i<TR[level].tr[t].size;i++)
	{
		better=	better_than(TR[level].tr[t].part[i].p_i.f,best.p_i.f,level);
		if(better==1) {best=TR[level].tr[t].part[i]; i_best=i;}
	}

	// Find a particle that is not an informer
	for(t2=0;t2<TR[level].size;t2++)
	{
		if (t2==t) continue; // All particles of the tribe t are informers for best

		label=TR[level].tr[t].part[i].label;


		for (i=0;i<best.ig.size;i++)
		{
			if (best.ig.label[i]!=label)
			{nearest=TR[level].tr[t].part[i];goto loop_swarm;}
		}

		printf("\nWARNING.link_reorg. Can't find a non informer for particle label %i",best.label);

		goto end;

	}

loop_swarm:
// ... loop on all particles of the swarm
	farest=best; far_ig_rank=-1;
	dist_far=0; // Distance farest informer<->best
	dist_near=infinite; // Distance nearest non informer <-> best
	not_inform=0;

	for(t2=0;t2<TR[level].size;t2++)
	{
		for (i=0;i<TR[level].tr[t2].size;i++)
		{

			label=TR[level].tr[t2].part[i].label;
			if(label==best.label) continue;

			// Find the particle whose label is "label"
			for(t3=0;t3<TR[level].size;t3++)
			{
				for (i3=0;i3<TR[level].tr[t3].size;i3++)
				{
					if(label!=TR[level].tr[t3].part[i3].label) continue;
						part=TR[level].tr[t3].part[i3];
						goto links;
				}
			}
			printf("\nERROR. link_reorg. Can't find a particle with label %i",label);

			display_swarm(1,level);
			goto end;

links:
			// Compute the distance
			dist=distance(best.p_i,part.p_i);
			inform=0;
//   if it is an informer, look for the farest one F

			for (i2=0;i2<best.ig.size;i2++)
			{
				if (best.ig.label[i2]==label)
				{
					inform=1;
//printf("\n dist,dist_far %f %f",dist,dist_far);
//printf("\n  label %i",best.label);display_position(best.p_i);
//printf("\n   label %i",part.label); display_position(part.p_i);

					if(dist>=dist_far)
					{farest=part;dist_far=dist; far_ig_rank=i2;}
				}
			}

			if(inform==1) continue;

			//   if it is not an informer, look for the nearest one N
			not_inform=1;
//printf("\n %f %f",dist,dist_near);
			if(dist<dist_near)
			{nearest=part;dist_near=dist;t_near=t3;i_near=i3;
//printf("\n t3,i3 %i %i", t3,i3);
			}
		} // end loop i on particles

	} // end loop t2 on tribes

//printf("\n not inform %i",not_inform);

if(far_ig_rank<0)
{
	printf("\nWARNING. link_reorg. No farer informer for particle %i \n",best.label);
	printf("\n ig size %i/ List:",best.ig.size);
	for (i3=0;i3<best.ig.size;i3++) printf(" %i",best.ig.label[i3]);
	goto end;
}
if(not_inform==1) // If there exist a non informer particle for "best"
{
	//   remove link B<-F and replace it by the link B<-N

//printf("\n best %i",best.label);
//printf(" ig %i",far_ig_rank);
//printf("  %i => %i",best.ig.label[far_ig_rank],nearest.label);

	TR[level].tr[t].part[i_best].ig.label[far_ig_rank]=nearest.label;

// Add symmetrical link
//printf("\n t_near i_near %i %i, ig size %i",t_near,i_near, nearest.ig.size);

	for (i4=0;i4<nearest.ig.size;i4++)
		if (nearest.ig.label[i4]==best.label) goto end_not_inform;
	// best.label not already on the informer list of nearest
		nearest.ig.label[nearest.ig.size]=best.label;
		nearest.ig.size=nearest.ig.size+1;
		TR[level].tr[t_near].part[i_near]=nearest;
end_not_inform:;
}
} // end loop t on tribes

//printf("\n");

end:;
}

// ----------------------------------------------------------------------------- LOCAL_IMPROVE
int	local_improv(struct tribe tr,int level)
{
// Count how many particles in the tribe have improved their position
int better;
int	i;
int	n;

n=0;
for (i=0;i<tr.size;i++)
{
better=better_than(tr.part[i].p_i.f, tr.part[i].prev_x.f,level);
   if (better==1) n=n+1;
}
return n;

}

// ----------------------------------------------------------------------------- MOVE_PARTICLE
struct particle	move_particle(struct particle part,struct particle partg,int level, int gener)
{ //  ********* CORE OF THE ALGORITHM ************************************************

 // * after a move option means "seems interesting"
// ?  means "questionable"

int                     better;
double				c0,c1,c2,c3;
struct	position	center_i,center_g;
int					choice;
double			coeff_rho;
int					d;
int					DD;
double               delta;
double               df;
int         dist;
double               error,error_g,error_i;
int           k;
double         lambda;
double         mu;
double          nonunif;

struct particle		partt;
double				pgd;
double				pid;
struct   position       pivot;
double				pr;
double        r;
struct	vector		rand_x;
struct	vector		rand_i,rand_g;
double				rho,rho1,rho2,rho3,rho_x,rho_i,rho_g;

double         sigma, sigma_d;
int               sign=1;
int                  status;
int					strategy;
int               sub_case;
double               xd;
double        vd;
double         z,zz;
 if(problem[level].printlevel>1)
 {
    if (gener==0)
	 printf("\n move_particle %i, status %i",part.label,part.status);
	else
		printf("\n move_particle. Generate a new particle");
	}


DD=part.x.p.size;
partt=part;

if (gener==1) {choice=1;goto strategy0;} // Not really a move but a generation


partt.prev_x=partt.x;

error=total_error(part.x.f);
error_i=total_error(part.p_i.f);
error_g=total_error(partg.p_i.f);
if(error_g<=0) return partt; // Solution has already been found!

if (part.status<=11) {status=0;   goto switch_s;}

if (part.status<35) {status=1;   goto switch_s;}
if (part.status<53) {status=2;   goto switch_s;}
status=3;

switch_s:

switch (status)
{
case 0:   // Bad particle
	strategy=(int)strategies[0][0];
      nonunif=strategies[1][0];
	break;

case 1:   // Neither bad nor good
	strategy=(int)strategies[0][1];
      nonunif=strategies[1][1];
	break;
case 2:     // Good particle
	strategy=(int)strategies[0][2];
      nonunif=strategies[1][2];
	break;
case 3:       // Very good particle
	strategy=(int)strategies[0][3];
      nonunif=  strategies[1][3];
	break;
}

status_count[status]= status_count[status]+1; // Just for information

switch (strategy)
{
case -1: // Pure random

for( d = 0;d<DD;d++)  //for each dimension
{
	partt.x.p.x[d] = alea_float(problem[level].H.min[d],problem[level].H.max[d]);

}
break;
 //------------------------------------------------------------------------------------------------- CASE 0
case 0: // ? Around p_i OR around p_g

// Method 0. Very good for some functions, very bad for some others

c1=1/error_i;
c2=1/error_g;
c3=c1+c2;
c1=c1/c3;
c2=c2/c3;
pr=alea_float(0,1);
if (pr<=c1) choice=0; else choice=1;    // More probably around p_g

//if (pr<=0.5) choice=0; else choice=1; // Same probability
//choice=1; // TEST only around p_g

strategy0:

rho1=distance(part.x,part.p_i);
rho2=distance(part.x,partg.p_i);
rho3=distance(part.p_i,partg.p_i);



if (choice==0)  // A point around p_i
{
      if (rho3>rho1) rho=rho3; else rho=rho1;
      rand_x=rand_in_hypersphere(DD, rho,nonunif);
      for (d = 0;d<DD;d++)
	{
		pid=part.p_i.p.x[d];
		partt.x.p.x[d]  = rand_x.x[d]+pid;
	}
}
else  // A point around p_g
{
   if (rho3>rho2) rho=rho3; else rho=rho2;


      rand_x=rand_in_hypersphere(DD, rho,nonunif);
   for (d = 0;d<DD;d++)
	{
		pgd=partg.p_i.p.x[d];
		partt.x.p.x[d]  = rand_x.x[d]+pgd;
	}
}

if (gener==1) return partt;
break;
  //------------------------------------------------------------------------------------------------- CASE 1

case 1: // * Pivot. Mixture "around p_i" "around p_g". The best compromise?
//  case1:
	lambda=1;
  case1a:
rho=distance(part.p_i,partg.p_i);
rho_i=rho; rho_g=rho;


//z=d_min(part.p_i,level); if (z<rho) rho_i=z;
//z=d_min(partg.p_i,level); if (z<rho) rho_g=z;

rand_i=rand_in_hypersphere(DD, rho_i,nonunif); // Prepare random points in hypersphere
rand_g=rand_in_hypersphere(DD, rho_g,nonunif);

    c1=1/error_i;
	c2=1/error_g;
	c3=c1+c2;

	c1=c1/c3;
	c2=c2/c3;

	for (d = 0;d<DD;d++) // for each dimension // HERE THE PARTICLE MOVES *****
	{
		pid=part.p_i.p.x[d];
		pgd=partg.p_i.p.x[d];
		partt.x.p.x[d]  = c1*(rand_i.x[d]+pid) + c2*(rand_g.x[d]+pgd);

		partt.x.p.x[d]=partt.x.p.x[d]*lambda;
	}
	break;

 //------------------------------------------------------------------------------------------------- CASE 2
case 2: // ? Mixture "around x" "around p_i" "around p_g"


//Method 2
c0=1/error;
c1=1/error_i;
c2=1/error_g;


c3=c0+c1+c2;
c0=c0/c3;
c1=c1/c3;
c2=c2/c3;

rho1=distance(part.x,part.p_i);
rho2=distance(part.x,partg.p_i);
rho3=distance(part.p_i,partg.p_i);

if (rho2<rho1) rho_x=rho2; else rho_x=rho1;
if (rho1<rho3) rho_i=rho1; else rho_i=rho3;
if (rho3<rho2) rho_g=rho3; else rho_g=rho2;


rand_x=rand_in_hypersphere(DD, rho_x,nonunif);
//printf("\n %f",rho);


rand_i=rand_in_hypersphere(DD, rho_i,nonunif);

//printf("\n %f",rho);

rand_g=rand_in_hypersphere(DD, rho_g,nonunif);
//printf("\n %f\n",rho);

for (d = 0;d<DD;d++) // for each dimension // HERE THE PARTICLE MOVES *****
{
	partt.x.p.x[d]  = c0*(part.x.p.x[d]+rand_x.x[d])+c1*(part.p_i.p.x[d]+rand_i.x[d])+c2*(partg.p_i.p.x[d]+rand_g.x[d]);
}
break;

  //------------------------------------------------------------------------------------------------- CASE 3
case 3: // ? Try to find promising areas (center and radius) and choose one at random
// Method 3.
center_i.p.size=DD;

center_g.p.size=DD;
c0=error-error_i;

if (c0<almostzero)
{
	center_i=part.p_i;
	rho_i=distance(part.p_i,partg.p_i);
}
else
{
	c0=error/c0;
	for (d=0;d<DD;d++)
		center_i.p.x[d]=part.x.p.x[d]+c0*(part.x.p.x[d]-part.p_i.p.x[d]);
	rho_i=distance(part.p_i,center_i);
}

c0=error_i-error_g;


if (c0<almostzero)
{
	center_g=partg.p_i;
	rho_g=distance(part.p_i,partg.p_i);

}
else
{
	c0=error_i/c0;
	for (d=0;d<DD;d++)
		center_g.p.x[d]=part.p_i.p.x[d]+c0*(part.p_i.p.x[d]-partg.p_i.p.x[d]);
		rho_g=distance(partg.p_i,center_g);
}


c1=1/error_i;
c2=1/error_g;
c3=c1+c2;
c1=c1/c3;
c2=c2/c3;
pr=alea_float(0,1);
if (pr<=c1) choice=0; else choice=1;


if (choice==0)
{
	rand_i=rand_in_hypersphere(DD, rho_i,nonunif);
	for (d = 0;d<DD;d++) // A point around center_i

	{
		pid=center_i.p.x[d];
		partt.x.p.x[d]  = rand_i.x[d]+pid;
	}
}
else
{
	rand_g=rand_in_hypersphere(DD, rho_g,nonunif);
	for (d = 0;d<DD;d++) // A point around center_g

	{
		pgd=center_g.p.x[d];
		partt.x.p.x[d]  = rand_g.x[d]+pgd;
	}
}
break;

  //------------------------------------------------------------------------------------------------- CASE 4
case 4: // ? Mixture "around p_i" "around p_g" and "particle previous move"

rho=distance(part.p_i,partg.p_i);   // Same radius for the two hyperspheres
rho_i=rho;
rho_g=rho;


c1=error_i;
c2=error_g;
c0=(c1-c2)/(c1+c2);  // Should, on the whole, decrease during the process

case4a:
    	c1=1/error_i;
	c2=1/error_g;
	c3=c1+c2;
	c1=c1/c3;
	c2=c2/c3;

rand_i=rand_in_hypersphere(DD, rho_i,nonunif); // Prepare random points in hypersphere
rand_g=rand_in_hypersphere(DD, rho_g,nonunif);

	for (d = 0;d<DD;d++) // Prepare a point "between" p_i and p_g
	{
		pid=part.p_i.p.x[d];
		pgd=partg.p_i.p.x[d];
		partt.x.p.x[d]  = c1*(rand_i.x[d]+pid) + c2*(rand_g.x[d]+pgd);

	}

 //case4b:
 	for (d = 0;d<DD;d++)  // Add a contribution of the "velocity"
	{
		partt.x.p.x[d]  =  partt.x.p.x[d] +sign*c0*( part.x.p.x[d]-part.prev_x.p.x[d]);
	}
	break;
  //------------------------------------------------------------------------------------------------- CASE 5
   case 5: // ? Mixture "around p_i" "around p_g" and "particle previous move"
 // Like case 4, but radiuses may be diffferent, and bigger
// case5:
rho1=distance(part.x,part.p_i);


rho2=distance(part.x,partg.p_i);

rho3=distance(part.p_i,partg.p_i);

//if (rho1>rho3) rho_i=rho1; else rho_i=rho3;
//if (rho3>rho2) rho_g=rho3; else rho_g=rho2;

rho_i=MAX(rho1,rho3);
rho_g=MAX(rho2,rho3);

c1=error_i;
c2=error_g;
c0=(c1-c2)/(c1+c2);  // Should, on the whole, decrease during the process

goto case4a;

   //------------------------------------------------------------------------------------------------- CASE 9  (HYPERPARALLELEPIDS)
   case 9: //* "Classical PSO", for comparison. Mixture "around p_i" "around p_g" and "particle previous move"
            // with hyperparallelepids


	for (d = 0;d<DD;d++)
	{
		pid=part.p_i.p.x[d];
		pgd=partg.p_i.p.x[d];
    xd=part.x.p.x[d];
    partt.x.p.x[d]=xd+khi*(xd-part.prev_x.p.x[d]) + alea_float(0,cmax/2)*(pid-xd) + alea_float(0,cmax/2)*(pgd-xd);
   }

	break;


    //------------------------------------------------------------------------------------------------- CASE 10
   case 10:  //?? Keep moving. Neighbours are not taken into account.
	   sign=1;
case10:
 c1=error_i;
c2=error_g;
c0=(c1-c2)/(c1+c2);  // Should, on the whole, decrease during the process


	for (d = 0;d<DD;d++)
	{
            xd=part.x.p.x[d];
            partt.x.p.x[d]  =xd+sign*c0*(xd-part.prev_x.p.x[d]) ;
            //partt.x.p.x[d]  =xd+alea_float(0,1)*(xd-part.prev_x.p.x[d]) ;
	}
   break;

     //------------------------------------------------------------------------------------------------- CASE 11
    case 11:   //?? Go back.   Neighbours are not taken into account.
    sign=-1;
    goto case10;
    break;
   //------------------------------------------------------------------------------------------------- CASE 12
case 12: //* Pivot+noise. Mixture "around p_i" "around p_g" + noise on position

    //6.1 Add decreasing noise (the same on each element)
        z=error_i;
        zz=error_g;
        c0=(z-zz)/(z+zz);
		lambda=1+alea_normal(0,c0);
		goto case1a;
	break;
 // ------------------------------------------------------------------------------------------------ CASE 13

 //------------------------------------------------------------------------------------------------- CASE 14

case 14: // Mixture "around p_i" "around p_g" , using Direct Gaussian distribution
rho=distance(part.p_i,partg.p_i);

rand_i=rand_in_hypersphere(DD, rho,-fabs(nonunif)); // Prepare random points in hypersphere

    c1=1/error_i;


	c2=1/error_g;
	c3=c1+c2;
	c1=c1/c3;
	c2=c2/c3;

	for (d = 0;d<DD;d++) // for each dimension // HERE THE PARTICLE MOVES *****
	{
		pid=part.p_i.p.x[d];
		pgd=partg.p_i.p.x[d];
		partt.x.p.x[d]  = c1*pid + c2*pgd+rand_i.x[d];
	}

	break;

  //------------------------------------------------------------------------------------------------- CASE 15
case 15: // Mixture "around p_i" "around p_g" , "around x",
            //using Direct Gaussian distribution (pivot method)
        mu=0.91;
     // Compute "gradients"
   rho=distance(part.x,part.p_i);
   rho_i=  distance(part.x,part.p_i);
   rho_g=distance(part.p_i,partg.p_i);
    if (rho>0) c1=(error-error_i)/rho;   else c1=0;   // "gradient" x-p
    if (rho_i>0) c2=(error-error_g)/rho_i;     else c2=0;        // "gradient" x-g
    if (rho_g>0) c3=(error_i-error_g)/rho_g;   else c3=0;     // "gradient p-g

   // Choose the pivot   (p or g)
      if (c3>=c1 && c3>=c2)
      {
       pivot=partg.p_i;
       rho=rho_g;
        //delta=1+2*(part.p_i.f.f[0] -partg.p_i.f.f[0]);

         sub_case=3;
       goto next_case15;
      }

        if (c2>=c1 && c2>=c3)
        {
         pivot=partg.p_i;
         rho=rho_i;
        // delta=1+2*(part.x.f.f[0] -partg.p_i.f.f[0]);
          sub_case=2;

         goto next_case15;
        }

         pivot=part.p_i;

         //delta=1+2*(part.x.f.f[0] -part.p_i.f.f[0]);
          sub_case=1;

next_case15:
//printf("\n rho %f",rho);
   if(rho<=0) {partt.x=part.x; break; } // No move

             lambda=1;

            for (d = 0;d<DD;d++)
            {
              switch (sub_case)
              {
              case 1:

              delta=fabs(part.x.p.x[d]-part.p_i.p.x[d]);
              break;
              case 2:
              delta=fabs(part.x.p.x[d]-partg.p_i.p.x[d]);
              break;
              case 3:
              delta=fabs(part.p_i.p.x[d]-partg.p_i.p.x[d]);
              break;
              }
              if (delta>0)
              lambda=lambda*delta;

            }
             sigma=mu*pow(rho/lambda,(double)(1/DD));

      for (d = 0;d<DD;d++) // for each dimension // HERE THE PARTICLE MOVES *****
	{
            switch (sub_case)
              {
              case 1:

               sigma_d=sigma*fabs(part.x.p.x[d]-part.p_i.p.x[d]);
              break ;
              case 2:
                sigma_d=sigma*fabs(part.x.p.x[d]-partg.p_i.p.x[d]);
              break;
              case 3:
                 sigma_d=sigma*fabs(part.p_i.p.x[d]-partg.p_i.p.x[d]);
              break;
              }
         //   delta=fabs(alea_normal(0,sigma_d));


// printf(" sigma_d %f  delta %f " ,sigma_d,delta);
		//partt.x.p.x[d]  = pivot.p.x[d]+delta;

         delta=fabs(alea_normal(0, fabs(part.p_i.p.x[d]-partg.p_i.p.x[d])));  // Test
           partt.x.p.x[d]  = (part.p_i.p.x[d]+part.p_i.p.x[d])/2+delta;      // Test
	}
               break;


  //------------------------------------------------------------------------------------------------- CASE 16
case 16: // Mixture "around p_i" "around p_g" , "around x",
            //pivot method (using Direct Gaussian distribution if nonunif<0)

c1=1/error;
c2=1/error_i;
c3=1/error_g;
c0=c1+c2+c3;
c1=c1/c0;
c2=c2/c0;
c3=c3/c0;

   for (d = 0;d<DD;d++) pivot.p.x[d]=c1*part.x.p.x[d]+c2*part.p_i.p.x[d]+c3*partg.p_i.p.x[d];

 pivot.p.size=DD;
  rho=distance(pivot,partg.p_i);
   //rho_i=  distance(part.x,part.p_i);
   //rho_g=distance(part.p_i,partg.p_i);
   rand_g=rand_in_hypersphere(DD, rho,nonunif);


    for (d = 0;d<DD;d++)
    {
       partt.x.p.x[d]  =pivot.p.x[d]+rand_g.x[d];

    }


            break;

   //------------------------------------------------------------------------------------------------- CASE 17

 case 17:  // Simplified mixture (x,p_i,p_g)

rho1=distance(part.x,part.p_i);
rho2=distance(part.x,partg.p_i);
if (rho2<rho1) rho2=rho1;

rand_x=rand_in_hypersphere(DD, rho2,nonunif);

for (d = 0;d<DD;d++) // for each dimension // HERE THE PARTICLE MOVES *****

{

	partt.x.p.x[d]  = part.x.p.x[d]+rand_x.x[d];
}

 break;

 //------------------------------------------------------------------------------------------------- CASE 18
 case 18:  // Pivot method, with hyperspheres
  pivot=pivot_choice(level);
  rho=distance(part.p_i,pivot);    //
  rand_x=rand_in_hypersphere(DD, rho,nonunif);

for (d = 0;d<DD;d++) // for each dimension // HERE THE PARTICLE MOVES *****
{
   partt.x.p.x[d]  = pivot.p.x[d]+rand_x.x[d];
}

  break;

  //------------------------------------------------------------------------------------------------- CASE 19
case 19: // Pivot method with mixture "around p_i" "around pivot"
  pivot=pivot_choice(level);
rho=distance(part.p_i,pivot);

rand_i=rand_in_hypersphere(DD, rho,nonunif); // Prepare random points in hypersphere
rand_g=rand_in_hypersphere(DD, rho,nonunif);
 error_g=total_error(pivot.f);

     	c1=1/error_i;
	c2=1/error_g;
	c3=c1+c2;

	c1=c1/c3;
	c2=c2/c3;

	for (d = 0;d<DD;d++) // for each dimension // HERE THE PARTICLE MOVES *****
	{
		pid=part.p_i.p.x[d];
		pgd=pivot.p.x[d];
		partt.x.p.x[d]  = c1*(rand_i.x[d]+pid) + c2*(rand_g.x[d]+pgd);
	}
	break;
 //----------------------------------------------------------------------------------- CASE 20
 case 20: // Semi-hypersphere in direction x(t-1)->x(t), centered on x
rho=distance(part.prev_x,part.x);
rand_x=rand_in_hypersphere(DD, rho,nonunif);

	for (d = 0;d<DD;d++)
	{
            xd=part.x.p.x[d];
            if (rand_x.x[d]*(xd-part.prev_x.p.x[d])<0)  rand_x.x[d]=0;
            partt.x.p.x[d]  =xd+rand_x.x[d];
	}
    break;
  //----------------------------------------------------------------------- CASE 21  (HYPERPARALLELEPIDS)

   case 21: // * For difficult problems (combinatorial, for example)
			// Gaussian on each dimension


	c1=2/0.97725; // With Gaussian distibrution,  0.97725 is the
	              //probability to have the result in [mean +- 2*standard_dev]

	c0=1/(c1-1+sqrt(c1*c1-2*c1));  // Constriction

	for (d = 0;d<DD;d++)
	{
		pid=part.p_i.p.x[d];
		pgd=partg.p_i.p.x[d];
		xd=part.x.p.x[d];

	rho1=fabs(pid-xd);
	rho2=fabs(pgd-xd);
	partt.x.p.x[d]  =xd+c0*((xd-part.prev_x.p.x[d]) + alea_normal(pid-xd,rho1/2)+ alea_normal(pgd-xd,rho2/2));

   }

	break;

  //---------------------------------------------------------------------- CASE 22  (HYPERPARALLELEPIDS)
   case 22: // For difficult problems (combinatorial?)
			// Intervals
	  c1=1/error_i;
	c2=1/error_g;
	c3=c1+c2;

	c1=c1/c3;
	c2=c2/c3;


	for (d = 0;d<DD;d++)
	{


		pid=part.p_i.p.x[d];
		pgd=partg.p_i.p.x[d];
		rho=fabs(pid-pgd);

		partt.x.p.x[d]  = c1*(rand_i.x[d]+pid) + c2*(rand_g.x[d]+pgd);
		partt.x.p.x[d]  = c1*alea_float(pid-rho,pid+rho) + c2*alea_float(pgd-rho,pgd+rho);

   }
	break;


case 23: // * Ellipso�dal positive sectors
	// TO DO
	// replace by the one of OEP5, which doesn't need the hard coded c1
//case_23:

	rho=coeff_S_C*cmax;
	rand_i=rand_in_hypersphere(DD, rho,nonunif);
	rand_g=rand_in_hypersphere(DD, rho,nonunif);

	 for (d=0;d<DD;d++)
     {
		pid=part.p_i.p.x[d];
		pgd=partg.p_i.p.x[d];
		xd=part.x.p.x[d];
		 partt.x.p.x[d]=xd+khi*(xd-part.prev_x.p.x[d])+fabs(rand_i.x[d])*(pid-xd);
        partt.x.p.x[d]=partt.x.p.x[d]+fabs(rand_g.x[d])*(pgd-xd);
     }
break;

//-------------------------------------------------------------------- CASE 21  (HYPERPARALLELEPIDS)
   case 24: // Gaussian on each dimension
	coeff_rho=1;
	c0=1;
	c1=1;
	c2=1;

	for (d = 0;d<DD;d++)
	{
		pid=part.p_i.p.x[d];
		pgd=partg.p_i.p.x[d];
		xd=part.x.p.x[d];

	rho1=fabs(pid-xd);
	rho2=fabs(pgd-xd);
	partt.x.p.x[d]  =c0*xd+ alea_normal(pid-xd,c1*rho1)+ alea_normal(pgd-xd,c2*rho2);
   }
	break;

   case 25: // Gaussian on each dimension
	       // Just "local" search around the best position
//	choice=alea(0,1);
	choice=0;
	if (choice==0)

		for (d = 0;d<DD;d++)
		{
			pgd=partg.p_i.p.x[d];
			xd=part.x.p.x[d];
			rho2=fabs(pgd-xd);
			partt.x.p.x[d]  =pgd+ alea_normal(pgd-xd,rho2);  // a
		}
	else
		for (d = 0;d<DD;d++)
		{
			pgd=partg.p_i.p.x[d];
			xd=part.x.p.x[d];

			rho2=fabs(pgd-xd);
			partt.x.p.x[d]  =pgd+ alea_normal(0,rho2); // b
		}

	break;

	case 26: // Uniform on each dimension
	       // Just "local" search around the best position

	for (d = 0;d<DD;d++)
	{
		pgd=partg.p_i.p.x[d];
		xd=part.x.p.x[d];
	rho2=fabs(pgd-xd);
	partt.x.p.x[d]  =pgd+ (2*alea(0,1)-1)*alea_float(0,rho2);
   }
	break;

  case 100: //  KE-C Method (Jim Kennedy and Russ Eberhart method, improved by Maurice Clerc)
	c0=1;
	c1=2.0;
	c2=c1;

	for (d=0;d<DD;d++)
    {
        xd=part.x.p.x[d];
        pid=part.p_i.p.x[d];
        pgd=partg.p_i.p.x[d];
        vd=part.v[d]; // "velocity"
        c3=c0*vd + alea_float(0,c1)*(pid-xd);  // c0*v + rand(0,c1)*(p-x)
        c3=c3+alea_float(0,c2)*(pgd-xd); //  + rand(0,c2)*(g-x)
                                      //Note that g is either the local best or the global best
                                      // depending on the neighbourhood size
        partt.v[d]=(int)c3; // New integer velocity
        c3=  xd+ partt.v[d]; // Add velocity
        c3=1/(1+exp(-c3)); // Keep position in ]0,1[

        pr=alea_float(0,1); // Probabilistic binary choice
        if (pr<c3) 	partt.x.p.x[d]=1; else 	partt.x.p.x[d]=0;
    }
    break;

   case 101: //  KE-C with taboo
	c0=1;
	c1=2.0;
	c2=c1;

	for (d=0;d<DD;d++)
    {
        if(part.mod[d]>0) continue; // "Taboo". Modify only if it has not just been modified
        xd=part.x.p.x[d];
        pid=part.p_i.p.x[d];
        pgd=partg.p_i.p.x[d];
        vd=part.v[d]; // "velocity"
        c3=c0*vd + alea_float(0,c1)*(pid-xd);  // c0*v + rand(0,c1)*(p-x)
        c3=c3+alea_float(0,c2)*(pgd-xd); //  + rand(0,c2)*(g-x)
                                      //Note that g is either the local best or the global best
                                      // depending on the neighbourhood size
        partt.v[d]=(int)c3; // New integer velocity
        c3=  xd+ partt.v[d]; // Add velocity
        c3=1/(1+exp(-c3)); // Keep position in ]0,1[

        pr=alea_float(0,1); // Probabilistic binary choice
        if (pr<c3) 	partt.x.p.x[d]=1; else 	partt.x.p.x[d]=0;
    }


    for (d=0;d<DD;d++)
    {
      if (partt.x.p.x[d]!=part.x.p.x[d])  partt.mod[d]=1; // Has been modified
          else partt.mod[d]=0; // Has not been modified
      }
    break;

  case 111: // Binary. Pivot method
   //from an idea of P. Serra, 1997
   // (seriously) modified and adapted to binary problems by M. Clerc, 2004
	// Works pretty well on some problems .. and pretty bad on some others
  // Note that there is no velocity, and the current position is not used
  // It works better when information links are modified at random
  // if there has been no improvement after the previous iteration

 partt.x=partg.p_i; // Best previous position g of the best informer

 // (for D>=2). Simplified formula.  Note that "dist" is an integer
 dist=(int)log(DD);
 if (dist<3) dist=3; if (dist>DD) dist=DD;  // 111a.

  r=alea_float(1,dist); //111. Note that r is never equal to dist
  //r=1+alea(1,dist-1);  //111b. Theoretically equivalent
 // r=1+alea(1,dist); //111c. A bit larger
 // r=alea(1,dist-1);  // 111d.   More uniform , on dist k values
                              // Sometimes far better but also sometimes far worse

 for (k=0;k<r;k++) // Define a position around g
 //    Note that even if there at least two k values
 //  there is a small probability (1/D^k)
 // that just one bit is modified
  {
    d=alea(0,DD-1);
    partt.x.p.x[d]=1- partt.x.p.x[d];

  }
 break;

 case 116: // Binary. Pivot method with "taboo"
 partt.x=partg.p_i; // Best previous position g of the best informer

 dist=(int)log(DD);
 if (dist<3) dist=3; if (dist>DD) dist=DD;  // 111a.

  r=alea_float(1,dist); //11. Note that r is never equal to dist

 for (k=0;k<r;k++)
  {
    d=alea(0,DD-1);
    if(part.mod[d]<1) // "Taboo". Modify only if it has not just been modified
      {
        partt.x.p.x[d]=1- partt.x.p.x[d];
      }
  }

    for (d=0;d<DD;d++)
    {
      if (partt.x.p.x[d]!=part.x.p.x[d])  partt.mod[d]=1; // Has been modified
          else partt.mod[d]=0; // Has not been modified
      }
 break;




 //------------------------------------------------------------------------------------------------- default
default:
   printf("\n ERROR. move_particle. Non defined strategy value %i",strategy);
     break;
}   // End switch(strategy)

// Adjustements -------
if (problem[level].printlevel>2)
{
	printf("\n  move_particle. Adjustements");
	display_position(partt.x);

}

partt=constrain(partt,level); // Keep it in the search space

// End od adjustements ------

 // ********** evaluation of the position
partt.x.f=MyFunction(partt.x,problem[level].funct,problem[level].target,level);

// Save best position
better=better_than(partt.x.f,partt.p_i.f,level);
if( better==1)
{
	partt.p_i = partt.x;
//	nb_improv=nb_improv+1;
}

//display_position(partt.x);
 //printf("\n %f %f",  total_error(partt.x.f),best_result.f.f[0]);

// Recompute the status
df=total_error(partt.x.f)-error;
if (df<0) partt.status=50;   // Improvement +
if (df==0) partt.status=30;   // No change  =
if (df >0) partt.status=10;     // Deterioration  -

df=error-total_error(part.prev_x.f);    // For the previous move...

if (df<0) partt.status=partt.status+5;   // Improvement   +
if (df==0) partt.status=partt.status+3;   // No change  =
if (df >0) partt.status=partt.status+1;     // Deterioration  -

/* Finally, we have the following table
55  <=> + +
53 <=> = +
51 <=> - +
35 <=> + =
33  <=>  = =
31 <=>  - =
15 <=> + -
13 <=> = -
11 <=> - -
*/

if (problem[level].printlevel>2) display_position(partt.x);
return partt;
}

 //---------------------------------------------------------------------------- PIVOT_CHOICE

struct position pivot_choice(int level)
{
/*
 Choose the pivot according to the memorized positions

 The better a position , the higher the probability it will be chosen as a pivot
Note: memo is a global variable

*/

int                     i;
double               r;

int                     rank_p;
int                     size;
double               sum,sum1;

size=memo[level].size;
rank_p=0;



if (size==1) goto end;

sum=0;
for (i=0;i<size;i++)
{
     sum=sum+1/memo[level].error[i];

}
 r=alea_float(0,sum);
  sum1=1/memo[level].error[0];

 check:
  if (r<=sum1) goto end;

     rank_p=  rank_p+1;
     sum1=sum1+1/memo[level].error[rank_p];

     goto check;

 end:

 //printf("\npivot_choice rank %i    size %i",rank_p,size);
  rank[level]=rank_p; // Global variable. Useful for some strategies
 return memo[level].x[rank_p];
}


// ----------------------------------------------------------------------------- REINIT_SWARM
void reinit_swarm(int level,int option)
{
  /*
  Warning. No special initialisation like TSP or QAP, for the moment
  TO DO later
	*/

	int better;
      int	d;
	int	i;
	int	i_best;
	int	init_option;
	int	j;
	int	j_best;
	struct position best;
// TR is a global variable
if(problem[level].printlevel>1) printf("\n reinit_swarm");
if(problem[level].printlevel>2) display_swarm(1,level);
best=best_result;

best_result.f.f[0]=infinite;

// Modify positions
init_option=1; // will re-init the position in init_particle, but not the links
//init_option=2; // Just recompute the position

for (i=0;i<TR[level].size;i++)
{
	for (j=0;j<TR[level].tr[i].size;j++)
	{

		TR[level].tr[i].part[j]=init_particle(init_option,level,TR[level].tr[i].part[j]);
         better=better_than(TR[level].tr[i].part[j].x.f, best_result.f,level);
       if (better==1)
        {
			i_best=i;
			j_best=j;
			best_result=TR[level].tr[i].part[j].x;
		}
	}
}

if (BIN==1)
	for (d=0;d<best_result.p.size;d++)
		best_result.p.x[d]=floor(best_result.p.x[d]+0.5);


if(option==0) // Eventually keep the best position found so far
{
	if (DYN==1) // The f value for the "best" may have changed
	{
		best.f=MyFunction(best,problem[level].funct,problem[level].target,level);
		eval_f_tot[level]+=-1; // We don't count this evaluation
	}

	better=better_than(best.f, best_result.f,level);

	if (better==1)
	{
		best_result=best;
		TR[level].tr[i_best].part[j_best].x=best;
		TR[level].tr[i_best].part[j_best].p_i=best;
	}
}
	if(problem[level].printlevel>1) display_swarm(1,level);

}

//================================================================  TOTAL_ERROR
 double    total_error(struct f err)
  {
   /*
     Compute the total error of a multiobjective result
   */
    double error;
    int  i;
    if (err.size==1) return err.f[0];


    error=0;
    for (i=0;i<err.size;i++) error=error+err.f[i]*err.f[i];
    error=sqrt(error);
    return error;
  }

