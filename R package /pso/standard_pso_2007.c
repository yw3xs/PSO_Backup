/*
Standard PSO 2007
 Contact for remarks, suggestions etc.:
 MauriceClerc@WriteMe.com

 Last update
 2007-12-10 Warning about rotational invariance (valid here only on 2D)
 2007-11-22 stop criterion (option): distance to solution < epsilon
            and log_progress evaluation
 2007-11-21 Ackley function

  -------------------------------- Contributors
 The works and comments of the following persons have been taken
 into account while designing this standard.  Sometimes this is for
 including a feature, and sometimes for leaving out one.

 Auger, Anne
 Blackwell, Tim
 Bratton, Dan
 Clerc, Maurice
 Croussette, Sylvain
 Dattasharma, Abhi
 Eberhart, Russel
 Hansen, Nikolaus
 Keko, Hrvoje
 Kennedy, James
 Krohling, Renato
 Langdon, William
 Liu, Hongbo
 Miranda, Vladimiro
 Poli, Riccardo
 Serra, Pablo
 Stickel, Manfred

 -------------------------------- Motivation
Quite often, researchers claim to compare their version of PSO
with the "standard one", but the "standard one" itself seems to vary!
Thus, it is important to define a real standard that would stay
unchanged for at least one year.
This PSO version does not intend to be the best one on the market
(in particular, there is no adaptation of the swarm size nor of the
coefficients). This is simply very near to the original version (1995),
with just a few improvements based on some recent works.
 --------------------------------- Metaphors
swarm: A team of communicating people (particles)
At each time step
    Each particle chooses a few informants at random, selects the best
    one from this set, and takes into account the information given by
    the chosen particle.
    If it finds no particle better than itself, then the "reasoning" is:
    "I am the best, so I just take my current velocity and my previous
    best position into account"
----------------------------------- Parameters/Options
clamping := true/false => whether to use clamping positions or not
randOrder:= true/false => whether to avoid the bias due to the loop
                on particles "for s = 1 to swarm_size ..." or not
rotation := true/false => whether the algorithm is sensitive
                to a rotation of the landscape or not
You may also modify the following ones, although suggested values
are either hard coded or automatically computed:
S := swarm size
K := maximum number of particles _informed_ by a given one
w := first cognitive/confidence coefficient
c := second cognitive/confidence coefficient
 ----------------------------------- Equations
For each particle and each dimension
Equation 1: v(t+1) = w*v(t) + R(c)*(p(t)-x(t)) + R(c)*(g(t)-x(t))
Equation 2: x(t+1) = x(t) + v(t+1)
where
v(t) := velocity at time t
x(t) := position at time t
p(t) := best previous position of the particle
g(t) := best position amongst the best previous positions
        of the informants of the particle
R(c) := a number coming from a random distribution, which depends on c
In this standard, the distribution is uniform on [0,c]
Note 1:
When the particle has no informant better than itself,
it implies p(t) = g(t)
Therefore, Equation 1 gets modified to:
v(t+1) = w*v(t) + R(c)*(p(t)-x(t))
Note 2:
When the "non sensitivity to rotation" option is activated
(p(t)-x(t)) (and (g(t)-x(t))) are replaced by rotated vectors,
so that the final DNPP (Distribution of the Next Possible Positions)
is not dependent on the system of co-ordinates.
 ----------------------------------- Information links topology
A lot of work has been done about this topic. The main result is this:
There is no "best" topology. Hence the random approach used here.
 ----------------------------------- Initialisation
Initial positions are chosen at random inside the search space
(which is supposed to be a hyperparallelepiped, and often even
a hypercube), according to a uniform distribution.
This is not the best way, but the one used in the original PSO.
Each initial velocity is simply defined as the half-difference of two
random positions. It is simple, and needs no additional parameter.
However, again, it is not the best approach. The resulting distribution
is not even uniform, as is the case for any method that uses a
uniform distribution independently for each component.
The mathematically correct approach needs to use a uniform
distribution inside a hypersphere. It is not very difficult,
and was indeed used in some PSO versions.  However, it is quite
different from the original one.
Moreover, it may be meaningless for some heterogeneous problems,
when each dimension has a different "interpretation".
------------------------------------ From SPSO-06 to SPSO-07
The main differences are:
1. option "non sensitivity to rotation of the landscape"
    Note: although theoretically interesting, this option is quite
        computer time consuming, and the improvement in result may
        only be marginal.
2. option "random permutation of the particles before each iteration"
    Note: same remark. Time consuming, no clear improvement
3. option "clamping position or not"
    Note: in a few rare cases, not clamping positions may induce an
    infinite run, if the stop criterion is the maximum number of
    evaluations

4. probability p of a particular particle being an informant of another
    particle. In SPSO-06 it was implicit (by building the random infonetwork)
    Here, the default value is directly computed as a function of (S,K),
    so that the infonetwork is exactly the same as in SPSO-06.
    However, now it can be "manipulated" ( i.e. any value can be assigned)

5. The search space can be quantised (however this algorithm is _not_
   for combinatorial problems)
Also, the code is far more modular. It means it is slower, but easier
to translate into another language, and easier to modify.
 ----------------------------------- Use
 Define the problem (you may add your own one in problemDef() and perf())
 Choose your options
 Run and enjoy!

 ------------------------------------ Some results
 So that you can check your own implementation, here are some results.

(clamping, randOrder, rotation, stop) = (1,1,1,0)
Function    Domain          error   Nb_of_eval  Nb_of_runs  Result
Parabola    [-100,100]^30   0.0001  6000        100          52% (success rate)
    shifted     ""          ""      ""          ""           7%
    shifted     ""          ""      7000        ""           100%
Griewank    [-100,100]^30   0.0001  9000        100          51%
Rosenbrock  [-10,10]^30     0.0001  40000       50           31
Rastrigin   [-10,10]^30     0.0001  40000       50           53 (error mean)
Tripod      [-100,100]^2    0.0001  10000       50           50%
(clamping, randOrder, rotation, stop) = (1,1,0,0)
Parabola    [-100,100]^30   0.0001  6000        100          0.69 (0%)
    shifted     ""          ""      ""          ""           2.16 (0%)
Griewank    [-100,100]^30   0.0001  9000        100          14%
Rosenbrock  [-10,10]^30     0.0001  40000       100          38.7
Rastrigin   [-10,10]^30     0.0001  40000       100          54.8 (error mean)
Tripod      [-100,100]^2    0.0001  10000       100          47%
Because of randomness, and also depending on your own pseudo-random
number generator, you just have to find similar values,
not exactly these ones.
 */

#include "stdio.h"
#include "math.h"
#include <stdlib.h>
#include <time.h>

#define D_max 100       // Max number of dimensions of the search space
#define R_max 200       // Max number of runs
#define S_max 100       // Max swarm size
#define zero  0         // 1.0e-30 // To avoid numerical instabilities
  // Structures
struct quantum
{
  double q[D_max];
  int size;
};

struct SS
{
    int D;
    double max[D_max];
    double min[D_max];
    struct quantum q;       // Quantisation step size. 0 => continuous problem
};

struct param
{
    double c;       // Confidence coefficient
    int clamping;   // Position clamping or not
    int K;          // Max number of particles informed by a given one
    double p;       // Probability threshold for random topology
                    // (is actually computed as p(S,K) )
    int randOrder;  // Random choice of particles or not
    int rotation;   // Sensitive to rotation or not
    int S;          // Swarm size
    int stop;       // Flag for stop criterion
    double w;       // Confidence coefficient
};

struct position
{
    double f;
    int improved;
    int size;
    double x[D_max];
};

struct velocity
{
    int size;
    double v[D_max];
};

struct problem
{
    double epsilon;     // Admissible error
    int evalMax;        // Maximum number of fitness evaluations
    int function;       // Function code
    double objective;   // Objective value
                        // Solution position (if known, just for tests)
    struct position solution;
    struct SS SS;       // Search space
};

struct swarm
{
    int best;                   // rank of the best particle
    struct position P[S_max];   // Previous best positions found by each particle
    int S;                      // Swarm size
    struct velocity V[S_max];   // Velocities
    struct position X[S_max];   // Positions
};

struct result
{
    double nEval;       // Number of evaluations
    struct swarm SW;    // Final swarm
    double error;       // Numerical result of the run
};
struct matrix   // Useful for "non rotation sensitive" option
{
    int size;
    double v[D_max][D_max];
};

  // Sub-programs
double alea (double a, double b);
int alea_integer (int a, int b);
double alea_normal (double mean, double stdev);
double distanceL(struct position x1, struct position x2, double L);
struct velocity aleaVector(int D,double coeff);
struct matrix matrixProduct(struct matrix M1,struct matrix M2);
struct matrix matrixRotation(struct velocity V);
struct velocity matrixVectProduct(struct matrix M,struct velocity V);
double normL (struct velocity v,double L);
double perf (struct position x, int function,struct SS SS); // Fitness evaluation
struct position quantis (struct position x, struct SS SS);
struct problem problemDef(int functionCode);
struct result PSO ( struct param param, struct problem problem);
int sign (double x);
  // Global variables
long double E;          // exp(1). Useful for some test functions
long double pi;         // Useful for some test functions
struct matrix reflex1;
long double sqrtD;

  // File(s);
FILE * f_run;
FILE * f_synth;

  // =================================================
int main ()
{

    int d;          // Current dimension
    double error;           // Current error
    double errorMean;       // Average error
    double errorMin;        // Best result over all runs
    double errorMeanBest[R_max];
    double evalMean;        // Mean number of evaluations
    int functionCode;
    int i,j;
    int nFailure;       // Number of unsuccessful runs
    double logProgressMean;
    struct param param;
    struct problem pb;
    int run, runMax;
    struct result result;
    double successRate;
    double variance;

    f_run = fopen ("f_run.txt", "w");
    f_synth = fopen ("f_synth.txt", "w");

    E = exp ((long double) 1);
    pi = acos ((long double) -1);


// ----------------------------------------------- PROBLEM
    functionCode = 0;
            /* (see problemDef( ) for precise definitions)
                 0 Parabola (Sphere)
                 1 Griewank
                 2 Rosenbrock (Banana)
                 3 Rastrigin
                 4 Tripod (dimension 2)
                 5 Ackley
                100 Shifted Parabola
                99 Test
             */

    runMax = 30; // Numbers of runs
    if (runMax > R_max) runMax = R_max;


    // -----------------------------------------------------
    // PARAMETERS
    // * means "suggested value"

    param.clamping =1;
            // 0 => no clamping AND no evaluation. WARNING: the program
            //              may NEVER stop (in particular with option move 20 (jumps)) *
            // *1 => classical. Set to bounds, and velocity to zero

    param.randOrder=1; // 0 => at each iteration, particles are modified
                        //     always according to the same order 0..S-1
                        //*1 => at each iteration, particles numbers are
                        //      randomly permutated
    param.rotation=0; // WARNING. Quite time consuming!
    // WARNING. Experimental code, completely valid only for dimension 2
    // 0 =>  sensitive to rotation of the system of coordinates
    //*1 => non sensitive (except side effects)
    // by using a rotated hypercube for the probability distribution

    param.stop = 0; // Stop criterion
                    // 0 => error < pb.epsilon
                    // 1 => eval < pb.evalMax
                    // 2 => ||x-solution|| < pb.epsilon


    // -------------------------------------------------------
    // Some information
    printf ("\n Function %i ", functionCode);
    printf("\n (clamping, randOrder, rotation, stop_criterion) = (%i, %i, %i, %i)",
        param.clamping, param.randOrder, param.rotation, param.stop);

    // ===========================================================
    // RUNs

    // Initialize some objects
    pb=problemDef(functionCode);

    // You may "manipulate" S, p, w and c
    // but here are the suggested values
    param.S = (int) (10 + 2 * sqrt(pb.SS.D));   // Swarm size
    if (param.S > S_max) param.S = S_max;

    printf("\n Swarm size %i", param.S);

    param.K=3;
    param.p=1-pow(1-(double)1/(param.S),param.K);

    // According to Clerc's Stagnation Analysis
    param.w = 1 / (2 * log ((double) 2)); // 0.721
    param.c = 0.5 + log ((double) 2); // 1.193

    // According to Poli's Sampling Distribution of PSOs analysis
    //param.w = ??; // in [0,1[
    //param.c =
    //    smaller than 12*(param.w*param.w-1)/(5*param.w -7);


    printf("\n c = %f,  w = %f",param.c, param.w);

    //---------------
    sqrtD=sqrt((long double) pb.SS.D);

    // Define just once the first reflexion matrix
    if(param.rotation>0)
    {

        reflex1.size=pb.SS.D;
        for (i=0;i<pb.SS.D;i++)
        {

            for (j=0;j<pb.SS.D;j++)
            {

                reflex1.v[i][j]=-2.0/pb.SS.D;

            }

        }

        for (d=0;d<pb.SS.D;d++)
        {

            reflex1.v[d][d]=1+reflex1.v[d][d];

        }

    }

    errorMean = 0;
    evalMean = 0;
    nFailure = 0;
    //------------------------------------- RUNS
    for (run = 0; run < runMax; run++)
    {

        //srand (clock () / 100);   // May improve pseudo-randomness
        result = PSO (param, pb);
        error = result.error;

        if (error > pb.epsilon) // Failure
        {

            nFailure = nFailure + 1;

}


            // Result display
        printf ("\nRun %i. Eval %f. Error %f ", run, result.nEval, error);
        printf("  x(0)= %f",result.SW.P[result.SW.best].x[0]);
            // Save result
            // fprintf( f_run, "\n%i %.0f %f ", run, result.nEval,  error );
            // for ( d = 0; d < SS.D; d++ ) fprintf( f_run, " %f",  bestResult.x[d] );

            // Compute/store some statistical information
        if (run == 0)
            errorMin = error;
        else if (error < errorMin)
            errorMin = error;
        evalMean = evalMean + result.nEval;
        errorMean = errorMean + error;
        errorMeanBest[run] = error;
        logProgressMean  = logProgressMean - log(error);

}       // End loop on "run"

    // ---------------------END
    // Display some statistical information

    evalMean = evalMean / (double) runMax;
    errorMean = errorMean / (double) runMax;
    logProgressMean = logProgressMean/(double) runMax;

    printf ("\n Eval. (mean)= %f", evalMean);
    printf ("\n Error (mean) = %f", errorMean);

    // Variance
    variance = 0;

    for (run = 0; run < runMax; run++)
                variance = variance + pow (errorMeanBest[run] - errorMean, 2);

    variance = sqrt (variance / runMax);
    printf ("\n Std. dev. %f", variance);
    printf("\n Log_progress (mean) = %f", logProgressMean);
    // Success rate and minimum value
    printf("\n Failure(s) %i",nFailure);
    successRate = 100 * (1 - nFailure / (double) runMax);
    printf ("\n Success rate = %.2f%%", successRate);

    if (run > 1)
        printf ("\n Best min value = %f", errorMin);

    // Save
    /*
        fprintf(f_synth,"\n"); for (d=0;d<SS.D;d++) fprintf(f_synth,"%f ",
            pb.offset[d]);
        fprintf(f_synth,"    %f %f %f %.0f%% %f",errorMean,variance,errorMin,
            successRate,evalMean);

     fprintf( f_synth, "\n%f %f %f %f %.0f%% %f ", shift,
            errorMean, variance, errorMin, successRate, evalMean );
    */
    fprintf (f_synth, "\n");
    fprintf (f_synth, "%f %f %.0f%% %f   ",
                    errorMean, variance, successRate, evalMean);

    return 0; // End of main program
}
// ===============================================================
// PSO
struct result PSO (struct param param, struct problem pb)
{

    struct velocity aleaV;
    int d;
    double error;
    double errorPrev;
    struct velocity expt1,expt2;
    int g;
    struct velocity GX;
    int index[S_max], indexTemp[S_max];
    int initLinks;  // Flag to (re)init or not the information links
    int iter;       // Iteration number (time step)
    int iterBegin;
    int length;
    int LINKS[S_max][S_max];    // Information links
    int m;
    int noEval;
    double normPX, normGX;
    int noStop;
    int outside;
    double p;
    struct velocity PX;
    struct result R;
    int rank;
    struct matrix RotatePX;
    struct matrix RotateGX;
    int s0, s,s1;
    int t;
    double zz;

    aleaV.size=pb.SS.D;
    RotatePX.size=pb.SS.D;
    RotateGX.size=pb.SS.D;

    // -----------------------------------------------------
    // INITIALISATION
    p=param.p; // Probability threshold for random topology
    R.SW.S = param.S; // Size of the current swarm

    // Position and velocity
    for (s = 0; s < R.SW.S; s++)
    {

        R.SW.X[s].size = pb.SS.D;
        R.SW.V[s].size = pb.SS.D;

        for (d = 0; d < pb.SS.D; d++)
        {

            R.SW.X[s].x[d] = alea (pb.SS.min[d], pb.SS.max[d]);

}

        for (d = 0; d < pb.SS.D; d++)
        {

            R.SW.V[s].v[d] =
            (alea( pb.SS.min[d], pb.SS.max[d] ) - R.SW.X[s].x[d])/2;

}

}

    // Take quantisation into account
    R.SW.X[s] = quantis (R.SW.X[s], pb.SS);

    // First evaluations
    for (s = 0; s < R.SW.S; s++)
    {

        R.SW.X[s].f =
            perf (R.SW.X[s], pb.function,pb.SS);

        R.SW.P[s] = R.SW.X[s];  // Best position = current one
        R.SW.P[s].improved = 0; // No improvement

}

    // If the number max of evaluations is smaller than
    // the swarm size, just keep evalMax particles, and finish
    if (R.SW.S>pb.evalMax) R.SW.S=pb.evalMax;
    R.nEval = R.SW.S;

    // Find the best
    R.SW.best = 0;
    switch (param.stop)
    {

        default:
        errorPrev = fabs(pb.epsilon-R.SW.P[R.SW.best].f);
        break;

        case 2:
        errorPrev=distanceL(R.SW.P[R.SW.best],pb.solution,2);
        break;

}

    for (s = 1; s < R.SW.S; s++)
    {

        switch (param.stop)
        {

            default:
            zz=fabs(pb.epsilon-R.SW.P[s].f);
            if (zz < errorPrev)
                R.SW.best = s;
                errorPrev=zz;
            break;

            case 2:
            zz=distanceL(R.SW.P[R.SW.best],pb.solution,2);
            if (zz<errorPrev)
                R.SW.best = s;
                errorPrev=zz;
            break;

}

}



/*
    // Display the best
    printf( " Best value after init. %f ", errorPrev );
    printf( "\n Position :\n" );
    for ( d = 0; d < SS.D; d++ ) printf( " %f", R.SW.P[R.SW.best].x[d] );
*/
    initLinks = 1;      // So that information links will beinitialized
            // Note: It is also a flag saying "No improvement"
    noStop = 0;

    // ---------------------------------------------- ITERATIONS
    iter=0; iterBegin=0;
    while (noStop == 0)
    {

        iter=iter+1;

        if (initLinks==1)   // Random topology
        {

            // Who informs who, at random
            for (s = 0; s < R.SW.S; s++)
            {

                for (m = 0; m < R.SW.S; m++)
                {

                    if (alea (0, 1)<p) LINKS[m][s] = 1;  // Probabilistic method
                    else LINKS[m][s] = 0;

}

}
            // Each particle informs itself
            for (m = 0; m < R.SW.S; m++)
            {

                    LINKS[m][m] = 1;

}

}

        // The swarm MOVES
        //printf("\nIteration %i",iter);
        for (s = 0; s < R.SW.S; s++)  index[s]=s;

        switch (param.randOrder)
        {

            default:
            break;

            case 1: //Define a random permutation
            length=R.SW.S;
            for (s=0;s<length;s++) indexTemp[s]=index[s];

            for (s=0;s<R.SW.S;s++)
            {

                rank=alea_integer(0,length-1);
                index[s]=indexTemp[rank];
                if (rank<length-1)   // Compact
                {

                    for (t=rank;t<length;t++)
                        indexTemp[t]=indexTemp[t+1];

}
                        length=length-1;

}
            break;

}
        for (s0 = 0; s0 < R.SW.S; s0++)  // For each particle ...
        {

            s=index[s0];
            // ... find the first informant
            s1 = 0;
            while (LINKS[s1][s] == 0)   s1++;
            if (s1 >= R.SW.S)    s1 = s;

        // Find the best informant
        g = s1;
        for (m = s1; m < R.SW.S; m++)
        {

            if (LINKS[m][s] == 1 && R.SW.P[m].f < R.SW.P[g].f)
                    g = m;

}

        //.. compute the new velocity, and move

        // Exploration tendency
        for (d = 0; d < pb.SS.D; d++)
        {

            R.SW.V[s].v[d]=param.w *R.SW.V[s].v[d];

}

        // Prepare Exploitation tendency
        for (d = 0; d < pb.SS.D; d++)
        {

            PX.v[d]= R.SW.P[s].x[d] - R.SW.X[s].x[d];

}
        PX.size=pb.SS.D;

        if(g!=s)
        {

            for (d = 0; d < pb.SS.D; d++)
            {

                GX.v[d]= R.SW.P[g].x[d] - R.SW.X[s].x[d];

}
            GX.size=pb.SS.D;

}

        // Option "non sentivity to rotation"
        if (param.rotation>0)
        {

            normPX=normL(PX,2);
            if (g!=s) normGX=normL(GX,2);
            if(normPX>0)
            {

                RotatePX=matrixRotation(PX);

}

            if(g!= s && normGX>0)
            {

                RotateGX=matrixRotation(GX);

}

}

        // Exploitation tendencies
        switch (param.rotation)
        {

            default:
            for (d = 0; d < pb.SS.D; d++)
            {

                R.SW.V[s].v[d]=R.SW.V[s].v[d] +
                +   alea(0, param.c)*PX.v[d];

}

            if (g!=s)
            {

                for (d = 0; d < pb.SS.D; d++)
                {

                    R.SW.V[s].v[d]=R.SW.V[s].v[d]
                    +   alea(0,param.c) * GX.v[d];

}

}
            break;

            case 1:
            // First exploitation tendency
            if(normPX>0)
            {

                zz=param.c*normPX/sqrtD;
                aleaV=aleaVector(pb.SS.D, zz);
                expt1=matrixVectProduct(RotatePX,aleaV);

                for (d = 0; d < pb.SS.D; d++)
                {

                    R.SW.V[s].v[d]=R.SW.V[s].v[d]+expt1.v[d];

}

}

            // Second exploitation tendency
            if(g!=s && normGX>0)
            {

                zz=param.c*normGX/sqrtD;
                aleaV=aleaVector(pb.SS.D, zz);
                expt2=matrixVectProduct(RotateGX,aleaV);

                for (d = 0; d < pb.SS.D; d++)
                {

                    R.SW.V[s].v[d]=R.SW.V[s].v[d]+expt2.v[d];

}

}
            break;

}

        // Update the position
        for (d = 0; d < pb.SS.D; d++)
        {

            R.SW.X[s].x[d] = R.SW.X[s].x[d] + R.SW.V[s].v[d];

}

        if (R.nEval >= pb.evalMax) goto end;

        // --------------------------
        noEval = 1;

        // Quantisation
        R.SW.X[s] = quantis (R.SW.X[s], pb.SS);

        switch (param.clamping)
        {

            case 0: // No clamping AND no evaluation
            outside = 0;

            for (d = 0; d < pb.SS.D; d++)
            {

                if (R.SW.X[s].x[d] < pb.SS.min[d] || R.SW.X[s].x[d] > pb.SS.max[d])
                                    outside++;

}

            if (outside == 0)   // If inside, the position is evaluated
            {

                R.SW.X[s].f =
                    perf (R.SW.X[s], pb.function, pb.SS);
                R.nEval = R.nEval + 1;

}
            break;

            case 1: // Set to the bounds, and v to zero
            for (d = 0; d < pb.SS.D; d++)
            {

                if (R.SW.X[s].x[d] < pb.SS.min[d])
                {

                    R.SW.X[s].x[d] = pb.SS.min[d];
                    R.SW.V[s].v[d] = 0;

}

                if (R.SW.X[s].x[d] > pb.SS.max[d])
                {

                    R.SW.X[s].x[d] = pb.SS.max[d];
                    R.SW.V[s].v[d] = 0;

}

}

            R.SW.X[s].f =perf(R.SW.X[s],pb.function, pb.SS);
            R.nEval = R.nEval + 1;
            break;

}

            // ... update the best previous position
            if (R.SW.X[s].f < R.SW.P[s].f)   // Improvement
            {

                R.SW.P[s] = R.SW.X[s];

                // ... update the best of the bests
                if (R.SW.P[s].f < R.SW.P[R.SW.best].f)
                {

                    R.SW.best = s;

}

}

}           // End of "for (s=0 ...  "

        // Check if finished
        switch (param.stop)
        {

            default:
            error = R.SW.P[R.SW.best].f;
            break;

            case 2:
            error=distanceL(R.SW.P[R.SW.best],pb.solution,2);
            break;

}
        error= fabs(error - pb.epsilon);

        if (error < errorPrev)   // Improvement
        {

            initLinks = 0;

}
        else            // No improvement
        {

            initLinks = 1;  // Information links will be    reinitialized

}

        errorPrev = error;
    end:
        switch (param.stop)
        {

            case 0:
            case 2:
            if (error > pb.epsilon && R.nEval < pb.evalMax)
                    noStop = 0; // Won't stop
            else
                noStop = 1; // Will stop
            break;

            case 1:
            if (R.nEval < pb.evalMax)
                    noStop = 0; // Won't stop
            else
                    noStop = 1; // Will stop
            break;

}


} // End of "while nostop ...

    // printf( "\n and the winner is ... %i", R.SW.best );
    // fprintf( f_stag, "\nEND" );
    R.error = error;
    return R;
}


    // ===========================================================
double alea (double a, double b)
{
                // random number (uniform distribution) in  [a b]
  double r;
  r = (double) rand ();
  r = r / RAND_MAX;
  return a + r * (b - a);
}

// ===========================================================
int alea_integer (int a, int b)
{
                // Integer random number in [a b]
  int ir;
  double r;

    r = alea (0, 1);
  ir = (int) (a + r * (b + 1 - a));

    if (ir > b)  ir = b;

    return ir;
}

    // ===========================================================
double alea_normal (double mean, double std_dev)
{

  /*
  Use the polar form of the Box-Muller transformation to obtain a pseudo
    random number from a Gaussian distribution
  */
  double x1, x2, w, y1;
      // double y2;

  do
  {

        x1 = 2.0 * alea (0, 1) - 1.0;
        x2 = 2.0 * alea (0, 1) - 1.0;
        w = x1 * x1 + x2 * x2;

}
  while (w >= 1.0);

    w = sqrt (-2.0 * log (w) / w);
  y1 = x1 * w;
  // y2 = x2 * w;
  if(alea(0,1)<0.5) y1=-y1;
  y1 = y1 * std_dev + mean;
  return y1;
}


// =============================================================
struct velocity aleaVector(int D,double coeff)
{

    struct velocity V;
    int d;
    int i;
    int K=2; // 1 => uniform distribution in a hypercube
            // 2 => "triangle" distribution
    double rnd;

    V.size=D;

    for (d=0;d<D;d++)
    {

        rnd=0;
        for (i=1;i<=K;i++) rnd=rnd+alea(0,1);
        V.v[d]=rnd*coeff/K;

}

    return V;
}

      // ===========================================================
double normL (struct velocity v,double L)
{
   // L-norm of a vector
    int d;
    double n;

    n = 0;

    for (d = 0; d < v.size; d++)
        n = n + pow(fabs(v.v[d]),L);

    n = pow (n, 1/L);
    return n;
}

// ===========================================================
double distanceL (struct position x1, struct position x2,double L)
{
  // Distance between two positions
    // L = 2 => Euclidean
    int d;
    double n;

    n = 0;

    for (d = 0; d < x1.size; d++)
        n = n + pow (fabs(x1.x[d] - x2.x[d]), L);

    n = pow (n, 1/L);
  return n;
}

//============================================================
struct matrix matrixProduct(struct matrix M1,struct matrix M2)
{

    // Two square matrices of same size
    struct matrix Product;
    int D;
    int i,j,k;
    double sum;
    D=M1.size;
    for (i=0;i<D;i++)
    {

        for (j=0;j<D;j++)
        {

            sum=0;
            for (k=0;k<D;k++)
            {

                sum=sum+M1.v[i][k]*M2.v[k][j];

}
            Product.v[i][j]=sum;

}

}
    Product.size=D;
    return Product;
}
//=========================================================
struct matrix matrixRotation(struct velocity V)
{

    /*
        Define the matrice of the rotation V' => V
    where V'=(1,1,...1)*normV/sqrt(D)  (i.e. norm(V') = norm(V) )

    */
    struct velocity B;
    int i,j,d, D;
    double normB,normV,normV2;
    //struct matrix reflex1; // Global variable
    struct matrix reflex2;
    struct matrix rotateV;
    double temp;

    D=V.size;
    normV=normL(V,2); normV2=normV*normV;
    reflex2.size=D;

    // Reflection relatively to the vector V'=(1,1, ...1)/sqrt(D)
    // norm(V')=1
    // Has been computed just once  (global matrix reflex1)

    //Define the "bisectrix" B of (V',V) as an unit vector
    B.size=D;
    temp=normV/sqrtD;

    for (d=0;d<D;d++)
    {

        B.v[d]=V.v[d]+temp;

}
    normB=normL(B,2);

    if(normB>0)
    {

        for (d=0;d<D;d++)
        {

            B.v[d]=B.v[d]/normB;

}

}

    // Reflection relatively to B
    for (i=0;i<D;i++)
    {

        for (j=0;j<D;j++)
        {

            reflex2.v[i][j]=-2*B.v[i]*B.v[j];

}

}

    for (d=0;d<D;d++)
    {

        reflex2.v[d][d]=1+reflex2.v[d][d];

}

    // Multiply the two reflections
    // => rotation
    rotateV=matrixProduct(reflex2,reflex1);
    return rotateV;

}
//==========================================================
struct velocity matrixVectProduct(struct matrix M,struct velocity V)
{

    struct velocity Vp;
    int d,j;
    int Dim;
    double sum;
    Dim=V.size;
    for (d=0;d<Dim;d++)
    {

        sum=0;
        for (j=0;j<Dim;j++)
        {

            sum=sum+M.v[d][j]*V.v[j];

}
        Vp.v[d]=sum;

}
    Vp.size=Dim;
    return Vp;
}
// ===========================================================
int sign (double x)
{

    if (x == 0) return 0;
    if (x < 0)   return -1;
            return 1;
}

// ===========================================================
struct position quantis (struct position x, struct SS SS)
{

    /*
      Quantisatition of a position
    Only values like x+k*q (k integer) are admissible
     */
  int d;
  double qd;
  struct position quantx;

    quantx = x;
    for (d = 0; d < x.size; d++)
    {

        qd = SS.q.q[d];

        if (qd > zero)   // Note that qd can't be < 0
      {

            qd = qd * (SS.max[d] - SS.min[d]) / 2;
            quantx.x[d] = qd * floor (0.5 + x.x[d] / qd);

}

}
    return quantx;
}
// ===========================================================
double perf (struct position x, int function, struct SS SS)
{
                // Evaluate the fitness value for the particle of rank s
    int d;
    double DD;
    int  k;
    double f, p, xd, x1, x2;
    double s11, s12, s21, s22;
    double sum1,sum2;
    double t0, tt, t1;
    struct position xs;
    // Shifted Parabola/Sphere (CEC 2005 benchmark)
    static double offset_0[30] =
{

-3.9311900e+001, 5.8899900e+001, -4.6322400e+001, -7.4651500e+001, -1.6799700e+001,
-8.0544100e+001, -1.0593500e+001, 2.4969400e+001, 8.9838400e+001, 9.1119000e+000,
-1.0744300e+001, -2.7855800e+001, -1.2580600e+001, 7.5930000e+000, 7.4812700e+001,
 6.8495900e+001, -5.3429300e+001, 7.8854400e+001, -6.8595700e+001, 6.3743200e+001,
 3.1347000e+001, -3.7501600e+001, 3.3892900e+001, -8.8804500e+001, -7.8771900e+001,
-6.6494400e+001, 4.4197200e+001, 1.8383600e+001, 2.6521200e+001, 8.4472300e+001
};
    xs = x;

    switch (function)
    {


        case 100:
        for (d = 0; d < xs.size; d++)
        {

            xs.x[d]=xs.x[d]-offset_0[d];

}

        case 0:     // Parabola (Sphere)
        f = 0;

        for (d = 0; d < xs.size; d++)
        {

            xd = xs.x[d];
            f = f + xd * xd;

}
        break;

        case 1:     // Griewank
      f = 0;
        p = 1;

        for (d = 0; d < xs.size; d++)
      {

            xd = xs.x[d];
            f = f + xd * xd;
            p = p * cos (xd / sqrt ((double) (d + 1)));

}
        f = f / 4000 - p + 1;
        break;

        case 2:     // Rosenbrock

      f = 0;
        t0 = xs.x[0]  + 1;  // Solution on (0,...0) when
                                                                // offset=0
      for (d = 1; d < xs.size; d++)
      {


            t1 = xs.x[d]  + 1;
            tt = 1 - t0;
            f += tt * tt;
            tt = t1 - t0 * t0;
            f += 100 * tt * tt;
            t0 = t1;

}
        break;

        case 3:     // Rastrigin
      k = 10;
        f = 0;

        for (d = 0; d < xs.size; d++)
      {

            xd = xs.x[d];
            f =f+ xd * xd - k * cos (2 * pi * xd);

}
        f =f+ xs.size * k;
        break;

        case 4:     // 2D Tripod function
      // Note that there is a big discontinuity right on the solution
      // point.

    x1 = xs.x[0] ;
    x2 = xs.x[1];
    s11 = (1.0 - sign (x1)) / 2;
    s12 = (1.0 + sign (x1)) / 2;
    s21 = (1.0 - sign (x2)) / 2;
    s22 = (1.0 + sign (x2)) / 2;

    //f = s21 * (fabs (x1) - x2); // Solution on (0,0)
    f = s21 * (fabs (x1) +fabs(x2+50)); // Solution on (0,-50)
    f = f + s22 * (s11 * (1 + fabs (x1 + 50) +
            fabs (x2 - 50)) + s12 * (2 +
            fabs (x1 - 50) +
            fabs (x2 - 50)));
    break;
    case 5:  // Ackley
    sum1=0;
    sum2=0;
    DD=x.size;
    pi=acos(-1);
    for (d=0;d<x.size;d++)
        {

        xd=xs.x[d];
        sum1=sum1+xd*xd;
        sum2=sum2+cos(2*pi*xd);

}
    f=-20*exp(-0.2*sqrt(  sum1/DD  ))-exp(sum2/DD)+20+exp(1);

    break;

    case 99: // Test
    xd=xs.x[0];
    //if (xd<9) f=10-xd; else f=10*xd-90;
    if (xd<=1) f=10*xd; else f=11-xd;
    break;


}

return f;
}

//===================================================
struct problem problemDef(int functionCode)
{

    int d;
    struct problem pb;

    pb.function=functionCode;
    pb.epsilon = 0.0001;    //0.0001 Acceptable error
    pb.objective = 0;       // Objective value

    // Define the solution point, for test
    // NEEDED when param.stop = 2
    // i.e. when stop criterion is distance_to_solution < epsilon
    for (d=0; d<30;d++)
    {

        pb.solution.x[d]=0;

}


    // ------------------ Search space
    switch (pb.function)
    {

        case 0:         // Parabola
        case 100:
        pb.SS.D =30;// 30   // Dimension                                                                                                // values

        for (d = 0; d < pb.SS.D; d++)
        {

            pb.SS.min[d] = -100; // -100
            pb.SS.max[d] = 100; // 100
            pb.SS.q.q[d] = 0;   // Relative quantisation, in [0,1].

}

        pb.evalMax = 6000;// 6000   // Max number of evaluations for each run
        break;

        case 1:     // Griewank
        pb.SS.D = 30;   // Dimension

         // Boundaries
         for (d = 0; d < pb.SS.D; d++)
         {

                pb.SS.min[d] = -100;
                pb.SS.max[d] = 100;
                pb.SS.q.q[d] = 0;

}

        pb.evalMax = 9000;
        break;

        case 2:     // Rosenbrock
        pb.SS.D = 30;   // Dimension

                    // Boundaries
        for (d = 0; d < pb.SS.D; d++)
        {

            pb.SS.min[d] = -10; pb.SS.max[d] = 10;
            pb.SS.q.q[d] = 0;

}

        pb.evalMax =40000; // 40000
        break;


        case 3:     // Rastrigin
        pb.SS.D = 30;   // Dimension

                    // Boundaries
        for (d = 0; d < pb.SS.D; d++)
        {

            pb.SS.min[d] =-10;
            pb.SS.max[d] =10;
            pb.SS.q.q[d] = 0;

}

        pb.evalMax = 40000;
        break;

        case 4:     // Tripod
        pb.SS.D = 2;    // Dimension

                    // Boundaries
        for (d = 0; d < pb.SS.D; d++)
        {

            pb.SS.min[d] = -100;
            pb.SS.max[d] = 100;
            pb.SS.q.q[d] = 0;

}

        pb.evalMax = 10000;
        break;
        case 5: // Ackley
        pb.SS.D = 10;   // Dimension
        // Boundaries
        for (d = 0; d < pb.SS.D; d++)
        {

            pb.SS.min[d] = -100; // -32
            pb.SS.max[d] = 100; // 32
            pb.SS.q.q[d] = 0;

}
        pb.evalMax = 6000;
        break;


    // case n: // Add here your own problem

    case 99: // Test
    pb.SS.D = 1;    // Dimension
        // Boundaries
        for (d = 0; d < pb.SS.D; d++)
        {

            pb.SS.min[d] = 0;
            pb.SS.max[d] = 10;
            pb.SS.q.q[d] = 0;

}
        pb.evalMax = 10000;
        break;


}

    pb.SS.q.size = pb.SS.D;
    return pb;
}
