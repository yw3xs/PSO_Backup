/* Moving Peaks Function --- 10/99 */

/* movpeaks.c
 * Copyright (C) 1999 Juergen Branke.

Modified and adapted for Particle Swarm Optimizer TRIBES by Maurice Clerc, 2004-07

 * This is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License.
 *
 * This module generates the Moving Peaks Evaluation Function, a
 * dynamic benchmark problem changing over time.
 */


//======================================================== PARAMETER SETTINGS
/*

T.M.Blackwell and J.Branke, "Multi-swarm optimization in dynamic
environments ", in Applications of Evolutionary Computing, ser LNCS,
G.R.Raidl, Ed., vol 3005. Springer, 2004, pp 489-599
...
The benchmarks parameter settings correspond
to the Scenario 2 as specified on the benchmark web site
http://www.aifb.uni-karlsruhe.de/~jbr/MovPeaks/
The search space has five dimensions [0, 100]^5 with 10 peaks. Peak heights can vary randomly in
the interval [30, 70], width can vary within [1, 12]. Scenario 2 specifies a family of
benchmark function since the initial location, height and width of the peaks, and their
subsequent development is determined by a pseudorandom number generator. Furthermore,
some of the parameters of Scenario 2 are only defined to within a range of
values. Of these, the following choices were made to facilitate comparisons with the
experiments of reference [7]: Change severity vlength = 1.0, correlation lambda =
0.0, and peak change frequency = 5000.
...
*/
int change_frequency = 0; /* number of evaluations between changes. change_frequency
		       =0 means that function never changes (or only if function change_peaks is called)*/


//int geno_size = 5;  // number of dimensions, or the number of double
                      //    valued genes

double vlength = 1; /* distance by which the peaks are moved, severity */

double height_severity=7.0; /* severity of height changes, larger numbers
                               mean larger severity */
double width_severity = 1; /* severity of width changes, larger numbers
                               mean larger severity */

/* lambda determines whether there is a direction of the movement, or whether
   they are totally random. For lambda = 1.0 each move has the same direction,
   while for lambda = 0.0, each move has a random direction */
double lambda=0.0;

int number_of_peaks = 50; /* number of peaks in the landscape */

int use_basis_function=0; /* if set to 1, a static landscape (basis_function) is included in the fitness evaluation */

int calculate_average_error=1; /* saves computation time if not needed and set to 0 */
int calculate_offline_performance = 1; /* saves computation time if not needed and set to 0 */
//int calculate_right_peak = 1; /* saves computation time if not needed and set to 0 */

/* minimum and maximum coordinate in each dimension */
double mincoordinate, maxcoordinate;

/* minimum and maximum height of the peaks
 height chosen randomly when standardheight <= 0.0
 If standardheight>0, the heigth is always equal to standardheight
*/
double minheight = 30.0, maxheight = 70.0, standardheight = 0.0;
/* width chosen randomly when standardwidth = 0.0
 If standardwidth>0, the width is always equal to standardwidth
 */
double minwidth = 1, maxwidth = 12.0, standardwidth = 0.0;


/* Functions */
/* evaluation function */
double eval_movpeaks (double *gen);

/* the following basis functions are provided :*/
double constant_basis_func(double *gen);
double five_peak_basis_func(double *gen);
/* the following peak functions are provided: */
double peak_function1(double *gen, int peak_number);
double peak_function_cone (double *gen, int peak_number);
double peak_function_hilly (double *gen, int peak_number);
double peak_function_twin (double  *gen, int peak_number);


/* allows to set the basis function used */
double (*basis_function) (double *gen)= constant_basis_func;
/* defines the form of a single peak */
double (*peak_function) (double *gen, int peak_number)=peak_function_cone; // You may change it


/****** END OF PARAMETER SECTION *******************/

void change_peaks();   /* preliminary declaration of function change_peaks()*/
//int recent_change; /*  1 indicates that a change has just ocurred */
//int current_peak;      /* peak on which the current best individual is located */
int maximum_peak;       /* number of highest peak */
//double current_maximum; /* fitness value of currently best individual */
//double offline_performance = 0.0;
//double offline_error = 0.0;
double avg_error=0;    /* average error so far */
//double current_error=0;/* error of the currently best individual */
double global_max;     /* absolute maximum in the fitness landscape */
//double global_min // absolute minimun
//int evals = 0;         /* number of evaluations so far */
double **peak;         /* data structure to store peak data */
double *shift;
double *coordinates;
//int *covered_peaks;    /* which peaks are covered by the population ? */

double **prev_movement;/* to store every peak's previous movement */
double dummy_eval (double *gen);

/* Basis Functions */

/* This gives a constant value back to the eval-function that chooses the max of them */
double constant_basis_func(double *gen)
{
  return 0.0;
}

double five_peak_basis_func(double *gen)
{
  int i,j;
  double maximum = -100000.0, dummy;
  static double basis_peak [5] [7] =
  {
    {8.0,  64.0,  67.0,  55.0,   4.0, 0.1, 50.0},
   {50.0,  13.0,  76.0,  15.0,   7.0, 0.1, 50.0},
    {9.0,  19.0,  27.0,  67.0,  24.0, 0.1, 50.0},
   {66.0,  87.0,  65.0,  19.0,  43.0, 0.1, 50.0},
   {76.0,  32.0,  43.0,  54.0,  65.0, 0.1, 50.0},
  };


  for(i=0; i<5; i++)
    {
      dummy = (gen[0]-basis_peak[i][0])*(gen[0]-basis_peak[i][0]);
      for (j=1; j< geno_size; j++)
	dummy += (gen[j]-basis_peak[i][j])*(gen[j]-basis_peak[i][j]);
      dummy = basis_peak[i][geno_size+1]-(basis_peak[i][geno_size]*dummy);
      if (dummy > maximum)
        maximum = dummy;
    }
  return maximum;
}


//======================================================= PEAK FUNCTIONS

/* sharp peaks */
double peak_function1 (double *gen, int peak_number)
{
  int j;
  double dummy;

  //dummy = (gen[0]-peak[peak_number][0])*(gen[0]-peak[peak_number][0]);
  //for (j=1; j< geno_size; j++)
  dummy=0;
	for (j=0; j< geno_size; j++)
    dummy += (gen[j]-peak[peak_number][j])*(gen[j]-peak[peak_number][j]);
  return peak[peak_number][geno_size+1]/(1+(peak[peak_number][geno_size])*dummy);
}

double peak_function_cone (double *gen, int peak_number)
{
  int j;
  double dummy;
  double peak_f;

  //dummy =  (gen[0]-peak[peak_number][0])*(gen[0]-peak[peak_number][0]);
  //for (j=1; j< geno_size; j++)
  //printf("\n");for (j=0;j<geno_size;j++) printf(" %f",gen[j]);

	dummy=0;
	for (j=0; j< geno_size; j++)
    dummy += (gen[j]-peak[peak_number][j])*(gen[j]-peak[peak_number][j]);
  peak_f=peak[peak_number][geno_size+1]-(peak[peak_number][geno_size]*sqrt(dummy));
//printf("\n peak_function_cone %f",peak_f);
  return peak_f ;
}

double peak_function_hilly (double *gen, int peak_number)
{
  int j;
  double dummy;

  dummy =  (gen[0]-peak[peak_number][0])*(gen[0]-peak[peak_number][0]);
  for (j=1; j< geno_size; j++)
    dummy += (gen[j]-peak[peak_number][j])*(gen[j]-peak[peak_number][j]);
  return peak[peak_number][geno_size+1]-(peak[peak_number][geno_size]*dummy)-0.01*sin(20.0*dummy);
}

double peak_function_twin (double  *gen, int peak_number) /* two twin peaks moving together */
{
  int j;
  double maximum = -100000.0, dummy;
  static double twin_peak [7] = /* difference to first peak */
  {
    1.0,  1.0,  1.0,  1.0,   1.0, 0.0, 0.0,
  };

  dummy = pow(gen[0]-peak[peak_number][0],2);
  for (j=1; j< geno_size; j++)
     dummy += pow(gen[j]-peak[peak_number][j],2);
  dummy = peak[peak_number][geno_size+1]-(peak[peak_number][geno_size]*dummy);
  maximum = dummy;
  dummy = pow(gen[j]-(peak[peak_number][0]+twin_peak[0]),2);
  for (j=1; j< geno_size; j++)
     dummy += pow(gen[j]-(peak[peak_number][j]+twin_peak[0]),2);
  dummy = peak[peak_number][geno_size+1]+twin_peak[geno_size+1]-((peak[peak_number][geno_size]+twin_peak[geno_size])*dummy);
  if (dummy > maximum)
    maximum = dummy;

  return maximum;
}

/* The following procedures may be used to change the step size over time */


void change_stepsize_random () /* assigns vlength a value from a normal distribution */
{
  vlength = movnrand();
}

void change_stepsize_linear() /* sinusoidal change of the stepsize, */
{
  static int counter = 1;
  static double frequency = 3.14159/20.0;  /* returns to same value after 20 changes */

  vlength = 1+ sin((double)counter*frequency);
  counter ++;
}



//================================================================= CHANGE_PEAKS
/* whenever this function is called, the peaks are changed */
void change_peaks()
{
  int i,j;
  double sum, sum2, offset, dummy;

  for(i=0; i<number_of_peaks; i++)
    {
      /* shift peak locations */
    sum = 0.0;
    for (j=0; j<geno_size; j++)
	{
      shift[j]=movrand()-0.5;
      sum += shift[j]*shift[j];
     }
    if(sum>0.0)
      sum = vlength/sqrt(sum);
    else                           /* only in case of rounding errors */
      sum = 0.0;
    sum2=0.0;
    for (j=0; j<geno_size; j++)
	{
      shift[j]=sum*(1.0-lambda)*shift[j]+lambda*prev_movement[i][j];
      sum2 += shift[j]*shift[j];
    }
    if(sum2>0.0)
      sum2 = vlength/sqrt(sum2);
    else                           /* only in case of rounding errors */
      sum2 = 0.0;
    for(j=0; j<geno_size; j++)
	{
      shift[j]*=sum2;
      prev_movement[i][j]= shift[j];
      if ((peak[i][j]+prev_movement[i][j]) < mincoordinate)
	  {
        peak[i][j] = 2.0*mincoordinate-peak[i][j]-prev_movement[i][j];
        prev_movement[i][j]*=-1.0;
      }
      else if ((peak[i][j]+prev_movement[i][j]) > maxcoordinate)
	  {
        peak[i][j]= 2.0*maxcoordinate-peak[i][j]-prev_movement[i][j];
        prev_movement[i][j]*=-1.0;
      }
      else
        peak[i][j] += prev_movement[i][j];
     }
    /* change peak width */
    j = geno_size;
    offset = movnrand()*width_severity;
    if ((peak[i][j]+offset) < minwidth)
      peak[i][j] = 2.0*minwidth-peak[i][j]-offset;
    else if ((peak[i][j]+offset) > maxwidth)
      peak[i][j]= 2.0*maxwidth-peak[i][j]-offset;
    else
      peak[i][j] += offset;
   /* change peak height */
    j++;
    offset = height_severity*movnrand();
    if ((peak[i][j]+offset) < minheight)
      peak[i][j] = 2.0*minheight-peak[i][j]-offset;
    else if ((peak[i][j]+offset) > maxheight)
      peak[i][j]= 2.0*maxheight-peak[i][j]-offset;
    else
      peak[i][j] += offset;
//printf("\n peak height %f", peak[i][j]);
  }

  //if(calculate_average_error)
  // Maximum peak value (global_max)
	{
    global_max = -100000.0;
    for (i=0;i<number_of_peaks; i++)
	{
      for (j=0; j<geno_size; j++)
        coordinates[j]=peak[i][j];
      dummy = dummy_eval(coordinates);
      if (dummy>global_max)
	  {
        global_max = dummy;
        maximum_peak = i;
      }
    }
	}

  recent_change = 1;
 //printf("\n eval %.0f Peaks have moved", eval_f_tot[0]);
 // printf("\n max peak %f at ", global_max);
 // for (j=0;j<geno_size;j++) printf(" %f",peak[maximum_peak][j]);

}

//================================================================= DISPLAY_PEAKS
void display_peaks()
{
  int i,j;
  printf("\n Height Width  Position");
	for (i=0; i< number_of_peaks; i++)
	{
		if(i==maximum_peak) printf("\n*"); else printf("\n ");
		printf(" %3.2f  %3.2f   ",peak[i][geno_size+1],peak[i][geno_size]);
		for (j=0; j< geno_size; j++)
		{
			printf(" %3.2f",peak[i][j]);
		}
	}
}

//=========================================================== DUMMY_EVAL
/* dummy evaluation function allows to evaluate without being counted */
double dummy_eval (double *gen)
{
  int i;
  double maximum = -100000.0, dummy;

  for(i=0; i<number_of_peaks; i++)
    {
    dummy = peak_function(gen, i);
    if (dummy > maximum)
      maximum = dummy;
    }

  if (use_basis_function) {

    dummy = basis_function(gen);
    /* If value of basis function is higher return it */
    if (maximum < dummy)
      maximum = dummy;
  }
  return(maximum);
}


//================================================================= EVAL_MOVPEAKS
/* evaluation function */
double eval_movpeaks (double *gen)
{
  int evals;
  double f;

	int i;
  double maximum = -100000.0, dummy;

  evals=eval_f_tot[0]; // for Tribes
  if ((change_frequency > 0)&&(evals%change_frequency==0) && recent_change==0)
  {
   change_peaks();// => recent_change=1
  }

  for(i=0; i<number_of_peaks; i++)
    {
    dummy = peak_function(gen, i);
    if (dummy > maximum)
      maximum = dummy;
    }

  if (use_basis_function)
  {
    dummy = basis_function(gen);
    /* If value of basis function is higher return it */
    if (maximum < dummy)
      maximum = dummy;
  }
 /*
    if (calculate_average_error)
	{
    avg_error+=global_max - maximum;
	}
*/

  /*if (calculate_offline_performance)
  {
    if(recent_change||(maximum > current_maximum))
	{
      current_error = global_max - maximum;

     // if (calculate_right_peak)
	      //current_peak_calc(gen);
     // current_maximum = maximum;
     // recent_change = 0;
    }
    //offline_performance += current_maximum;
    offline_error+= current_error;
  }
*/
  //evals ++;     /* increase the number of evaluations by one */
  f=global_max-maximum; // So that optimisation is minimisation
						// Notice that maximum may be <0

  return(f);
}

//================================================================= INIT_PEAKS
/* initialize all variables at the beginning of the program */

void init_peaks ()
{
  int i,j;
  double dummy;


  shift = (double *) calloc(geno_size, sizeof(double));
  coordinates = (double *) calloc(geno_size, sizeof(double));
//	covered_peaks = (int *) calloc(number_of_peaks, sizeof(int));
  peak = (double **) calloc(number_of_peaks, sizeof(double*));
  prev_movement = (double **) calloc(number_of_peaks, sizeof(double*));
  for (i=0; i< number_of_peaks; i++)
  {
    peak[i]= (double *) calloc(geno_size+2, sizeof(double));
    prev_movement[i] = (double *) calloc(geno_size, sizeof(double));
  }

  for (i=0; i< number_of_peaks; i++)
    for (j=0; j< geno_size; j++)
	{
      //peak[i][j] = 100.0*movrand();
		peak[i][j] =mincoordinate+ (maxcoordinate-mincoordinate)*movrand();
      prev_movement[i][j] = movrand()-0.5;
    }
  if (standardheight <= 0.0)
    for (i=0; i< number_of_peaks; i++)
      peak[i][geno_size+1]= (maxheight-minheight)*movrand()+minheight;
  else
    for (i=0; i< number_of_peaks; i++)
      peak[i][geno_size+1]= standardheight;
  if (standardwidth <= 0.0)
    for (i=0; i< number_of_peaks; i++)
      peak[i][geno_size]= (maxwidth-minwidth)*movrand()+minwidth;
  else
    for (i=0; i< number_of_peaks; i++)
      peak[i][geno_size]= standardwidth;

     // if(calculate_average_error)
	  {
		global_max = -100000.0;
		for (i=0;i<number_of_peaks; i++)
		{
			for (j=0; j<geno_size; j++)
				coordinates[j]=peak[i][j];
			dummy = dummy_eval(coordinates);
			 if (dummy>global_max)
			 { maximum_peak=i;
					global_max = dummy;
			 }
		}
	}

//printf("\n max peak %f , number %i, at \n", global_max,maximum_peak+1);
//for (j=0;j<geno_size;j++) printf(" %f",peak[maximum_peak][j]);
display_peaks();

}



double movrand() { return alea_float(0,1);} // Compatibility with Tribes
double movnrand() { return alea_normal(0,1);}
