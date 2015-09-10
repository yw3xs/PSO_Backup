#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define ulong unsigned long // For pseudo random numbers generation
#define RAND_MAX_KISS 4294967295.0


#define	Max_C			70   // Max line number for a matrix
#define Max_DD 			30 // Max number of dimensions
                              // Warning. For TSP problems, be sure Max_C=Max_DD>number of nodes
#define	Max_discrete	50 // Max number of discrete possible values for a discrete variable
#define Max_f        10 // Maximum number of objective functions (for multiobjective problems)
#define Max_nodes		50 // Maximum number of nodes in a level (for Neural Network Training)

#define	Max_M			10 // Max vector size for RtoR+ problem
#define  Max_Memo    10 // Max number of memorized positions
#define  Max_run     1000 // Max number of times you can run the same problem
#define Max_status 4
#define Max_swarmsize 	101 // Theoretically infinite, but my computer does not like infinite...
#define nb_strateg   150 // Max number of defined strategies (see move_particle())

struct search_space	{int dim;float min[Max_DD];float max[Max_DD];float dx[Max_DD];float volume;};
struct vector		{int size;double x[Max_DD];};
struct f            {int size;double f[Max_f];};
struct position		{struct vector p;struct f f;};
struct matrix		{int size;double val[Max_C][Max_DD];int nb_cont;double max_cont;};
struct problem 		{int constrain;int nb_f;struct matrix P;int funct;int DD;struct search_space H;float target;float
eps;int printlevel;int Max_Eval; int Max_Eval_2;int Max_Eval_delta;int save;int N;int times;int mine;int all_diff;char data[20];float used_constr;int
init_file;char init_swarm_r[20];char init_swarm_w[20];float max_fr;int init_size;int fuzzy;int peaks;};
struct	i_group		{int size;int label[Max_swarmsize];}; // List of particle labels
// A given particle may belong to several i-groups
struct particle		{int label;struct position x;struct position p_i;struct i_group ig;int status;struct position
prev_x;int mod[Max_DD];double v[Max_DD];};
struct	tribe		{int size;struct particle part[Max_swarmsize];}; // Set of particles
// A given particle belongs to just _one_ tribe
struct	tribe_list	{int size;struct tribe tr[Max_swarmsize];};
struct	vector_c	{int size;double x[Max_C];};
struct	vector_i	{int size; int x[Max_DD];}; // For special integer problem (binary, in particular)
struct memo       {int size;struct position x[Max_Memo];double error[Max_Memo];};
struct	discrete	{int d; int size; double v[Max_discrete];}; // d= rank of the discrete variable
												// size= nb of values, v = value list
struct roots {double z[3]; int status;};// Roots of a polynom a*x^3+...=0

// Subroutines (declarations)
#include "read_display_save.h"
#include "tools.h"
#include "extra_tools.h"
#include "myconstrain.h"
 #include "TSP.h"
 #include "QAP.h"
#include "movpeaks_mc.h"

void              add_memo(struct position x, int level);
double			ANNCOLORCUBE(struct vector weights);
double          ANNParity4(struct vector weights);
double			ANNPIMA(struct vector weights);
double			ANNSERVO(struct vector weights);
double          ANNSINSIMP(struct vector weights);
double          ANNXor(struct vector weights);
double			apple_trees(struct position pos); // For "apple trees" example
struct particle	best_informer(struct particle part,int no_best,int level);
int				better_than(struct f f1,struct f f2,int level);
void clean_run(FILE *f_run,FILE *f_run_clean, int nb_f,int DD,int level);
double coeff_SC(int D);
struct particle complete_part(struct particle par,int option,int level);
struct position	homogen_to_carte(struct position pos) ;
struct particle	init_particle(int option,int level, struct particle part0);
void			link_reorg(int level);
int				local_improv(struct tribe tr,int level);
double 			max_comp(struct position pos);
double          MINLP(struct position pos,int option);
double 			min_comp(struct position pos);

struct particle	move_particle(struct particle part,struct particle partg,int level, int gener);
struct f			MyFunction(struct position pos,int funct,double target,int level);
struct position pivot_choice(int level);
struct position			PSO(int level, float Max_Eval);
ulong	rand_kiss();
void			reinit_swarm(int level,int option);
void	seed_rand_kiss(ulong seed);
double			tot_fifty_fifty(struct position pos);
double               total_error(struct f err);


// Global variables
float			a[Max_M]; // Vector for RtoR+ problem
int				adapt;
double		almostzero=0.000000001;  //to avoid overflow by dividing by too small value;
int                  circular_hood; // Flag for option. Read as data. Useful just to have the "classical
                                             // constricted PSO" for comparison. Usually equal to 0. If not, the value
                                             // is the size of the circular neighbourhood.
int				AS=-3; // Arbitrary value. Means "Assigned"
struct position BEST;
struct position best_result;
int				BIN; // Flag for binary problem
clock_t			clock_tick;
double			cmax;
struct position coeff;
double			coeff_S_C;
int       confin_interv;
struct	discrete discrete[Max_DD];
int	discrete_nb; // Number of special discrete variables
struct particle	dummy_part; // Just as empty parameter
int     DYN; // Flag for dynamic optimisation
double E;
float 			e[Max_DD][Max_M]; // Matrix for RtoR+ problem
double         eval_f[2];
double         eval_f_tot[2];
char			functions[100][100]; // Function names
int				geno_size; // for Moving Peaks. = search space dimension
int				H; // Flag for Coloring problem. Indicates what kind of projection to do
int       HIDDEN; // For Neural Network Training
double			infinite=99999999999999999.;
int       INPUT; // For Neural Network Training
double			khi;
int				label[2]; // To labellize particles
int       landscape[100][Max_DD]; // For binary multimodal problem
int             lexico; // Flag for using or not lexicographical order for multiobjective optim.
int				linkreorg; // Flag for using (1) or not (0) link reorganization
struct memo   memo[2]; // To memorize positions (mainly for pivot method)
 int				NA=-32000;// Arbitrary value. Means "non assigned"
int				NO=-2; // Arbitrary value. Means "no change"

int				max_rand;
//int				Max_swarm; // Just for info
double			Mean_swarm; // Just for info
int               MEMO; // Flag. 1 if positions have to be memorized (pivot methods)
int				MM;
int				nb_pb_max;
double			n_change; // For dynamic optimisation. Number of change
int				no_best;

int       nprogr; // For progressive approach
double			offline_error; // Error for dynamic optimisation
								// modified before each change
double			offline_error_cont; // Error for dynamic optimisation,
									// but continuously computed
int       OUTPUT; // For Neural Network Training
double			phi;
double			pi;
struct problem problem[2];
int   PROGR; // Flag for progressive approach
int            QAP; // Flag for QAP. Use it only for level 0
int               rand_hood; // Flag for using random i-groups or not (>0 =yes, and it gives the size, 0=no)
int               rank[2];
int				recent_change; // Move Peaks.  1 indicates that a change has just ocurred
int				recurs; // Flag for recursive call (particularly for reinit_swarm)
double				retp[Max_nodes];
double   status_count[Max_status]; // Just for information about particle status
float strategies[2][Max_status]; //See  Adaptations and move_particle(). See also explanation in problem.txt
int				times;
struct tribe_list	TR[2];
int               TSP; // Flag. 1 if TSP. Use it only for level 0
double			two_pi;
struct vector	Xmax,Xmin;



//--------
FILE *f_discrete; // Possibles values for discrete variables
				// (if they can't be computed just by giving min, max and granularity)
FILE *f_energy;
FILE *f_init_r; // (optional) to Read an initial swarm
FILE *f_init_w; // to Write an initial swarm for future use
FILE *f_data; // For additional data, defining the problem (see Coloring/Frequency Assignment, RtoR+)
FILE *f_functs; // Function names
FILE *f_matrix;
FILE *f_problem; // Problem file
FILE *f_run; // To save the run
FILE *f_run_c; // To temporarily save the run for multiobjective problem
FILE *f_run_clean; // "Cleaned runs" for a multiobjective problem
FILE *f_vector;
FILE *f_swarm;  // To save swarm
FILE *f_synth; // Summary (mean values) if a given problem is ran several times
FILE *f_trace; // Some additional information (see parameter "save" in problems.txt)
FILE *f_trace_run; // save the number of evaluations and the best result after each iteration
			      //(in order to plot the convergence curve)






