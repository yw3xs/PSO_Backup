//======================================================== CHECK_SEQ
 int check_seq(struct position s)
 // Check how many times there is two times the same value
 {
  int check;
  int i,j;

  check=0;

  for (i=0;i<s.p.size-1;i++)
 {
     for (j=i+1;j<s.p.size;j++)
     {
        if ((int)s.p.x[i]==(int)s.p.x[j]) check=check+1;
     }
 }
 return check;
 }

 //======================================================== F_TSP
 double f_tsp(struct position pos,int option, int penalty,int rank1,int rank2,int level)
 {
 /*
   The graph matrix P   is a global variable, in problem
   problem is a global variable
   If rank1>=0, positions in rank1 and rank2 have just been switched
 */
 int        d;
 double  delta_eval;
 int        DD;
 int        i,j;
 double total;

   DD=pos.p.size;
 /*
   When there has been just a switch between positions, the f_value could be (and will be in next version)
   evaluated by updating  at most 4 arcs values
 */
 if (rank1<0) delta_eval=1; else delta_eval=4/(double)DD;

 eval_f[level]=eval_f[level]+delta_eval;
eval_f_tot[level]=eval_f_tot[level]+delta_eval;

 // Option 1. Classical (discrete) way
  total=0;
for (d=0;d<DD-1;d++)
{
    i=(int)pos.p.x[d];
    j=(int)pos.p.x[d+1];
   total=total+problem[level].P.val[i][j];
 }
  i=(int)pos.p.x[DD-1];
  j=(int)pos.p.x[0];
 total=total+problem[level].P.val[i][j];

  if (penalty==0) goto end;

  // ... else, add penalty if several times the same node
 d=check_seq(pos);
 if (d>0) total=total+d*problem[level].P.max_cont;

 end:
 return total;
}

  //======================================================== INIT_PARTICLE_TSP
 struct particle init_particle_tsp(int size,double target,int level)
{
int                     d;

struct particle part;
struct position   x0;

x0.p.size=size;
for( d = 0;d<size;d++)  	x0.p.x[d] = d;// initial sequence

 x0.p=random_permut(x0.p);

 x0.f.f[0]=f_tsp(x0,1,0,-1,-1,level);   // Evaluate. We are _sure_ there is no penalty to add
 x0.f.f[0]=fabs(target-x0.f.f[0]);
 x0.f.size=problem[level].nb_f;
 part.x=x0;
 part.p_i=part.x; // Best previous = current
 part.prev_x=part.x; // Previous=current

 part.label=label[level];
 label[level]=label[level]+1;

 return part;

 }
 //======================================================== INIT_SWARM_TSP
struct tribe init_swarm_tsp(int size,double target,int trace,int level)
{
/* Initialization of the swarm */
int 				i;
int                     init_level;
int                     init_option;
struct   position     s1;
struct tribe   T;
double            total;


	printf("\n Swarm initialization for TSP");
	if (size>10) printf(" Please wait");

init_option=2; // *** Hard coded option
T.size=size;
 //s1.f.size=problem[level].nb_f;   // MUST be equal to 1, for the moment
  s1.f.size=1;

switch (init_option)
{
    case 1:  // Deterministic initialization, level 1
      for (i=0;i<size;i++)  // Define a tour beginning on i, using systematically
                                    // the smallest arc still free
      {
        if (trace>1) printf("\n init particle %i",i+1);

           T.part[i].x=min_tour_tsp(i,level);    // Tour, built from i, with a simple greedy algorithm
         total=f_tsp(T.part[i].x,1,0,-1,-1,level);   // Evaluate. We are _sure_ there is no penalty to add
         T.part[i].x.f.f[0]=fabs(target-total);
           //T.part[i].x.f.size=problem[level].nb_f; // Number of functions to optimize. MUST be equal to 1, for the moment
           T.part[i].x.f.size=1;
         T.part[i].p_i=T.part[i].x; // Best previous = current
         T.part[i].prev_x=T.part[i].x; // Previous = current
      }
      break;

 case 2:  // Deterministic initialization, level 2
 printf("\n init level 2. Size %i", size);
 init_level=1;
 case2:
      for (i=0;i<size;i++)  // Define a tour beginning on i, using systematically
                                    // the smallest 2-arcs path still free
      {
          // Tour, built from i, with a simple greedy algorithm
         s1=min_tour_tsp(i,level);
         total=f_tsp(s1,1,0,-1,-1,level);  // Evaluate. We are _sure_ there is no penalty to add
           s1.f.f[0]=fabs(target-total);
           //s1.f=T.part[i].x.f;

         T.part[i].x=min_tour_2_tsp(i,init_level,level);
         T.part[i].x.f.f[0]=f_tsp(T.part[i].x,1,0,-1,-1,level);
         T.part[i].x.f.f[0]=fabs(target-T.part[i].x.f.f[0]);

         if (T.part[i].x.f.f[0]>s1.f.f[0]) T.part[i].x=s1; // Choose the best
         //T.part[i].x.f.size=problem[level].nb_f;
         T.part[i].x.f.size=1;
          T.part[i].p_i=T.part[i].x; // Best previous = current
         T.part[i].prev_x=T.part[i].x; // Previous = current
      }
      break;

 case 3:
 printf("\n init level 2a");
 init_level=2;
 goto case2;

 }

 for (i=0;i<size;i++)    // Set the labels
 {
  T.part[i].label=label[level];
 label[level]=label[level]+1;
}

 printf("\n End of swarm initialization for TSP ");
 printf("\n");
return T;

}

//============================================================ LOCAL_SEARCH_TSP
struct position local_search_tsp(struct position pos,double target,int level)
{
/*
 Check the immediate neighbours, as long as there is some improvement
 Note: pos is supposed to be an integer position
*/
int                     DD;
int                     d1,d2;
double               f;
double               node;
struct position x1,x2;

 DD=pos.p.size;
 x1=pos;
 x2=pos;

 look_around:
 for (d1=0;d1<DD-1;d1++)
 {
      for (d2=d1+1;d2<DD;d2++)
      {
         // Transpose
         node=x2.p.x[d1];
         x2.p.x[d1]= x2.p.x[d2];
          x2.p.x[d2]=node;
        // Evaluate
        f=f_tsp(x2,1,1,d1,d2,level);
        f=fabs(target-f);
        if (f<x1.f.f[0])  // if improvement, memorize and try again
        {
 //printf("\nlocal_search_tsp. x1.f %f, f %f",x1.f.f[0],f);
          x1=x2;
          x1.f.f[0]=f;
          if((problem[level].Max_Eval<0 && eval_f_tot[level]>-problem[level].Max_Eval) || problem[level].Max_Eval>clock_tick)
            return x1;
          goto look_around;
        }
        else  // Transpose back and continue
        {
           node=x2.p.x[d1];
           x2.p.x[d1]= x2.p.x[d2];
           x2.p.x[d2]=node;
        }
      }
 }
//printf("\nlocal_search_tsp. %f => %f",pos.f.f[0],x1.f.f[0]);
return x1;
}
  //============================================================ MIN_TOUR _TSP
 struct   position   min_tour_tsp(int i,int level)
{
    //    Start from node i, and systematically add the smallest arc still free

 double        big_value;
  int             j;
  double       min;
  int             next;
  int             rank;
  double       val;
  struct position x;
  int          used[Max_C]={0};
 if(problem[level].printlevel>1) printf("\n min_tour_tsp node %i",i);

  big_value= problem[level].P.max_cont+1;
  x.p.size=problem[level].P.size;
  rank=0;
  x.p.x[rank]=i;
  used[i]=1;

  loop:
  min=big_value+1;

  for (j =0;j<x.p.size;j++) // Check all neighbours of x.s[rank]
  {
     if (used[j]==1) continue;
     val= problem[level].P.val[ (int)x.p.x[rank]][j];
      if (val<0) val=big_value;    // For non existent arc

      if (val<min)
      {
         next=j;
         min=val;
       }
   }
   rank=rank+1;
   x.p.x[rank]=next;
   used[next]=1;
   if (rank<x.p.size-1)  goto loop;

  x.p.x[x.p.size]=x.p.x[0]; // To complete the tour
  return x;
}

 //============================================================ MIN_TOUR_2_TSP
 struct   position   min_tour_2_tsp(int i,int init_level,int level)
{
    //    Start from node i, and systematically add the smallest 2-arcs path still free

 double        big_value;
  int             j,k;
  double       min;
  int             next,next2;
  int             rank;
  double       val_j,val_k;
  struct position x;
  int          used[Max_C]={0};
 if(problem[level].printlevel>1) printf("\n min_tour_2_tsp node %i",i);
   big_value= problem[level].P.max_cont+1;
  x.p.size=problem[level].P.size;
  rank=0;
  x.p.x[rank]=i;
  used[i]=1;

  loop:
 //printf("\n rank %i",rank);
  min=2*big_value+1;

  for (j =0;j<x.p.size;j++) // Check all neighbours of x.s[rank]
  {
 //printf("\n %i, %i, %f /",j,used[j],min);
     if (used[j]==1) continue;
     val_j= problem[level].P.val[(int) x.p.x[rank]][j];
      if (val_j<0) val_j=big_value;    // For non existent arc

      if (rank==problem[level].P.size-2)
      {
         x.p.x[x.p.size-1]=j;
         goto end;
      }

      for(k=0;k<x.p.size;k++)
      {
 //printf("\n %i",k);
      if (k==j) continue;
      if (used[k]==1) continue;
        val_k= problem[level].P.val[ j][k];
        if (val_k<0) val_k=big_value;

         if (val_j+val_k<min)
         {
            next=j;
            next2=k;
            min=val_j+val_k;
         }
      }
  }
  //printf("\n next %i",next);
  if(init_level==1)
  {
   rank=rank+1;
       x.p.x[rank]=next;
       used[next]=1;
       if (rank<x.p.size-1) goto loop;
  }

  if (init_level==2)
  {
       rank=rank+1;
       x.p.x[rank]=next;
       used[next]=1;

       rank=rank+1;
       x.p.x[rank]=next2;
       used[next2]=1;
       if (rank<x.p.size-2) goto loop;
  }


 end:
  x.p.x[x.p.size]=x.p.x[0]; // To complete the tour
  return x;
}

 //========================================================  NEAREST_INT
struct position nearest_int(struct position pos,double target,int level)
{
int                     d;
 struct position p;
 //printf("\n nearest int");
 //display_position(pos);
 p.p.size=pos.p.size;

 for (d=0;d<p.p.size;d++)
 {
    p.p.x[d]=floor(pos.p.x[d]+0.5);
    if (p.p.x[d]>=p.p.size) p.p.x[d]=p.p.size-1;
 }

 p.f.f[0]=f_tsp(p,1,1,-1,-1,level);
 p.f.f[0]=fabs(target-p.f.f[0]);
 return p;
}
 //========================================================  READ_GRAPH_TSP
struct matrix read_graph_tsp(int trace)
{
/* TSPLIB MATRIX format. One value/line. First value = number of nodes
each value = arc value. If <0 => no arc
*/
char			bidon[50];
char			comment[100];
char              graph[80];
char			edge_weight_format[30];
char			edge_weight_type[30];
struct matrix	Gt;
int 			i,j;
char			name[20];
double         total;
char			type[20];
float         zzz;

FILE *file;

printf("\nFile name for the graph, please: ");
gets(graph);

file=fopen(graph,"r");

printf("\nI am reading the graph");

fscanf(file," %s %s\n",bidon,name);
fscanf(file," %s %s\n",bidon,type);
	fscanf(file," %s %s\n",bidon,comment);
	fscanf(file,"%s %i\n",bidon,&Gt.size); // dimension
	fscanf(file,"%s %s\n",bidon,edge_weight_type);
	fscanf(file,"%s %s\n",bidon,edge_weight_format);
	fscanf(file,"%s\n",bidon);


	printf("\n Name: %s, type: %s, (%s)",name,type,comment);
	printf("\n Number of nodes: %i",Gt.size);
	printf("\n %s %s\n",edge_weight_type,edge_weight_format);

	if (edge_weight_type[0]=='E' && edge_weight_format[0]=='F')
		{

		for (i=0;i<Gt.size;i++)
			{
                for (j=0;j<Gt.size;j++)
				{
                        fscanf(file,"%e ",&zzz);
                            Gt.val[i][j]=zzz;
                            // A negative value is for "no edge". See below
					if (trace>2) printf(" %f",Gt.val[i][j]);
				}
			Gt.val[i][i]=0;
			}
            	printf("\n Weigths read. %i lines %i columns",Gt.size,Gt.size);

 Gt.max_cont=0;   //  Memorize the max value for future penalty
  		for (i=0;i<Gt.size;i++)
			{
                for (j=0;j<Gt.size;j++)
				{
                        if (Gt.val[i][j]>Gt.max_cont) Gt.max_cont=Gt.val[i][j];
                     }
                  }

 printf(" Max arc value %f",Gt.max_cont);

// Modify P to take unexistent edges into account
total=  1.1*Gt.max_cont*Gt.size; // Big value
  for (i=0;i<Gt.size;i++)
	{
		for (j=0;j<Gt.size;j++)
		{
            if(Gt.val[i][j]<0) Gt.val[i][j]=total;
         }
      }
		goto end;
		}
   printf("\nERROR by reading graph");
 end:
 fclose(file);
  return Gt;
}
