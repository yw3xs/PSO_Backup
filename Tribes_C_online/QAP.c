
 //======================================================== F_QAP
 double f_qap(struct position pos,int option, int penalty,int rank1,int rank2,int level)
 {
 /*
   The matrix P   is a global variable. It contains the flow matrix  and the distance matrix
   problem is a global variable
   If rank1>=0, positions in rank1 and rank2 have just been switched
 */
 int        d;
 double  delta_eval;
  int        DD;
 int        i,j;
 int                    p_i,p_j;
 double total;

   DD=problem[level].P.size;
 /*
   When there has been just a switch between positions, the f_value could be (and will be in next version)
   evaluated by updating  at most 4 arcs values
 */
 if (rank1<0) delta_eval=1; else delta_eval=4/(double)DD;

eval_f[level]=eval_f[level]+delta_eval;
eval_f_tot[level]=eval_f_tot[level]+delta_eval;

// if (option==0) goto continuous;   // Not yet written

 // Option 1. Classical (discrete) way
  total=0;
for (i=0;i<DD;i++)
{
	p_i=(int)pos.p.x[i]+DD;
    for(j=0;j<DD;j++)
    {
       p_j=(int)pos.p.x[j];
       total=total+problem[level].P.val[i][j]*problem[level].P.val[p_i][p_j];
   }
 }


  if (penalty==0) goto end; // the penalty will be taken into account in continuous mode

  // ... else, add penalty if several times the same node
 d=check_seq(pos);
 if (d>0) total=total+d*problem[level].P.max_cont;
 goto end;

 /*
 //----------------------------------------- (to rewrite)
 continuous:  // Option 0
 test=1;
  switch(test)
 {
    case 1:
    p1=nearest_int(pos);
     d1=distance(pos,p1);
     big_value=1.1*DD*P.max_cont;
     max_d=sqrt(DD)/2;
     f1=p1.f+problem[recursive].target;
     total=  f1+d1*(big_value-f1)/max_d;
 //printf("\nf_tsp. f1 %f, d1 %f,  total %f",f1, d1,total);
     break;

 case 2:
  p1=nearest_int(pos);
  p2.p.size=DD;
   for (d=0;d<DD;d++)
 {
     if (pos.p.x[d]>p1.p.x[d]) {p2.p.x[d]=p1.p.x[d]+1;continue;}
     if (pos.p.x[d]<p1.p.x[d]) {p2.p.x[d]=p1.p.x[d]-1;  continue;}
     // pos.p.x[d] is already an integer
     if (pos.p.x[d]>1) p2.p.x[d]=p1.p.x[d]-1; else p2.p.x[d]=p1.p.x[d]+1;
 }
 break;

 }

  //printf("\n f: %f %f, d: %f %f, total %f",p1.f, p2.f,d1,d2,total);

  //Add penalty if values are not different enough
  for (i=0;i<DD-1;i++)
  {
    for (j=i+1;j<DD;j++)
    {
     d1=fabs(pos.p.x[i]-pos.p.x[j]);
      if (d1<1) total=total+P.max_cont*(1-d1);
    }
  }
 */

 end:
 //printf("\nf_qap. total %f",total);
 return total;
}

  //======================================================== INIT_PARTICLE_QAP
 struct particle init_particle_qap(int size,double target,int level)
{
int                     d;
struct particle part;
struct position   x0;

x0.p.size=size;
for( d = 0;d<size;d++)  // initial sequence
{
	x0.p.x[d] = d;
}

 x0.p=random_permut(x0.p);

 x0.f.f[0]=f_qap(x0,1,0,-1,-1,level);   // Evaluate. We are _sure_ there is no penalty to add
 x0.f.f[0]=fabs(target-x0.f.f[0]);
 part.x=x0;
 part.p_i=part.x; // Best previous = current
 part.prev_x=part.x; // Previous = current

 part.label=label[level];
 label[level]=label[level]+1;

 return part;

 }
 //======================================================== INIT_SWARM_QAP
struct tribe init_swarm_qap(int size,double target,int trace,int level)
{
/* Initialization of the swarm */

int 				i;
int                     init_level;
int                     init_option;
struct   position     s1;
struct tribe   T;
double            total;


	printf("\n Swarm initialization for QAP");
	if (size>10) printf(" Please wait");

init_option=2; // *** Hard coded
                        // 1 => use min_tour
                        // 2 => use best of min_tour, min_tour_2
T.size=size;
s1.f.size=problem[level].nb_f;

switch (init_option)
{
    case 1:  // Deterministic initialization, level 1
      for (i=0;i<size;i++)  // Define a permutation beginning on i, using systematically
                                    // the best next node still free
      {
        if (trace>1) printf("\n init particle %i",i+1);

           T.part[i].x=min_tour_qap(i,level);
         total=f_qap(T.part[i].x,1,0,-1,-1,level);   // Evaluate. We are _sure_ there is no penalty to add
         T.part[i].x.f.f[0]=fabs(target-total);
         T.part[i].x.f.size=problem[level].nb_f;
         T.part[i].p_i=T.part[i].x; // Best previous = current
         T.part[i].prev_x=T.part[i].x; // Previous = current
//display_position( T.part[i].x);
      }
      break;

 case 2:  // Deterministic initialization, level 2
 printf("\n init level 2");
 init_level=1;
 case2:
      for (i=0;i<size;i++)  // Define a permatutation beginning on i, using systematically
                                    // the smallest 2-arcs path still free
      {
         s1=min_tour_qap(i,level);
         total=f_qap(s1,1,0,-1,-1,level);  // Evaluate. We are _sure_ there is no penalty to add
           s1.f.f[0]=fabs(target-total);

         T.part[i].x=min_tour_2_qap(i,init_level,level);
         T.part[i].x.f.f[0]=f_qap(T.part[i].x,1,0,-1,-1,level);
         T.part[i].x.f.f[0]=fabs(target-T.part[i].x.f.f[0]);

         if (T.part[i].x.f.f[0]>s1.f.f[0]) T.part[i].x=s1; // Choose the best
                  T.part[i].x.f.size=problem[level].nb_f;
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
 printf("\n End of swarm initialization for QAP ");

return T;

}

//============================================================ LOCAL_SEARCH_QAP
struct position local_search_qap(struct position x,double target,int level)
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

 DD=x.p.size;
 x1=x;
 x2=x;

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
        f=f_qap(x2,1,1,d1,d2,level);
        f=fabs(target-f);
        if (f<x1.f.f[0])  // if improvement, memorize and try again
        {
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
//printf("\nlocal_search_qap. %f => %f",pos.f,x1.f);
return x1;
}
  //============================================================ MIN_TOUR _QAP
 struct   position   min_tour_qap(int i,int level)
{
    //    Start from node i, and systematically add the smallest arc still free

 double        big_value;
 double        delta_val;
 int              DD;
  int             j;
  double       min;
  int             next;
  int             r;
  int             rank;
  double       val;
  struct position x;
  int          used[Max_C]={0};

  DD=problem[level].P.size;
  big_value= DD*DD*problem[level].P.max_cont+1;
  x.p.size=DD;
  rank=0;
  x.p.x[rank]=i;
  used[i]=1;

  loop:
  min=big_value+1;

  for (j =0;j<DD;j++) // Check all neighbours of x.p.x[rank]
  {
     if (used[j]==1) continue;
     val=0;
     for (r=0;r<rank;r++)
     {
        delta_val=problem[level].P.val[r][rank+1]*problem[level].P.val[ (int)x.p.x[r]+DD][j];
        val= val+delta_val;
     }

      if (val<min)
      {
         next=j;
         min=val;
       }
   }

   rank=rank+1;
   x.p.x[rank]=next;
   used[next]=1;
   if (rank<DD-1)  goto loop;
 //display_position(x);
  return x;
}

 //============================================================ MIN_TOUR_2_QAP
 struct   position   min_tour_2_qap(int i,int init_level,int level)
{
    //    Start from node i, and systematically add the smallest 2-arcs path still free

 double        big_value;
 int              DD;
 double        delta_val;
  int             j,k;
  double       min;
  int             next,next2;
  int             r;
  int             rank;
  double       val_k;
  struct position x;
  int          used[Max_C]={0};
  double    val;

   DD=problem[level].P.size;
  big_value= DD*DD*problem[level].P.max_cont+1;
  x.p.size=DD;
  rank=0;
  x.p.x[rank]=i;
  used[i]=1;

  loop:
 //printf("\n rank %i",rank);
  min=2*big_value+1;

  for (j =0;j<DD;j++) // Check all neighbours of x.p.x[rank]
  {
 //printf("\n %i, %i, %f /",j,used[j],min);
     if (used[j]==1) continue;
        val=0;
     for (r=0;r<rank;r++)
     {
        delta_val=problem[level].P.val[r][rank+1]*problem[level].P.val[ (int)x.p.x[r]+DD][j];
        val= val+delta_val;
     }


      if (rank==DD-2)
      {
         x.p.x[DD-1]=j;
         goto end;
      }

      for(k=0;k<DD;k++)
      {
 //printf("\n %i",k);
      if (k==j) continue;
      if (used[k]==1) continue;
        val_k= problem[level].P.val[ j][k];
        if (val_k<0) val_k=big_value;

         if (val+val_k<min)
         {
            next=j;
            next2=k;
            min=val+val_k;
         }
      }
  }
  //printf("\n next %i",next);
  if(level==1)
  {
   rank=rank+1;
       x.p.x[rank]=next;
       used[next]=1;
       if (rank<DD-1) goto loop;
  }

  if (level==2)
  {
       rank=rank+1;
       x.p.x[rank]=next;
       used[next]=1;

       rank=rank+1;
       x.p.x[rank]=next2;
       used[next2]=1;
       if (rank<DD-2) goto loop;
  }


 end:
  x.p.x[x.p.size]=x.p.x[0]; // To complete the tour
  return x;
}


 //========================================================  READ_GRAPH_QAP
struct matrix read_graph_qap(int trace)
{
/* Read data for Quadratic Assignment Problem
<dimension d>
<d*d flow matrix>
<d*d distance matrix>
Diagonal values are supposed to be null
*/

double            big_value;
char              graph[80];
struct matrix	Gt;  // Contains the two matrices
int 			i,j;
float         zzz;

FILE *file;

printf("\nFile name for the QAP matrices, please: ");
gets(graph);

file=fopen(graph,"r");

printf("\nI am reading the matrices");

fscanf(file,"%i\n",&Gt.size); // dimension

		for (i=0;i<2*Gt.size;i++)
			{
                for (j=0;j<Gt.size;j++)
				{
                        fscanf(file,"%f",&zzz);
                            Gt.val[i][j]=zzz;
					if (trace>2) printf(" %f",Gt.val[i][j]);
				}
			}
            	printf("\n Flows and distances read. %i lines %i columns",2*Gt.size,Gt.size);

   big_value=0;

  		for (i=Gt.size;i<2*Gt.size;i++)
			{
                for (j=0;j<Gt.size;j++)
				{
                        if (Gt.val[i][j]>big_value) big_value=Gt.val[i][j];
                     }
                  }
   // Replace the negative arc values   (unexistent arcs)

    		for (i=Gt.size;i<2*Gt.size;i++)
			{
                for (j=0;j<Gt.size;j++)
				{
                        if (Gt.val[i][j]<0) Gt.val[i][j]=1.1*(Gt.size*big_value);
                     }
                  }

Gt.max_cont=big_value;
big_value=0;   //  Memorize a big value value for future penalty
    		for (i=0;i<Gt.size;i++)
			{
                for (j=0;j<Gt.size;j++)
				{
                        if (Gt.val[i][j]>big_value) big_value=Gt.val[i][j];
                     }
                  }

Gt.max_cont=big_value*  Gt.max_cont;
 printf(" Max flow*distance value %f",Gt.max_cont);

 fclose(file);
  return Gt;
}
