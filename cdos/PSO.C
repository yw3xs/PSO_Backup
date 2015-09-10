/*
  This program optimizes five benchmark functions using swarm algorithm
  Asynchronous version
  Yuhui Shi, May 15, 1998
*/

#include "headfile.h"
#include "global.h"
#include "mem_loc.h"
#include "myfun.h"



int c_break(void);

/* ******************** main() ******************** */
main (int argc, char *argv[])
{
  int  NUMBER_OF_AGENTS;

  float E_CUTOFF,  MAXV, MAXX;

  float weight, weight_up;
  int  MAXITER;
  int run_no;  //numbre of runs

  FILE *fp, *frun;
  char runfile[60], resfile[60];
  char temp[10];
  char tempfile[60];

  float minval=0.0;
  int DIMENSION;
  int fun_type;	//0:Schaffer f6	1:sphere		2:Rosenbrock	3:generalized Rastrigrin
						//4:generalized Griewank
  float IRang_L, IRang_R;  // initialization rang: left and right range
  int a,b;

  int i;

  int iter;
  int gbest;
  int firsttime;
  int tmp,finish;

  time_t tt;

	/* *********************************************************
		Open runfile
	********************************************************* */
	if (argc<2)
	{
		printf("Need to specify runfile. For Example: pso pso.run");
		exit(1);
	}
	strcpy(runfile,argv[1]);

	if ((frun=fopen(runfile,"r"))==NULL)
	{
		printf("Cant read file");
		exit(1);
	}

	/*                N  G  S  E  V  X  M  U  Res */
	fscanf(frun, "%d  %f  %f %f %f %f  %d  %s %f %d %d %d",
	 &NUMBER_OF_AGENTS,  &E_CUTOFF,
	 &MAXV, &MAXX,&IRang_L,&IRang_R, &MAXITER,  tempfile,
	 &weight, &fun_type,&DIMENSION,&run_no);
	fclose(frun);

	FVectorAllocate(&pbest,            NUMBER_OF_AGENTS);
	FVectorAllocate(&maxx,             DIMENSION);
	FMatrixAllocate(&vx,     DIMENSION,NUMBER_OF_AGENTS);
	FMatrixAllocate(&xx,     DIMENSION,NUMBER_OF_AGENTS);
	FMatrixAllocate(&tx,     DIMENSION,NUMBER_OF_AGENTS);
	FMatrixAllocate(&pbestx, DIMENSION,NUMBER_OF_AGENTS);

	for (a=0;a<DIMENSION;a++)
	{
		maxx[a]=MAXX;          /* range of xx[]  */
	}
	clrscr();
	randomize();

	//set control break
	ctrlbrk(c_break);


	time(&tt);
	printf("begin time: %s\n",ctime(&tt));

	//loop for runs
	for (i=0;i<run_no;i++)
	{
		firsttime=1;           //first iteration of this run
		iter=0;
		gbest=0;               //initialy assume the first particle as the gbest

		 /* **********************************************
			This loop initializes the individual agents  for each run
		********************************************** */
		for (a=0;a<NUMBER_OF_AGENTS;a++)
		{
			for (b=0;b<DIMENSION;b++)
			{
				xx[b][a] = (float) ((IRang_R - IRang_L)*(rand()/32767.0) + IRang_L);
				pbestx[b][a]=xx[b][a];
				vx[b][a] = MAXV*(rand()/32767.0);

				if ((rand()/32767.0) > 0.5) vx[b][a]=-vx[b][a];
			}
		}
		/* *******************************************************
			Main Work Loop for each run here
		******************************************************** */

		/********************************************************
			get the output result file name
		**********************************************************/
		if (i>=10)
		{
			int temdec=i/10;
			temdec=temdec+48;
			strcpy(temp,(char*)&temdec);
			tmp=i%10 +48;
			strcat(temp,(char*)&tmp);
		}
		else
		{
			int tmp=i+48;
			strcpy(temp,(char*)&tmp);
		}
		strcpy(resfile,tempfile);
		strcat(resfile,temp);
		strcat(resfile,".txt");

		/*****************************
			open file for output best agent index vs iteration
		*****************************************/
		if ((fp=fopen(resfile,"w"))==NULL)
		{
			printf("Cant write file");
			exit(1);
		}

		finish=0;

		do
		{
			iter++;
			if (iter >32760) iter=0;   /* so it doesnt crash the data type */

			//update inertia weight
			weight_up = (weight-0.4) * (MAXITER - iter) /MAXITER +0.4;    //time variant weight, linear from weight to 0.4

			//weight_up=weight;		//constant inertia weight

			for (a=0;a<NUMBER_OF_AGENTS;a++)
			{
				/* *********************************************
				eval(a) error is returned by function routines 
				********************************************* */
				switch (fun_type)
				{
					case 0:
						minval=f6(a);
						break;
					case 1:
						minval=sphere(a,DIMENSION);
						break;
					case 2:
						minval=rosenbrock(a,DIMENSION);
						break;
					case 3:
						minval=rastrigrin(a,DIMENSION);
						break;
					case 4:
						minval=griewank(a,DIMENSION);
						break;
					default:
						printf("\n Not a valid function type\n");
						exit(1);
				}

				if (firsttime==1) pbest[a]=minval;

				if (minval < pbest[a])
				{
					pbest[a]=minval;
					for (b=0;b<DIMENSION;b++) pbestx[b][a]=xx[b][a];
					if (pbest[a] < pbest[gbest])
					{
						gbest=a;

					}
				}

				/* asynchronous version */
				for (b=0;b<DIMENSION;b++)
				{
					vx[b][a] = weight_up*vx[b][a] + 2*(rand()/32767.0)*(pbestx[b][a]-xx[b][a]) +
					2*(rand()/32767.0)*(pbestx[b][gbest]-xx[b][a]);
					if (vx[b][a]>MAXV)
						vx[b][a]=MAXV;
					else if (vx[b][a]<-MAXV)
						vx[b][a]=-MAXV;
				}

				/*********************************
				Tx allows simultaneous updates
				*********************************/
				for (b=0; b<DIMENSION;b++)
				{
					tx[b][a]=xx[b][a]+vx[b][a];
				}
			}     	/* ******* END OF a LOOP ************************* */

			/* *******************************************************
				Update positions
			******************************************************** */
			for (a=0;a<NUMBER_OF_AGENTS;a++)
			{
				/* *******************************************************
					Define new coordinates
				******************************************************** */
				for (b=0;b<DIMENSION;b++)
				{
					xx[b][a] =tx[b][a];
				}

			}   /* end a loop  */

			/* *******************************************************
				In case iterations become greater than 32767
			******************************************************** */

			if (firsttime!=1)
			{
				if (iter==32766) iter=0;
			}
			/* ******************************************************
				Terminate on criterion
			********************************************************* */
			//output best index vs iteration
			if (fun_type==0)
			{
				fprintf(fp,"%f\n",1.0-pbest[gbest]);

			}
			else
			{
				fprintf(fp,"%f\n",pbest[gbest]);

			}

			if ((pbest[gbest] <= E_CUTOFF) || (iter >= MAXITER))
			{
				fclose(fp);
				printf("%d run finished!\n",i);
				finish=1;
			}
			firsttime=0;

		}     /* **************** End of do-loop *************** */
		while (! finish);
	}
	time(&tt);
	printf("end time: %s\n",ctime(&tt));
	return 0;
}




int c_break(void)
{
	closegraph();
	return 0;
}




