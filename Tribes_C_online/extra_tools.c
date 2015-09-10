// ----------------------------------------------------------------------------- D_MIN
double d_min(struct position pos, int level)
{
/*
Compute the minimal distance of the position pos to the "box"
(search space of the problem, global variable)
*/
int		d;
double	dmin;
double	md;
double	z;

dmin=fabs(pos.p.x[0]-problem[level].H.min[0]);

for (d=0;d<pos.p.size;d++)
{
	md=fabs(pos.p.x[d]-problem[level].H.min[d]);
	z=fabs(pos.p.x[d]-problem[level].H.max[d]);
	if (z<md) md=z;
	if(md<dmin) dmin=md;
}

return dmin;
}

// ----------------------------------------------------------------------------- DISTANCE
double	distance(struct position pos1,struct position pos2)
{
double 	d,dx;
int		i;
d=0;
for (i=0;i<pos1.p.size;i++)
{
	dx=pos1.p.x[i]-pos2.p.x[i];
	d=d+dx*dx;
}
d=sqrt(d);
return d;
}

// ----------------------------------------------------------------------------- ENERGY_K
double energy_k(int level)
{
/*
Compute the kinetic energy of the swarm
The tribe list TR is a global variable
*/
int      d;
int      DD;
double Ek;
int      i;
int      j;
double vd;

DD=TR[level].tr[0].part[0].x.p.size;

 Ek=0;
 for (i=0;i<TR[level].size;i++)  // For each tribe
 {
      for (j=0;j<TR[level].tr[i].size;j++)  // for each particle in the tribe
      {
         for (d=0;d<DD;d++)  // Velocity component
         {
            vd= TR[level].tr[i].part[j].x.p.x[d]- TR[level].tr[i].part[j].prev_x.p.x[d] ;
            Ek=Ek+vd*vd;
         }
      }

 }

return Ek/2 ;
}


// ----------------------------------------------------------------------------- ENERGY_P
double energy_p(int level)
{
/*
Compute the potential energy of the swarm
The tribe list TR and "problem" are global variable
*/
int      d;
int      DD;
double Ep;
int      i;
int      j;
double vd;

DD=TR[level].tr[0].part[0].x.p.size;

 Ep=0;
 for (i=0;i<TR[level].size;i++)  // For each tribe
 {
      for (j=0;j<TR[level].tr[i].size;j++)  // for each particle in the tribe
      {
		Ep=Ep+total_error(TR[level].tr[i].part[j].x.f);
      }

 }

return Ep/ problem[level].eps;
}

// ----------------------------------------------------------------------------- DISTANCES_PAIRS
void distances_pairs(int level)
{
/*
    Histogram of distances between pairs of particules.  If they are random, it is "flat"
    TO COMPLETE
*/

}

// ----------------------------------------------------------------------------- INFORMATIVE_LINKS
void informative_links (int level)
{
 /*
 Compute some charactristics of the informative links graph:
 - diameter
 - aggregation coefficient
 etc.

 In particular, try to see if the graph is indeed a "small world".
   TO COMPLETE
 */


}


