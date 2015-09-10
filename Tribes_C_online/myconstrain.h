struct position	all_different(struct position pos, float eps,int level ) ;
struct particle constrain(struct particle part,int level);
double         constrain_positive(double A);
double         constrain_null(double A);
struct position dichotom(struct position p0, struct position p1, int level);
struct position domain_pos(struct position pos,struct search_space H);
double   Myconstrain(struct position pos, int level);


