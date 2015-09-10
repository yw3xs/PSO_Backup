int check_seq(struct position s);
double f_tsp(struct position pos, int option, int penalty,int rank1,int rank2,int level);
struct particle init_particle_tsp(int size,double target,int level);
struct tribe init_swarm_tsp(int size,double target,int trace,int level);
struct position local_search_tsp(struct position x,double target,int level);
struct   position   min_tour_tsp(int i,int level);
struct   position   min_tour_2_tsp(int i,int init_level,int level);
struct position nearest_int(struct position pos,double target,int level);
struct matrix  read_graph_tsp(int trace);

