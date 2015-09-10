double f_qap(struct position pos, int option, int penalty,int rank1,int rank2,int level);
struct particle init_particle_qap(int size,double target,int level);
struct tribe init_swarm_qap(int size,double target,int trace,int level);
struct position local_search_qap(struct position x,double target,int level);
struct   position   min_tour_qap(int i,int level);
struct   position   min_tour_2_qap(int i,int init_level,int level);
struct matrix  read_graph_qap(int trace);

