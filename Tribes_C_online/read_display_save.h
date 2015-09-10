void			display_i_group(FILE *f_trace,struct particle part,int option);
void              display_memo(int level);
void			display_problem(int level);
void			display_position(struct position pos);
void			display_swarm(int type,int level);
void			display_tribe(FILE *f_trace,struct tribe tr, int option,int level);
void 			print_N_run(int level,struct position best_result);
struct problem  read_problem(FILE *f_param, FILE *f_data,int level);
void			save_swarm(FILE *f_swarm,int level);

