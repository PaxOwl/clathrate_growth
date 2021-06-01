void distance(double *p1, double *p2, double *box, double *vec);
void periodic_conditions(double *delta, double *box);
void norm_vec(double *vec);
void angle(double *v1, double *v2, double *theta);
void closest_atom(double *names, double *oxygens, double *hydrogens,
                  double *box, double *out_name);
void nearest_neighbours(double *center, double (*neighbours)[3], double *box,
                        double h_lim, double l_lim,
                        int n_size, double *vec, long *out);
void neighbours(double (*centers)[3], double (*neighbours)[3], double *box,
                double h_lim, double l_lim, int c_size, int n_size,
                double *vec, long (*out)[n_size]);
void aop(double *center, double (*neighbours)[3]);
