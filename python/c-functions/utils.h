void distance(double* p1, double* p2, double* box, double* out);
void periodic_conditions(double* delta, double* box);
void norm_vec(double* vec);
void angle(double* v1, double* v2, double* theta);
void closest_atom(double* names, double* oxygens, double* hydrogens,
                  double* box, double* out_name);
