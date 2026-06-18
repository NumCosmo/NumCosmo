int
pdip_initial_point_feasible_x(const gsl_matrix * A, const gsl_vector *b, gsl_vector *x);

int
pdip_initial_point_feasible_s(const gsl_matrix * C, const gsl_vector *d, const gsl_vector *x, gsl_vector *s);

int
pdip_initial_point_y(const gsl_matrix *Q, const gsl_vector *q, const gsl_matrix *A, const gsl_vector *x, gsl_vector *y);

int
pdip_initial_point_z(gsl_vector *z);

int
pdip_initial_point_strict_feasible(gsl_vector *x, gsl_vector *s);
