void read_vtk(struct Field fldo, const int n, double t);
void remap_output_back(	double wri[], 
					const double t);
void remap_output(	double wri[], 
					const double t);
void write_vtk(double complex *wi, 
			   const char filename[], 
			   const char fieldname[], 
			   const double t);
double read_timevar(const int n);
void init_output();
void finish_output();