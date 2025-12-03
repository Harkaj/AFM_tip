#ifndef AFM_TIP_H
#define AFM_TIP_H

#define N 500
#define PI 3.14159
#define F0 0.04
#define Fg 0.00003
#define Fo 0.005
#define RR 15
#define k 0.5
#define mu 0.02
#define dt 0.005
#define STUCK 1
#define FREE 0

extern int itmax, i_out, create_new, stable_state, it_start;
extern double c, p, l, S, c_x, c_y, c_y0, c_R, F_c;
extern double N_x[N], N_y[N], v_x[N], v_y[N], F_x[N], F_y[N], N_c[2], F[2];
extern int stuck[N];

void Read_file(const char *name, double N_x[], double N_y[]);
void Create_points(double N_x[], double N_y[], double c);
int Determine_neighbour(int i, int direction);
double Collision_correction(double N_x[], double N_y[], double v_x[], 
                            double v_y[], double c_x, double c_y, 
                            double c_R, double F_x[], double F_y[], int output);
void Calculate_mass_centre(double N_x[], double N_y[], double N_c[]);
double Calculate_distance(double x_1, double y_1, double x_2, double y_2);
double Calculate_vector_angle(double x_1, double y_1, double x_2, double y_2);
double Calculate_surface_circumference(double N_x[], double N_y[]);
double Calculate_surface_area(double N_x[], double N_y[]);
void Calculate_forces(double N_x[], double N_y[], double F_x[], double F_y[], double S, double l, int i);
void Calculate_forces_damped(double N_x[], double N_y[], double v_x[], double v_y[], double F_x[], double F_y[], double S, double l, int i);
void Calculate_total_force(double F_x[], double F_y[], double F[]);
double Calculate_next_step(double r[], double v[], double F[]);
void POVRAY_output_ini(int itmax, int i_out);
void POVRAY_output(double N_x[], double N_y[], int it, int itmax);
void GNUplot_output_no_border(char name[], int column1, int column2);
void GNUplot_output_border(char name[], char out_name[], int column1, int column2);
void Text_output(double N_x[], double N_y[], double F_x[], double F_y[], char name[]);

#endif
