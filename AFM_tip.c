#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "AFM_tip.h"

int itmax, i_out, create_new, stable_state, it_start, snaps;
double c, p, l, S, c_x, c_y, c_y0, c_R, F_cx, F_cy;
double N_x[N], N_y[N], v_x[N], v_y[N], F_x[N], F_y[N], N_c[2], F[2];
int stuck[N];
char snap_name[100];

FILE *out;

int main(int argc, char** argv[])
{
    c = 15.1;
    c_x = 0.0;
    c_y0 = 40.0;
    c_y = c_y0;
    c_R = 5;
    itmax = 5000000;
    i_out = 10000;
    create_new = 0;
    stable_state = 1;
    snaps = 1;

    if (create_new) 
    { 
        Create_points(N_x, N_y, c); 
        it_start = 2000000; 
    }
    else 
    { 
        Read_file("Shape_unstuck.txt", N_x, N_y); 
        it_start = 0; 
    }

    POVRAY_output_ini(itmax, i_out);

    out = fopen("Evolution.txt", "w");
    fprintf(out, "       i       Fx       Fy        S        l     dy_c     F_cx     F_cy\n");

    GNUplot_output_no_border_loop(1, 2);
    system("rm Snapshots/*.txt");
    system("rm Snapshots/GNUplot/*.png");

    for (int i = 0; i < itmax; i++)
    {
        if (i > it_start && i % 20 == 0 && c_y > 10)
        {
            c_y -= 0.0001;//= c_y0 - (i - it_start) * 0.000002;// -= 0.0001;
            if (c_y < 35.0 && i % 20000 == 0) 
            {
                sprintf(snap_name, "Snapshots/Shape%d.txt", snaps); 
                Text_output(snap_name, N_x, N_y, F_x, F_y);
                snaps++;
            }
        }
        
        Calculate_mass_centre(N_x, N_y, N_c);
        l = Calculate_surface_circumference(N_x, N_y);
        S = Calculate_surface_area(N_x, N_y);

        for (int j = 0; j < N; j++) { Calculate_forces_damped(N_x, N_y, v_x, v_y, F_x, F_y, S, l, j); }

        Calculate_next_step(N_x, v_x, F_x);
        Calculate_next_step(N_y, v_y, F_y);
        F_cy = Collision_correction(N_x, N_y, v_x, v_y, c_x, c_y, c_R, F_x, F_y, 1);
        
        Calculate_total_force(F_x, F_y, F);

        if (i % 1000 == 0)
        {
            fprintf(out, "%8d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", i, F[0], F[1], S, l, c_y0 - c_y, F_cx, F_cy);
        }

        if (i % i_out == 0) { /*POVRAY_output(N_x, N_y, i / i_out, itmax / i_out);*/ }
    }
    
    fclose(out);

    POVRAY_output(N_x, N_y, itmax / i_out, itmax / i_out);
    Text_output("Shape.txt", N_x, N_y, F_x, F_y);
    GNUplot_output_no_border("Shape.txt", 1, 2, c_x, c_y, c_R);
    GNUplot_output_border("Evolution.txt", "Evolution", 1, 4);
    GNUplot_output_border("Evolution.txt", "Forces", 6, 8);

    system("gnuplot -p Shape_GNU.gp");
    system("gnuplot -p Evolution_GNU.gp");
    system("gnuplot -p Forces_GNU.gp");
}


