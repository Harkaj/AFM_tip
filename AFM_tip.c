#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "AFM_tip.h"

int itmax, i_out, create_new, stable_state, it_start;
double c, p, l, S, c_x, c_y, c_y0, c_R, F_c;
double N_x[N], N_y[N], v_x[N], v_y[N], F_x[N], F_y[N], N_c[2], F[2];

FILE *out;

int main(int argc, char** argv[])
{
    c = 15.1;
    c_x = 0.0;
    c_y0 = 30.0;
    c_y = c_y0;
    c_R = 3;
    itmax = 2000000;
    i_out = 10000;
    create_new = 0;
    stable_state = 1;

    if (create_new) 
    { 
        Create_points(N_x, N_y, c); 
        it_start = 1000000; 
    }
    else 
    { 
        Read_file("Shape_unstuck.txt", N_x, N_y); 
        it_start = 0; 
    }

    POVRAY_output_ini(itmax, i_out);

    out = fopen("Evolution.txt", "w");
    fprintf(out, "i Fx     Fy    S       l      dy_c  F_c\n");

    for (int i = 0; i < itmax; i++)
    {
        Calculate_mass_centre(N_x, N_y, N_c);
        l = Calculate_surface_circumference(N_x, N_y);
        S = Calculate_surface_area(N_x, N_y);

        for (int j = 0; j < N; j++)
        {
            //if (N_y[j] <= 0.0) { continue; }
            //Calculate_forces(N_x, N_y, F_x, F_y, S, l, j);
            Calculate_forces_damped(N_x, N_y, v_x, v_y, F_x, F_y, S, l, j);
        }

        Calculate_next_step(N_x, v_x, F_x);
        Calculate_next_step(N_y, v_y, F_y);

        Collision_correction_floor(N_y, v_y);
        F_c = Collision_correction_circle(N_x, N_y, v_x, v_y, c_x, c_y, c_R, F_y);
        if (i > it_start && i % 100 == 0 && c_y > 10)
        {
            c_y -= 0.001;
            if (stable_state) 
            {
                Text_output(N_x, N_y, F_x, F_y, "Shape_stable.txt");
                stable_state = 0;
            }
        }

        Calculate_total_force(F_x, F_y, F);
        if (i % 1000 == 0)
        {
            fprintf(out, "%d %.3f %.3f %.3f %.3f %.3f %.3f\n", i, F[0], F[1], S, l, c_y0 - c_y, F_c);
        }

        if (i % i_out == 0)
        {
            POVRAY_output(N_x, N_y, i / i_out, itmax / i_out);
        }
    }
    
    fclose(out);

    POVRAY_output(N_x, N_y, itmax / i_out, itmax / i_out);
    Text_output(N_x, N_y, F_x, F_y, "Shape.txt");
    GNUplot_output_border("Evolution.txt", "Evolution", 1, 4);
    GNUplot_output_border("Evolution.txt", "Forces", 6, 7);

    system("gnuplot -p Shape_GNU.gp");
    system("gnuplot -p Evolution_GNU.gp");
    system("gnuplot -p Forces_GNU.gp");
}


