#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "AFM_tip.h"

void Read_file(const char *name, double N_x[], double N_y[])
{
    char line[100], *stop;
    FILE *in = fopen(name, "r");

    int i = 0;
    while (fgets(line, sizeof(line), in))
    {
        char *token = strtok(line, " ");
        N_x[i] = strtod(token, &stop);
        token = strtok(NULL, " ");
        N_y[i] = strtod(token, &stop);
        i++;
    }
    fclose(in);
}

void Create_points(double N_x[], double N_y[], double c)
{
    for (int i = 0; i < N; i++)
    {
        N_x[i] = - RR * cos(2.0 * PI * i / N);
        N_y[i] = c - RR * sin(2.0 * PI * i / N);
    }
}

int Determine_neighbour(int i, int direction)
{    
    int minus = i - 1;
    int plus = i + 1;
    if (i == 0) { minus = N - 1; }
    if (i == N - 1) { plus = 0; }

    if (direction < 0) { return minus; } 
    else if (direction > 0) { return plus; }
}

double Collision_correction(double N_x[], double N_y[], double v_x[], 
                            double v_y[], double c_x, double c_y, 
                            double c_R, double F_x[], double F_y[], int output)
{
    double dr, phi, result;

    result = 0.0;
    for (int i = 0; i < N; i++)
    {
        if (N_y[i] <= 0.0)
        {
            N_y[i] = 0.0;
            F_y[i] = 0.0;
            v_x[i] = 0.0;
            v_y[i] = 0.0;
            stuck[i] = STUCK;
            continue;
        }

        dr = Calculate_distance(c_x, c_y, N_x[i], N_y[i]);
        if (dr > c_R + 0.00001) { stuck[i] = FREE; continue; }

        phi = Calculate_vector_angle(c_x, c_y, N_x[i], N_y[i]);
        v_x[i] = 0.0;
        v_y[i] = 0.0;
        
        while (dr < c_R)
        {
            if (N_x[i] < c_x && F_x[i] <= 0.0) { N_x[i] -= 0.000002; }
            else if (F_x[i] > 0.0) { N_x[i] += 0.000002; }
            
            N_y[i] -= 0.00002;
            dr = Calculate_distance(c_x, c_y, N_x[i], N_y[i]);
        }
        stuck[i] = STUCK;
        if (output == 0) { result += F_x[i];}
        else { result += F_y[i]; }
    }

    return result;
}

void Calculate_mass_centre(double N_x[], double N_y[], double N_c[])
{
    N_c[0] = 0.0;
    N_c[1] = 0.0;

    for (int i = 0; i < N; i++)
    {
        N_c[0] += N_x[i];
        N_c[1] += N_y[i];
    }
    
    N_c[0] += N_c[0] / N;
    N_c[1] += N_c[1] / N;
}

double Calculate_distance(double x_1, double y_1, double x_2, double y_2)
{
    double dx, dy, result;

    dx = x_2 - x_1;
    dy = y_2 - y_1;

    result = dx * dx + dy * dy;
    result = sqrt(result);

    return result;
}

double Calculate_vector_angle(double x_1, double y_1, double x_2, double y_2)
{
    double dx, dy;
    double result = 0.0;

    dx = x_2 - x_1;
    dy = y_2 - y_1;

    result = atan2(dy, dx);
    
    return result;
}

double Calculate_surface_circumference(double N_x[], double N_y[])
{
    double dr;
    double result = 0.0;

    for (int i = 1; i < N; i++)
    {
        dr = Calculate_distance(N_x[i - 1], N_y[i - 1], N_x[i], N_y[i]);
        result += dr;
    }

    dr = Calculate_distance(N_x[N - 1], N_y[N - 1], N_x[0], N_y[0]);
    result += dr;

    return result;
}

double Calculate_surface_area(double N_x[], double N_y[])
{
    int min_i, max_i, im, ip, which_sum;
    double min, max, sum1, sum2, ds, result;

    min = N_x[0];
    max = N_x[0];
    min_i = 0;
    max_i = 0;
    for (int i = 1; i < N; i++)
    {
        if (N_x[i] < min) { min = N_x[i]; min_i = i; }
        if (N_x[i] > max) { max = N_x[i]; max_i = i; }
    }

    sum1 = 0.0;
    sum2 = 0.0;
    which_sum = 1;
    for (int i = 0; i < N; i++)
    {
        if (i == min_i) { which_sum *= -1; }
        if (i == max_i) { which_sum *= -1; } 
        
        ds = (N_y[i + 1] + N_y[i]) * (N_x[i + 1] - N_x[i]) / 2.0;
        if (i == N - 1) { ds = (N_y[0] + N_y[i]) * (N_x[0] - N_x[i]) / 2.0; }

        if (which_sum > 0) { sum1 += abs(ds); } 
        if (which_sum < 0) { sum2 += abs(ds); } 
    }

    result = abs(sum1 - sum2);
    return result;
}

void Calculate_forces(double N_x[], double N_y[], double F_x[], double F_y[], double S, double l, int i)
{
    double d_m, d_p, phi_m, phi_p, phi_n;
    double result = 0.0;

    int im = Determine_neighbour(i, -1);
    int ip = Determine_neighbour(i, 1);

    d_m = Calculate_distance(N_x[i], N_y[i], N_x[im], N_y[im]);
    d_p = Calculate_distance(N_x[i], N_y[i], N_x[ip], N_y[ip]);
    phi_m = Calculate_vector_angle(N_x[i], N_y[i], N_x[im], N_y[im]);
    phi_p = Calculate_vector_angle(N_x[i], N_y[i], N_x[ip], N_y[ip]);

    phi_n = Calculate_vector_angle(N_x[im], N_y[im], N_x[ip], N_y[ip]);
    phi_n -= PI / 2.0;

    F_x[i] = k * (d_p * cos(phi_p) + d_m * cos(phi_m)) + (F0 * l * cos(phi_n))/ S;
    F_y[i] = k * (d_p * sin(phi_p) + d_m * sin(phi_m)) + (F0 * l * sin(phi_n))/ S -  Fg;
}

void Calculate_forces_damped(double N_x[], double N_y[], double v_x[], double v_y[], double F_x[], double F_y[], double S, double l, int i)
{
    double d_m, d_p, phi_m, phi_p, phi_n;
    double result = 0.0;

    int im = Determine_neighbour(i, -1);
    int ip = Determine_neighbour(i, 1);

    d_m = Calculate_distance(N_x[i], N_y[i], N_x[im], N_y[im]);
    d_p = Calculate_distance(N_x[i], N_y[i], N_x[ip], N_y[ip]);
    phi_m = Calculate_vector_angle(N_x[i], N_y[i], N_x[im], N_y[im]);
    phi_p = Calculate_vector_angle(N_x[i], N_y[i], N_x[ip], N_y[ip]);

    phi_n = Calculate_vector_angle(N_x[im], N_y[im], N_x[ip], N_y[ip]);
    phi_n -= PI / 2.0;

    F_x[i] = k * (d_p * cos(phi_p) + d_m * cos(phi_m)) + (F0 * l * cos(phi_n))/ S - (Fo * cos(phi_n)) - v_x[i] * mu;
    F_y[i] = k * (d_p * sin(phi_p) + d_m * sin(phi_m)) + (F0 * l * sin(phi_n))/ S - (Fo * sin(phi_n)) - v_y[i] * mu - Fg;
}

void Calculate_total_force(double F_x[], double F_y[], double F[])
{
    F[0] = 0.0;
    F[1] = 0.0;
    for (int i = 0; i < N; i++)
    {
        F[0] += F_x[i];
        F[1] += F_y[i];
    }
}

double Calculate_next_step(double r[], double v[], double F[])
{
    for (int i = 0; i < N; i++)
    {
        if (stuck[i] == STUCK) { continue; }
        v[i] += F[i] * dt;// - 0.01 * abs(v[i]) / v[i];
        r[i] += v[i] * dt;
    }
}

void POVRAY_output_ini(int itmax, int i_out)
{
    FILE *out;
    char name[100] = { };

    sprintf(name, "POVRAY/Shape.ini");
    out = fopen(name, "w");

    fprintf(out, "Input_File_Name=Shape.pov\n\n");
    fprintf(out, "; these are the default values\n");
    fprintf(out, "Initial_Clock=0.000\n");
    fprintf(out, "Final_CLock=1.000\n");
    fprintf(out, "Antialias=On\n");
    fprintf(out, "Antialias_Threshold=0.05\n\n");

    fprintf(out, "Initial_Frame=0\n");
    fprintf(out, "Final_Frame=%d\n\n", itmax/i_out);

    fprintf(out, "Height=1024\n");
    fprintf(out, "Width=1280\n\n");

    fprintf(out, "Output_to_File=true\n");
    fprintf(out, "Output_File_Type=N8\n");
    fprintf(out, "Output_File_Name=Images/\n");

    fclose(out);
}

void POVRAY_output(double N_x[], double N_y[], int it, int itmax)
{
    FILE *out;
    char name[100] = { };

    if (it < 10) { sprintf(name, "POVRAY/Shape00%d.pov", it); }
    else if (it < 100) { sprintf(name, "POVRAY/Shape0%d.pov", it); }
    else { sprintf(name, "POVRAY/Shape%d.pov", it); }

    out = fopen(name, "w");

    fprintf(out, "#include \"colors.inc\"\n");
    fprintf(out, "#include \"textures.inc\"\n");
    fprintf(out, "#include \"shapes.inc\"\n\n");
    fprintf(out, "background { color White }\n\n");

    fprintf(out, "camera { orthographic\n");
    fprintf(out, "  location <0, 20, -40>\n");
    fprintf(out, "  look_at  <0, 20, 0>\n}\n");
    
    fprintf(out, "light_source { <0, 20, -50> color White shadowless\n");
    fprintf(out, "               area_light <100, 0, 0>, <0, 100, 0>, 5, 5\n");
    fprintf(out, "               adaptive 1 jitter }\n\n");

    fprintf(out, "object{\n");
    fprintf(out, "  Round_Cylinder\n");
    fprintf(out, "   (<-200,0,0>,<200,0,0>, 0.2, 0.1, 1)\n");
    fprintf(out, "  texture{ pigment{ color Red}\n");
    fprintf(out, "    finish { reflection 0.05 phong 1}\n");
    fprintf(out, "  }\n}\n");
    
    fprintf(out, "blob {\n");
    fprintf(out, "  threshold 0.99\n");
    for (int i = 0; i < N; i++)
    {
        fprintf(out, "  sphere {\n");
        fprintf(out, "           <%.3f,%.3f,0>, 0.4, 1.0\n", N_x[i], N_y[i]);
        fprintf(out, "         }\n");
    }
    fprintf(out, "   scale 1\n");
    fprintf(out, "   pigment {rgb <0,0,0>}\n");
    fprintf(out, "   finish { phong 0.8 }\n}\n");

    fclose(out);
}

void GNUplot_output_no_border(char name[], int column1, int column2)
{
    char filename[100];
    strcpy(filename, name);
    filename[strlen(filename) - 4] = '\0';
    sprintf(filename, "%s_GNU.gp", filename); 

    FILE *out;
    out = fopen(filename, "w");
    filename[strlen(filename) - 7] = '\0';

    fprintf(out, "reset\n");
    fprintf(out, "set terminal png size 800,600 lw 3\n");
    fprintf(out, "set output 'GNUplot/%s.png'\n", filename);
    fprintf(out, "set title 'Liposome shape'\n");
    fprintf(out, "unset border\n");
    fprintf(out, "unset xtics\n");
    fprintf(out, "unset ytics\n");
    fprintf(out, "unset key\n");
    fprintf(out, "plot [-25:25] [-1:35] '%s' using %d:%d with lines, 0 with lines\n", name, column1, column2);

    fclose(out);
}

void GNUplot_output_border(char name[], char out_name[], int column1, int column2)
{
    char filename[100] = { };
    sprintf(filename, "%s_GNU.gp", out_name); 

    FILE *out;
    out = fopen(filename, "w");
    filename[strlen(filename) - 7] = '\0';

    fprintf(out, "reset\n");
    fprintf(out, "set terminal png size 800,600 lw 3\n");
    fprintf(out, "set output 'GNUplot/%s.png'\n", out_name);
    fprintf(out, "set title '%s'\n", out_name);
    fprintf(out, "plot '%s' using %d:%d with points\n", name, column1, column2);

    fclose(out);
}

void Text_output(double N_x[], double N_y[], double F_x[], double F_y[], char name[])
{
    FILE *out;
    out = fopen(name, "w");
    for (int i = 0; i < N; i++)
    {
        fprintf(out, "%.6f %.6f %.6f %.6f\n", N_x[i], N_y[i], F_x[i], F_y[i]);
    }

    fclose(out);

    GNUplot_output_no_border(name, 1, 2);
}
