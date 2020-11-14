#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

#define SIGMA_X(x_max) (0.1 * x_max)
#define SIGMA_Y(y_max) (0.1 * y_max)

double rho(double x, double y, double x_max, double y_max)
{
    double rho1 = exp(-pow((x - 0.35 * x_max) / (SIGMA_X(x_max)), 2)
            - pow((y - 0.5 * y_max) / (SIGMA_Y(y_max)), 2));
    double rho2 = -exp(-pow((x - 0.65 * x_max) / (SIGMA_X(x_max)), 2)
            - pow((y - 0.5 * y_max) / (SIGMA_Y(y_max)), 2));
    return rho1 + rho2;
}

double calc_sum(double** V, int max_x, int max_y, double delta, double** rho_values)
{
    double sum = 0.0;
    for (int i = 0; i < max_x; i++)
    {
        for (int j = 0; j < max_y; j++)
        {
            double temp1 = (V[i+1][j] - V[i][j]) / delta;
            double temp2 = (V[i][j+1] - V[i][j]) / delta;
            sum += (delta * delta) * (0.5 * temp1 * temp1 + 0.5 * temp2 * temp2 - rho_values[i][j] * V[i][j]);
        }
    }
    return sum;
}

void calculate_global(equation_params params, FILE* file_sum, FILE* file_result, FILE* file_sigma)
{
    double** Vs = alloc_matrix(params.nx+1, params.ny+1);
    double** Vn = alloc_matrix(params.nx+1, params.ny+1);
    double** rho_values = alloc_matrix(params.nx+1, params.ny+1);

    for (int i = 0; i < params.nx + 1; i++) {
        Vs[i][0] = Vn[i][0] = params.V1;
        Vs[i][params.ny] = Vn[i][params.ny] = params.V2;
        for (int j = 0; j < params.ny + 1; j++) {
            rho_values[i][j] = rho(i * params.delta, j * params.delta, params.x_max, params.y_max);
        }
    }

    int iter_count = 0;
    double prev_S, curr_S;
    curr_S = calc_sum(Vs, params.nx, params.ny, params.delta, rho_values);
    do {
        prev_S = curr_S;
        for (int i = 1; i < params.nx; i++)
        {
            for (int j = 1; j < params.ny; j++)
            {
                Vn[i][j] = 0.25 * (Vs[i+1][j] + Vs[i-1][j] + Vs[i][j+1] + Vs[i][j-1] + params.delta * params.delta
                        / params.eps * rho_values[i][j]);
            }
        }

        for (int j = 0; j < params.ny + 1; j++)
        {
            Vn[0][j] = Vn[1][j];
            Vn[params.nx][j] = Vn[params.nx-1][j];
        }

        for (int i = 0; i < params.nx + 1; i++)
        {
            for (int j = 0; j < params.ny + 1; j++)
            {
                Vs[i][j] = (1.0 - params.omega) * Vs[i][j] + params.omega * Vn[i][j];
            }
        }
        curr_S = calc_sum(Vs, params.nx, params.ny, params.delta, rho_values);
        fprintf(file_sum, "%d %f\n", iter_count++, curr_S);
    }
    while (fabs((curr_S - prev_S) / prev_S) > params.TOL);

    for (int i = 0; i < params.nx + 1; ++i) {
        for (int j = 0; j < params.ny + 1; ++j) {
            fprintf(file_result, "%d %d %f\n", i, j, Vs[i][j]);
        }
        fprintf(file_result, "\n");
    }

    for (int i = 1; i < params.nx; ++i) {
        for (int j = 1; j < params.ny; ++j) {
            double temp1 = (Vs[i+1][j] - 2.0*Vs[i][j] + Vs[i-1][j]) / (params.delta * params.delta);
            double temp2 = (Vs[i][j+1] - 2.0*Vs[i][j] + Vs[i][j-1]) / (params.delta * params.delta);
            fprintf(file_sigma, "%d %d %f\n", i, j, temp1 + temp2 + rho_values[i][j] / params.eps);
        }
        fprintf(file_sigma, "\n");
    }

    free_matrix(Vs, params.nx + 1);
    free_matrix(Vn, params.nx + 1);
    free_matrix(rho_values, params.nx + 1);
}

void calculate_local(equation_params params, FILE* file)
{
    double **V = alloc_matrix(params.nx + 1, params.ny + 1);
    double **rho_values = alloc_matrix(params.nx + 1, params.ny + 1);
    for (int i = 0; i < params.nx + 1; i++) {
        V[i][0] = params.V1;
        V[i][params.ny] = params.V2;
        for (int j = 0; j < params.ny + 1; j++) {
            rho_values[i][j] = rho(i * params.delta, j * params.delta, params.x_max, params.y_max);
        }
    }

    int iter_count = 0;
    double prev_S, curr_S;
    curr_S = calc_sum(V, params.nx, params.ny, params.delta, rho_values);
    do {
        for (int i = 1; i < params.nx; i++) {
            for (int j = 1; j < params.ny; j++) {
                V[i][j] = (1.0 - params.omega) * V[i][j] + params.omega / 4.0 * (V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1]
                                                                   + params.delta * params.delta / params.eps * rho_values[i][j]);
            }
        }
        for (int j = 0; j < params.ny + 1; j++) {
            V[0][j] = V[1][j];
            V[params.nx][j] = V[params.nx - 1][j];
        }
        prev_S = curr_S;
        curr_S = calc_sum(V, params.nx, params.ny, params.delta, rho_values);
        fprintf(file, "%d %f\n", ++iter_count, curr_S);
    } while (fabs((curr_S - prev_S) / prev_S) > params.TOL);

    free_matrix(V, params.nx + 1);
    free_matrix(rho_values, params.nx + 1);
}

int main(int argc, const char* argv[])
{
    const double eps = 1.0;
    const double delta = 0.1;
    const int nx = 150.0;
    const int ny = 100.0;
    const double V1 = 10.0;
    const double V2 = 0.0;
    const double x_max = delta * nx;
    const double y_max = delta * ny;
    const double TOL = 1e-8;

    double omegas_g[] = {0.6, 1.0};
    double omegas_l[] = {1.0, 1.4, 1.8, 1.9};

    // Relaksacja globalna
    FILE* f_global_sum_1 = fopen("global_sum_omega1.dat", "w");
    FILE* f_global_sum_2 = fopen("global_sum_omega2.dat", "w");
    FILE* f_global_solution_1 = fopen("global_solution_omega1.dat", "w");
    FILE* f_global_solution_2 = fopen("global_solution_omega2.dat", "w");
    FILE* f_global_sigma_1 = fopen("global_sigma_omega1.dat", "w");
    FILE* f_global_sigma_2 = fopen("global_sigma_omega2.dat", "w");

    equation_params params = {nx, ny, x_max, y_max, V1, V2, delta, omegas_g[0], eps, TOL};
    calculate_global(params, f_global_sum_1, f_global_solution_1, f_global_sigma_1);
    params.omega = omegas_g[1];
    calculate_global(params, f_global_sum_2, f_global_solution_2, f_global_sigma_2);

    fclose(f_global_sum_1);
    fclose(f_global_sum_2);

    // Relaksacja lokalna
    FILE* f_local_sum_1 = fopen("local_sum_omega1.dat", "w");
    FILE* f_local_sum_2 = fopen("local_sum_omega2.dat", "w");
    FILE* f_local_sum_3 = fopen("local_sum_omega3.dat", "w");
    FILE* f_local_sum_4 = fopen("local_sum_omega4.dat", "w");

    params.omega = omegas_l[0];
    calculate_local(params, f_local_sum_1);
    params.omega = omegas_l[1];
    calculate_local(params, f_local_sum_2);
    params.omega = omegas_l[2];
    calculate_local(params, f_local_sum_3);
    params.omega = omegas_l[3];
    calculate_local(params, f_local_sum_4);

    fclose(f_local_sum_1);
    fclose(f_local_sum_2);
    fclose(f_local_sum_3);
    fclose(f_local_sum_4);
}
