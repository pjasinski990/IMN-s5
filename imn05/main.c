#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "util.h"

#define NX 128
#define NY 128

double V_B1(double y, double y_max);
double V_B2(double x, double x_max);
double V_B3(double y, double y_max);
double V_B4(double x, double x_max);
double calc_V(int i, int j, int k, double** V_mat);
double calc_sum(double** V_mat, double delta, int k);
void fill_mesh(double** V_mat, double delta, int k, double tol);
void thicken_mesh(double** V_mat, double delta, int old_k);
void write_to_file(double** V_mat, double delta, int k, FILE* file);

int main(int argc, const char* argv[])
{
    const double delta = 0.2;
    const double x_max = delta * NX;
    const double y_max = delta * NY;
    const double tol = 1e-8;
    const int k_max = 16;

    double** V_mat = alloc_matrix(NX + 1, NY + 1);
    for (int i = 0; i < NX + 1; ++i) {
        V_mat[i][NY] = V_B2(delta * i, y_max);
        V_mat[i][0] = V_B4(delta * i, y_max);
    }

    for (int i = 0; i < NY + 1; ++i) {
        V_mat[0][i] = V_B1(delta * i, x_max);
        V_mat[NX][i] = V_B3(delta * i, x_max);
    }

    FILE* fk16 = fopen("k16.dat", "w");
    FILE* fk8 = fopen("k8.dat", "w");
    FILE* fk4 = fopen("k4.dat", "w");
    FILE* fk2 = fopen("k2.dat", "w");
    FILE* fk1 = fopen("k1.dat", "w");
    FILE* f_sum = fopen("sum.dat", "w");

    FILE* files[] = {fk16, fk8, fk4, fk2, fk1};

    int k = k_max;
    int i = 0;
    while (k > 0) {
        fill_mesh(V_mat, delta, k, tol);
        write_to_file(V_mat, delta, k, files[i++]);
        thicken_mesh(V_mat, delta, k);
        k /= 2;
    }

    fclose(f_sum);
    fclose(fk1);
    fclose(fk2);
    fclose(fk4);
    fclose(fk8);
    fclose(fk16);
    free_matrix(V_mat, NX + 1);
}


double V_B1(double y, double y_max)
{
    return sin(M_PI * y / y_max);
}

double V_B2(double x, double x_max)
{
    return -1.0 * sin(2.0 * M_PI * x / x_max);
}

double V_B3(double y, double y_max)
{
    return sin(M_PI * y / y_max);
}

double V_B4(double x, double x_max)
{
    return sin(2.0 * M_PI * x / x_max);
}

double calc_V(int i, int j, int k, double** V_mat)
{
    return 0.25 * (V_mat[i+k][j] + V_mat[i-k][j] + V_mat[i][j+k] + V_mat[i][j-k]);
}

double calc_sum(double** V_mat, double delta, int k)
{
    double sum = 0.0;
    for (int i = 0; i <= NX - k; i+=k) {
        for (int j = 0; j <= NY - k; j+=k) {
            double temp1 = (V_mat[i+k][j] - V_mat[i][j]) / (2.0 * k * delta) + (V_mat[i+k][j+k] - V_mat[i][j+k]) / (2.0 * k * delta);
            double temp2 = (V_mat[i][j+k] - V_mat[i][j]) / (2.0 * k * delta) + (V_mat[i+k][j+k] - V_mat[i+k][j]) / (2.0 * k * delta);
            sum += (k * k * delta * delta) / 2.0 * (temp1 * temp1 + temp2 * temp2);
        }
    }
    return sum;
}

void fill_mesh(double** V_mat, double delta, int k, double tol)
{
    static int iter_count = 0;
    double prev_S, curr_S;
    curr_S = calc_sum(V_mat, delta, k);

    do {
        prev_S = curr_S;
        for (int i = k; i <= NX - k; i += k) {
            for (int j = k; j <= NY - k; j += k) {
                V_mat[i][j] = calc_V(i, j, k, V_mat);
            }
        }
        curr_S = calc_sum(V_mat, delta, k);
        iter_count++;
    } while(fabs((curr_S - prev_S) / prev_S) >= tol);
}

void thicken_mesh(double** V_mat, double delta, int old_k)
{
    int new_k = old_k / 2;
    for (int i = 0; i < NX; i += old_k) {
        for (int j = 0; j < NY; j += old_k) {
            V_mat[i][j + new_k] = 0.5 * (V_mat[i][j] + V_mat[i][j + old_k]);
            V_mat[i + new_k][j + new_k] = 0.25 * (V_mat[i][j] + V_mat[i+old_k][j] + V_mat[i][j+old_k] + V_mat[i+old_k][j+old_k]);
            V_mat[i + new_k][j] = 0.5 * (V_mat[i][j] + V_mat[i + old_k][j]);
        }
    }
}

void write_to_file(double** V_mat, double delta, int k, FILE* file)
{
    for (int i = 0; i < NX + 1; i += k) {
        for (int j = 0; j < NY + 1; j += k) {
            fprintf(file, "%f %f %f\n", i * delta, j * delta, V_mat[i][j]);
        }
    }
}