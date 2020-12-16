#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

#define DELTA 0.01
#define SIGMA (10.0 * DELTA)
#define IT_MAX 10000
#define NX 400
#define NY 90
#define I1 200
#define I2 210
#define J1 50
#define XA 0.45
#define YA 0.45
#define MU 20
#define SNAP_COUNT 5

void fill_v(double **vx, double **vy, double **psi) {
    for (int i = 1; i < NX; ++i) {
        for (int j = 1; j < NY; ++j) {
            vx[i][j] = (psi[i][j+1] - psi[i][j-1]) / (2*DELTA);
            vy[i][j] = -1.0 * (psi[i+1][j] - psi[i-1][j]) / (2*DELTA);
        }
    }
    for (int i = I1; i <= I2; ++i) {
        for (int j = 0; j <= J1; ++j) {
            vx[i][j] = vy[i][j] = 0.0;
        }
    }
    for (int i = 1; i < NX; ++i) {
        vy[i][0] = vy[i][NY] = 0.0;
    }
    for (int j = 0; j <= NY; ++j) {
        vx[0][j] = vx[1][j];
        vx[NX][j] = vx[NX-1][j];
    }
}

void initialize_u(double **u0) {
    double x, y;
    for (int i = 0; i <= NX; ++i) {
        for (int j = 0; j <= NY; ++j) {
            x = i * DELTA;
            y = j * DELTA;
            u0[i][j] = 1.0 / (2.0 * M_PI * SIGMA*SIGMA) * exp(-1.0 * (pow(x-XA, 2) + pow(y-YA, 2)) / (2.0 * SIGMA*SIGMA));
        }
    }
}

int is_border(int i, int j) {
    return j <= J1 && (i >= I1 && i <= I2);
}

void picard_next(double **u0, double **u1, double **vx, double **vy, double D, double dt) {
    double u_new[NX+1][NY+1];
    for (int i = 0; i <= NX; ++i) {
        for (int j = 1; j < NY; ++j) {
            if (is_border(i, j)) {
                continue;
            }
            else if (i == 0) {
                double multiplier = 1.0 / (1.0 + (2.0 * D * dt) / (DELTA*DELTA));
                double comp1 = u0[i][j] - dt / 2.0 * vx[i][j] * ((u0[i+1][j] - u0[NX][j]) / (2.0 * DELTA) + (u1[i+1][j] - u1[NX][j]) / (2.0 * DELTA));
                double comp2 = -1.0 * dt / 2.0 * vy[i][j] * ((u0[i][j+1] - u0[i][j-1]) / (2.0 * DELTA) + (u1[i][j+1] - u1[i][j-1]) / (2.0 * DELTA));
                double comp3 = dt / 2.0 * D * (u0[i+1][j] + u0[NX][j] + u0[i][j+1] + u0[i][j-1] - 4.0 * u0[i][j]) / (DELTA*DELTA);
                double comp4 = dt / 2.0 * D * (u1[i+1][j] + u1[NX][j] + u1[i][j+1] + u1[i][j-1]) / (DELTA*DELTA);
                u_new[i][j] = multiplier * (comp1 + comp2 + comp3 + comp4);
                u1[i][j] = multiplier * (comp1 + comp2 + comp3 + comp4);
            }
            else if (i == NX) {
                double multiplier = 1.0 / (1.0 + (2.0 * D * dt) / (DELTA*DELTA));
                double comp1 = u0[i][j] - dt / 2.0 * vx[i][j] * ((u0[0][j] - u0[i-1][j]) / (2.0 * DELTA) + (u1[0][j] - u1[i-1][j]) / (2.0 * DELTA));
                double comp2 = -1.0 * dt / 2.0 * vy[i][j] * ((u0[i][j+1] - u0[i][j-1]) / (2.0 * DELTA) + (u1[i][j+1] - u1[i][j-1]) / (2.0 * DELTA));
                double comp3 = dt / 2.0 * D * (u0[0][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4.0 * u0[i][j]) / (DELTA*DELTA);
                double comp4 = dt / 2.0 * D * (u1[0][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1]) / (DELTA*DELTA);
                u_new[i][j] = multiplier * (comp1 + comp2 + comp3 + comp4);
                u1[i][j] = multiplier * (comp1 + comp2 + comp3 + comp4);
            }
            else {
                double multiplier = 1.0 / (1.0 + (2.0 * D * dt) / (DELTA*DELTA));
                double comp1 = u0[i][j] - dt / 2.0 * vx[i][j] * ((u0[i+1][j] - u0[i-1][j]) / (2.0 * DELTA) + (u1[i+1][j] - u1[i-1][j]) / (2.0 * DELTA));
                double comp2 = -1.0 * dt / 2.0 * vy[i][j] * ((u0[i][j+1] - u0[i][j-1]) / (2.0 * DELTA) + (u1[i][j+1] - u1[i][j-1]) / (2.0 * DELTA));
                double comp3 = dt / 2.0 * D * (u0[i+1][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4.0 * u0[i][j]) / (DELTA*DELTA);
                double comp4 = dt / 2.0 * D * (u1[i+1][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1]) / (DELTA*DELTA);
                u_new[i][j] = multiplier * (comp1 + comp2 + comp3 + comp4);
                u1[i][j] = multiplier * (comp1 + comp2 + comp3 + comp4);
            }
        }
    }
    for (int i = 0; i <= NX; ++i) {
        for (int j = 1; j < NY; ++j) {
            u1[i][j] = u_new[i][j];
        }
    }
}

void ad_calculate(double **u0, double **u1, double **vx, double **vy, double D, double dt, FILE* file, FILE** ss_files) {
    int l = 0;
    int it_snap_step = IT_MAX / SNAP_COUNT;

    for (int it = 0; it < IT_MAX; ++it) {
        for (int i = 0; i <= NX; ++i) {
            for (int j = 0; j <= NY; ++j) {
                u1[i][j] = u0[i][j];
            }
        }

        for (int k = 0; k < MU; ++k) {
            picard_next(u0, u1, vx, vy, D, dt);
        }

        double c = 0.0;
        double xsr = 0.0;
        for (int i = 0; i <= NX; ++i) {
            for (int j = 1; j < NY; ++j) {
                u0[i][j] = u1[i][j];
                c += u0[i][j] * DELTA*DELTA;
                xsr += i * DELTA * u0[i][j] * DELTA*DELTA;
            }
        }
        fprintf(file, "%f %f %f\n", it * dt, c, xsr);

        if (it > l * it_snap_step + it_snap_step / 2) {
            for (int i = 0; i <= NX; ++i) {
                for (int j = 0; j <= NY; ++j) {
                    fprintf(ss_files[l], "%f %f %f\n", i * DELTA, j * DELTA, u0[i][j]);
                }
                fprintf(ss_files[l], "\n");
            }
            ++l;
        }
    }
}

int main(int argc, const char* argv[]) {
    double** u0 = alloc_matrix(NX+1, NY+1);
    double** u1 = alloc_matrix(NX+1, NY+1);
    double** vx = alloc_matrix(NX+1, NY+1);
    double** vy = alloc_matrix(NX+1, NY+1);
    double** psi = alloc_matrix(NX+1, NY+1);

    FILE* psi_data = fopen("psi.dat", "r");
    for (int i = 0; i <= NX; ++i) {
        for (int j = 0; j <= NY; ++j) {
            fscanf(psi_data, "%*d %*d %lf\n", &psi[i][j]);
        }
    }
    fill_v(vx, vy, psi);

    double v_max = 0.0;
    double temp;
    for (int i = 0; i <= NX; ++i) {
        for (int j = 0; j <= NY; ++j) {
            temp = sqrt(vx[i][j] * vx[i][j] + vy[i][j] * vy[i][j]);
            if (v_max < temp) {
                v_max = temp;
            }
        }
    }
    double dt = DELTA / (4.0 * v_max);

    FILE* file_vx= fopen("vx.dat", "w");
    FILE* file_vy = fopen("vy.dat", "w");

    for (int i = 0; i <= NX; ++i) {
        for (int j = 0; j <= NY; ++j) {
            fprintf(file_vx, "%f %f %f\n", i * DELTA, j * DELTA, vx[i][j]);
            fprintf(file_vy, "%f %f %f\n", i * DELTA, j * DELTA, vy[i][j]);
        }
        fprintf(file_vx, "\n");
        fprintf(file_vy, "\n");
    }

    FILE* file_no_diff= fopen("no_diff.dat", "w");
    FILE* file_with_diff= fopen("with_diff.dat", "w");
    FILE* dist_no_diff[SNAP_COUNT];
    FILE* dist_with_diff[SNAP_COUNT];

    for (int i = 0; i < SNAP_COUNT; ++i) {
        char buff1[24];
        char buff2[24];
        sprintf(buff1, "dist_no_diff_%d.dat", i+1);
        sprintf(buff2, "dist_with_diff_%d.dat", i+1);
        dist_no_diff[i] = fopen(buff1, "w");
        dist_with_diff[i] = fopen(buff2, "w");
    }

    initialize_u(u0);
    ad_calculate(u0, u1, vx, vy, 0.0, dt, file_no_diff, dist_no_diff);
    initialize_u(u0);
    ad_calculate(u0, u1, vx, vy, 0.1, dt, file_with_diff, dist_with_diff);

    fclose(file_no_diff);
    fclose(file_with_diff);
    for (int i = 0; i < SNAP_COUNT; ++i) {
        fclose(dist_no_diff[i]);
        fclose(dist_with_diff[i]);
    }

    free_matrix(u0, NX+1);
    free_matrix(u1, NX+1);
    free_matrix(vx, NX+1);
    free_matrix(vy, NX+1);
    free_matrix(psi, NX+1);
}
