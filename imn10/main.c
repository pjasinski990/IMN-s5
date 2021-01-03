#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NT 1000
#define NX 150
#define DELTA 0.1
#define DT 0.05
#define XA 7.5
#define XF 2.5
#define SIGMA 0.5
#define T_MAX (DT * NT)

double aF(double x, double t) {
    double delta = fabs(x - XF) < 10e-8 ? 1.0 : 0.0;
    return cos(50.0 * t / T_MAX) * delta;
}

void update_a(double* a, double* u, double* u0, double alpha, double beta, double t) {
    for (int i = 1; i < NX; ++i) {
        double x =  DELTA * i;
        a[i] = (u[i+1] - 2.0 * u[i]  + u[i-1]) / (DELTA * DELTA) - beta * (u[i] - u0[i]) / DT + alpha * aF(x, t);
    }
}

void solve(double alpha, double beta) {
    double* u0 = calloc(NX+1, sizeof(double));
    double* u = calloc(NX+1, sizeof(double));
    double* v = calloc(NX+1, sizeof(double));
    double* vp = calloc(NX+1, sizeof(double));
    double* a = calloc(NX+1, sizeof(double));

    char fname_v[20];
    char fname_E[20];
    sprintf(fname_v, "v_%.1f_%.1f.dat", alpha, beta);
    sprintf(fname_E, "E_%.1f_%.1f.dat", alpha, beta);
    FILE* file_v = fopen(fname_v, "w");
    FILE* file_E = fopen(fname_E, "w");

    for (int i = 0; i <= NX; ++i) {
        double x = i * DELTA;
        double u_val = exp(-1.0 * pow(x - XA, 2) / (2.0 * SIGMA*SIGMA));
        u[i] = u_val;
        u0[i] = u_val;
    }

    update_a(a, u, u0, alpha, beta, 0.0);

    for (int n = 0; n <= NT; ++n) {
        double t = n * DT;
        for (int i = 0; i <= NX; ++i) {
            vp[i] += DT / 2.0 * a[i];
            u0[i] = u[i];
            u[i] += DT * vp[i];
        }
        update_a(a, u, u0, alpha, beta, t);
        for (int i = 1; i < NX; ++i) {
            v[i] = vp[i] + DT / 2.0 * a[i];
        }

        double E = DELTA / 4.0 * (pow((u[1] - u[0]) / DELTA, 2) + pow((u[NX] - u[NX-1]) / DELTA, 2));
        double temp = 0.0;
        for (int i = 1; i < NX; ++i) {
            temp += v[i]*v[i] + pow((u[i+1] - u[i-1]) / (2.0 * DELTA), 2);
        }
        E += DELTA / 2.0 * temp;

        for (int i = 0; i <= NX; ++i) {
            fprintf(file_v, "%g %g %g\n", t, DELTA * i, v[i]);
            fprintf(file_E, "%g %g\n", t, E);
        }
        fprintf(file_v, "\n");
    }
    fclose(file_v);
    fclose(file_E);

    free(u0);
    free(u);
    free(v);
    free(vp);
    free(a);
}

int main(int argc, const char* argv[]) {
    solve(0.0, 0.0);
    solve(0.0, 0.1);
    solve(0.0, 1.0);
    solve(1.0, 1.0);
}