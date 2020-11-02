#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DELTA (1e-10)

typedef void (*num_method)(double, double, double, double, double*, double*);
typedef struct
{
    const double x0;
    const double v0;
    const double dt0;
    const double S;
    const int p;
    const double tmax;
    const int alpha;
} equation_parameters;

double f(double v) { return v; }
double g(double x, double v, double alpha) { return alpha * (1 - x*x) * v - x;}

void trap_calculate(double xn, double vn, double dt, double alpha, double* x_res, double* v_res)
{
    double xn1 = xn;
    double vn1 = vn;
    double a11 = 1.0;
    double a12 = -dt / 2.0;
    double a21, a22;
    double F, G;
    double dx, dv;

    do {
        F = xn1 - xn - dt / 2.0 * (f(vn) + f(vn1));
        G = vn1 - vn - dt / 2.0 * (g(xn, vn, alpha) + g(xn1, vn1, alpha));
        a21 = -dt / 2.0 * (-2.0 * alpha * xn1 * vn1 - 1.0);
        a22 = 1.0 - dt / 2.0 * alpha * (1.0 - xn1*xn1);
        dx = (-F * a22 + G * a12) / (a11 * a22 - a12 * a21);
        dv = (-G * a11 + F * a21) / (a11 * a22 - a12 * a21);
        xn1 += dx;
        vn1 += dv;
    } while (fabs(dx) >= DELTA || fabs(dv) >= DELTA);
    *x_res = xn1;
    *v_res = vn1;
}

void rk2_calculate(double xn, double vn, double dt, double alpha, double* x_res, double* v_res)
{
    double k1x = f(vn);
    double k1v = g(xn, vn, alpha);
    double k2x = f(vn + dt * k1v);
    double k2v = g(xn + dt * k1x, vn + dt * k1v, alpha);

    *x_res = xn + dt / 2.0 * (k1x + k2x);
    *v_res = vn + dt / 2.0 * (k1v + k2v);
}

void solve_with_step_correction(equation_parameters params, double tolerance, num_method f, FILE* results_file)
{
    double t = 0.0;
    double dt = params.dt0;

    double xn = params.x0;
    double vn = params.v0;
    double xn1_two_step, xn2_two_step, vn1_two_step, vn2_two_step;
    double xn2_one_step, vn2_one_step;
    do {
        f(xn, vn, dt, params.alpha, &xn1_two_step, &vn1_two_step);
        f(xn1_two_step, vn1_two_step, dt, params.alpha, &xn2_two_step, &vn2_two_step);

        f(xn, vn, 2 * dt, params.alpha, &xn2_one_step, &vn2_one_step);

        double ex = (xn2_two_step - xn2_one_step) / (pow(2, params.p) - 1);
        double ev = (vn2_two_step - vn2_one_step) / (pow(2, params.p) - 1);

        double max_val = fmax(fabs(ex), fabs(ev));
        if (max_val < tolerance)
        {
            t += 2 * dt;
            xn = xn2_two_step;
            vn = vn2_two_step;
            fprintf(results_file, "%g %g %g %g\n", t, dt, xn, vn);
        }
        dt *= pow(params.S * tolerance / max_val, 1.0 / (params.p + 1.0));
    } while(t < params.tmax);
}

int main(int argc, const char* argv[])
{
    const double x0 = 0.01;
    const double v0 = 0.0;
    const double dt0 = 1.0;
    const double S = 0.75;
    const int p = 2;
    const double tmax = 40.0;
    const int alpha = 5;
    equation_parameters params = {x0, v0, dt0, S, p, tmax, alpha};

    const double tol1 = 1e-2;
    const double tol2 = 1e-5;

    FILE* file1 = fopen("trap_tol1.dat", "w");
    FILE* file2 = fopen("trap_tol2.dat", "w");
    FILE* file3 = fopen("rk2_tol1.dat", "w");
    FILE* file4 = fopen("rk2_tol2.dat", "w");

    solve_with_step_correction(params, tol1, trap_calculate, file1);
    solve_with_step_correction(params, tol2, trap_calculate, file2);
    solve_with_step_correction(params, tol1, rk2_calculate, file3);
    solve_with_step_correction(params, tol2, rk2_calculate, file4);

    fclose(file1);
    fclose(file2);
    fclose(file3);
    fclose(file4);

    return 0;
}
