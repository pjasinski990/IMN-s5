#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ALPHA (beta * N - gamma)

#define BUTCH_A11 (0.25)
#define BUTCH_A12 (0.25 - sqrt(3) / 6)
#define BUTCH_A21 (0.25 + sqrt(3) / 6)
#define BUTCH_A22 (0.25)
#define BUTCH_B1 (0.5)
#define BUTCH_B2 (0.5)

double f(double u, double beta, int N, double gamma)
{
    return (beta*N - gamma) * u - beta * u*u;
}

double m11(double dt, double alpha, double beta, double U1)
{
    return 1 - dt * BUTCH_A11 * (alpha - 2 * beta * U1);
}

double m12(double dt, double alpha, double beta, double U2)
{
    return -dt * BUTCH_A12 * (alpha - 2 * beta * U2);
}

double m21(double dt, double alpha, double beta, double U1)
{
    return -dt * BUTCH_A21 * (alpha - 2 * beta * U1);
}

double m22(double dt, double alpha, double beta, double U2)
{
    return 1 - dt * BUTCH_A22 * (alpha - 2 * beta * U2);
}

double F1(double U1, double U2, double dt, double un, double alpha, double beta)
{
    return U1 - un - dt * (BUTCH_A11 * (alpha * U1 - beta * U1*U1) + BUTCH_A12 * (alpha * U2 - beta * U2*U2));
}

double F2(double U1, double U2, double dt, double un, double alpha, double beta)
{
    return U2 - un - dt * (BUTCH_A21 * (alpha * U1 - beta * U1*U1) + BUTCH_A22 * (alpha * U2 - beta * U2*U2));
}

int main(int argc, char const *argv[])
{
    const double beta = 1e-3;
    const int N = 500;
    const double tmax = 100.0;
    const double dt = 0.1;
    const double gamma = 0.1;
    const double u0 = 1.0;

    // Czesc pierwsza - metoda trapezow
    const double TOL = 1e-6;
    const int mu_max = 20;

    // Metoda Picarda
    FILE* file_pic = fopen("picard.dat", "w");
    double t = 0.0;
    double un = u0;
    while (t <= tmax) 
    {
        fprintf(file_pic, "%f\t%f\t%f\n", t, un, N - un);
        int iter_count = 0;
        double last_approx = un;
        double curr_approx = 0.0;
        while (fabs(last_approx - curr_approx) >= TOL && iter_count++ <= mu_max) 
        {
            double temp = curr_approx;
            curr_approx = un + dt / 2 * 
                    ((ALPHA * un - beta * un*un) + (ALPHA * last_approx - beta * last_approx*last_approx));
            last_approx = temp;
        }
        un = curr_approx;
        curr_approx = 0.0;
        t += dt;
    }

    // Iteracja Newtona
    FILE* file_newt = fopen("newton.dat", "w");
    t = 0.0;
    un = u0;
    while (t <= tmax)
    {
        fprintf(file_newt, "%f\t%f\t%f\n", t, un, N - un);
        int iter_count = 0;
        double last_approx = un;
        double curr_approx = 0.0;
        while (fabs(last_approx - curr_approx) >= TOL && iter_count++ <= mu_max) 
        {
            double temp = curr_approx;
            curr_approx = last_approx - (last_approx - un - dt/2 
                    * ((ALPHA * un - beta * un*un) + (ALPHA * last_approx - beta * last_approx*last_approx)))
                    / (1 - dt/2 * (ALPHA - 2 * beta * last_approx));
            last_approx = temp;
        }
        un = curr_approx;
        curr_approx = 0.0;
        t += dt;
    }



    // Czesc druga - niejawna metoda RK2
    t = 0.0;
    un = u0;
    FILE* file_rk2 = fopen("rk2.dat", "w");

    while (t <= tmax)
    {   
        fprintf(file_rk2, "%f %f %f\n", t, un, N - un);
        double u1_prev = un;
        double u2_prev = un;
        double u1_curr = 0.0;
        double u2_curr = 0.0;
        int iter_count = 0; 

        while ((fabs(u1_curr - u1_prev) >= TOL || fabs(u2_curr - u2_prev) >= TOL) && iter_count++ <= mu_max)
        {               
            if (fabs(u1_curr - u1_prev) >= TOL)
            {
                double dU1 = (F2(u1_prev, u2_prev, dt, un, ALPHA, beta) * m12(dt, ALPHA, beta, u2_prev) 
                        - F1(u1_prev, u2_prev, dt, un, ALPHA, beta) * m22(dt, ALPHA, beta, u2_prev)) 
                        / (m11(dt, ALPHA, beta, u1_prev) * m22(dt, ALPHA, beta, u2_prev) 
                        - m12(dt, ALPHA, beta, u2_prev) * m21(dt, ALPHA, beta, u1_prev));
                double temp = u1_curr;
                u1_curr = u1_prev + dU1;
                u1_prev = temp;
            }
            if (fabs(u2_curr - u2_prev) >= TOL)
            {
                double dU2 = (F1(u1_prev, u2_prev, dt, un, ALPHA, beta) * m21(dt, ALPHA, beta, u1_prev) 
                        - F2(u1_prev, u2_prev, dt, un, ALPHA, beta) * m11(dt, ALPHA, beta, u1_prev)) 
                        / (m11(dt, ALPHA, beta, u1_prev) * m22(dt, ALPHA, beta, u2_prev) 
                        - m12(dt, ALPHA, beta, u2_prev) * m21(dt, ALPHA, beta, u1_prev));
                double temp = u2_curr;
                u2_curr = u2_prev + dU2;
                u2_prev = temp;
            }
        }

        iter_count = 0;
        un += dt * (BUTCH_B1 * f(u1_curr, beta, N, gamma) + BUTCH_B2 * f(u2_curr, beta, N, gamma));
        u1_curr = 0.0;
        u2_curr = 0.0;
        t += dt;
    }
    
    fclose(file_pic);
    fclose(file_newt);
    fclose(file_rk2);
    return 0;
}
