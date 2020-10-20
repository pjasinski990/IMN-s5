#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double exact(double lambda, double t)
{
    return exp(t * lambda);
}

int main(int argc, char const *argv[])
{
    // Zadanie 1sze 

    const double y0 = 0.0;
    const double fy0 = 1.0;
    const double lambda = -1.0;
    const double tmin = 0;
    const double tmax = 5.0;
    const double ts[3] = {0.01, 0.1, 1.0};

    // Euler
    FILE* f1 = fopen("1a.dat", "w");
    for (int i = 0; i < 3; i++)
    {
        double t = tmin;
        double y = fy0;
        while (t <= tmax)
        {
            fprintf(f1, "%f %f %f\n", t, y, fabs(y - exact(lambda, t)));
            t += ts[i];
            y += y*lambda*ts[i];
        }
        fprintf(f1, "\n\n");
    }
    fclose(f1);

    // jawna RK2 (trapezow)
    FILE* f2 = fopen("1b.dat", "w");
    for (int i = 0; i < 3; i++)
    {
        double t = tmin;
        double y = fy0;
        double k1, k2;
        while (t <= tmax)
        {
            k1 = lambda * y;
            k2 = lambda * (y + ts[i] * k1);
            fprintf(f2, "%f %f %f\n", t, y, fabs(y - exact(lambda, t)));
            t += ts[i];
            y += ts[i]/2 * (k1 + k2);
        }
        fprintf(f2, "\n\n");
    }
    fclose(f2);

    // jawna RK4
    FILE* f3 = fopen("1c.dat", "w");
    for (int i = 0; i < 3; i++)
    {
        double t = tmin;
        double y = fy0;
        double k1, k2, k3, k4;
        while (t <= tmax)
        {
            k1 = lambda * y;
            k2 = lambda * (y + ts[i]/2 * k1);
            k3 = lambda * (y + ts[i]/2 * k2);
            k4 = lambda * (y + ts[i] * k3);
            fprintf(f3, "%f %f %f\n", t, y, fabs(y - exact(lambda, t)));
            t += ts[i];
            y += ts[i]/6 * (k1 + 2*k2 + 2*k3 + k4);
        }
        fprintf(f3, "\n\n");
    }
    fclose(f3);


    // Zadanie 2 - RRZ 2go rzedu

    double V(double omegaV, double t) 
    {
        return 10 * sin(omegaV * t);
    }

    const double dt = 10e-4;
    const double R = 100;
    const double L = 0.1;
    const double C = 0.001;
    const double omega0 = 1 / sqrt(L*C);
    const double T0 = 2 * M_PI / omega0;
    const double tmin2 = 0.0;
    const double tmax2 = 4.0 * T0;
    double k1q, k2q, k3q, k4q, k1i, k2i, k3i, k4i;

    const double Q0 = 0;
    const double I0 = 0;

    double omegas[] = {0.5*omega0, 0.8*omega0, omega0, 1.2*omega0};

    FILE* f4 = fopen("2.dat", "w");
    for (int i = 0; i < 4; i++)
    {
        double t = tmin2;
        double Q = Q0;
        double I = I0;
        while (t <= tmax2) 
        {
            fprintf(f4, "%f %f %f\n", t, Q, I);
            k1q = I;
            k1i = V(omegas[i], t);
            k2q = I + dt/2 * k1i;
            k2i = V(omegas[i], t + dt/2) / L + 1/(L*C) * (Q + dt/2 * k1q) - R/L * (I + dt/2 * k1i);
            k3q = I + dt/2 * k2i;
            k3i = V(omegas[i], t + dt/2) / L - 1/(L*C) * (Q + dt/2 * k2q) - R/L * (I + dt/2 * k2i);
            k4q = I + dt * k3i;
            k4i = V(omegas[i], t + dt) / L - 1/(L*C) * (Q + dt * k3q) - R/L * (I + dt * k3i);

            Q += dt/6 * (k1q + 2*k2q + 2*k3q + k4q);
            I += dt/6 * (k1i + 2*k2i + 2*k3q + k4i);
            t += dt;
        }
        fprintf(f4, "\n\n");
    }
    fclose(f4);
    return 0;
}
