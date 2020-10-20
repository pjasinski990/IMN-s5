#include <stdio.h>
#include <stdlib.h>

int main(int argc, char const *argv[])
{
    const double y0 = 0.0;
    const double fy0 = 1.0;
    const double lambda = -1.0;
    const double tmin = 0;
    const double tmax = 5.0;
    const double ts[3] = {0.01, 0.1, 1.0};

    FILE* f1 = fopen("a.dat", "w");
    for (int i = 0; i < 3; i++)
    {
        double t = tmin;
        double y = fy0;
        while (t <= tmax)
        {
            fprintf(f1, "%f %f\n", t, y);
            t += ts[i];
            y += y*lambda*ts[i];
        }
        fprintf(f1, "\n\n");
    }
    fclose(f1);

    return 0;
}
