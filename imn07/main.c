#include <stdio.h>
#include "boundary_conditions.h"
#include "constants.h"
#include "util.h"

int is_on_edge(int i, int j) {
    if (j == 0 || j == NY) {
        return 1;
    }
    else if (i <= I1 && j <= J1) {
        return 1;
    }
    return 0;
}

double calc_psi(int i, int j, double** psi, double** zeta) {
    return 0.25 * (psi[i + 1][j] + psi[i - 1][j] + psi[i][j + 1] + psi[i][j - 1] - DELTA * DELTA * zeta[i][j]);
}


double calc_zeta(int i, int j, double** psi, double** zeta, int omega) {
    double temp1 = 0.25 * (zeta[i+1][j] + zeta[i-1][j] + zeta[i][j+1] + zeta[i][j-1]);
    if (omega) {
        double temp2 = (psi[i][j+1] - psi[i][j-1]) * (zeta[i+1][j] - zeta[i-1][j]) - (psi[i+1][j] - psi[i-1][j]) * (zeta[i][j+1] - zeta[i][j-1]);
        return temp1 - omega * RHO / (16.0 * MU) * temp2;
    }
    else {return temp1;}
}

double calc_gamma(double** psi, double** zeta) {
    double res = 0.0;
    int j2 = J1 + 2;
    for (int i = 1; i < NX; ++i) {
        res += psi[i+1][j2] + psi[i-1][j2] + psi[i][j2+1] + psi[i][j2-1] - 4.0 * psi[i][j2] - DELTA*DELTA * zeta[i][j2];
    }
    return res;
}

void ns_relax(double** psi, double** zeta, double q_in) {
    int omega = 0;
    for (int it = 1; it <= IT_MAX; ++it) {
        if (it >= 2000) {omega = 1;}

        for (int i = 1; i < NX; ++i) {
            for (int j = 1; j < NY; ++j) {
                if (!is_on_edge(i, j)) {
                    psi[i][j] = calc_psi(i, j, psi, zeta);
                    zeta[i][j] = calc_zeta(i, j, psi, zeta, omega);
                }
            }
        }
        update_zeta_bounds(psi, zeta, q_in);
        printf("gamma is: %f\n", calc_gamma(psi, zeta));
    }
}

int main(int argc, const char* argv[]) {
    double** psi = alloc_matrix(NX+1, NY+1);
    double** zeta = alloc_matrix(NX+1, NY+1);
    double q_in = -4000.0;

    FILE* file_psi = fopen("psi.dat", "w");
    FILE* file_zeta = fopen("zeta.dat", "w");

    ns_relax(psi, zeta, q_in);
    matrix_to_file(psi, NX+1, NY+1, file_psi);
    matrix_to_file(zeta, NX+1, NY+1, file_zeta);

    fclose(file_psi);
    fclose(file_zeta);
    free_matrix(psi, NX+1);
    free_matrix(zeta, NX+1);
}
