#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_linalg.h>

#define NX 40
#define NY 40
#define N ((NX+1) * (NY+1))
#define DELTA 1.0
#define DT 1.0
#define TA 40.0
#define TB 0.0
#define TC 30.0
#define TD 0.0
#define KB 0.1
#define KD 0.6
#define IT_MAX 2000

int get_l(int i, int j) {
    return i + j * (NX + 1);
}

int get_i(int j, int l) {
    return l - j * (NX + 1);
}

int get_j(int l) {
    return l / (NX + 1);
}

void vector_to_file(gsl_vector* vector, FILE* file) {
    for (int i = 0; i <= NX; ++i) {
        for (int j = 0; j <= NY; ++j) {
            int l = get_l(i, j);
            fprintf(file, "%lf %lf %g\n", i * DELTA, j * DELTA, fabs(gsl_vector_get(vector, l)));
        }
        fprintf(file, "\n");
    }
}

int main(int argc, const char* argv[]) {
    gsl_matrix *mA = gsl_matrix_calloc(N, N);
    gsl_matrix *mB = gsl_matrix_calloc(N, N);
    gsl_vector *c = gsl_vector_calloc(N);
    gsl_vector *d = gsl_vector_calloc(N);
    gsl_vector *T = gsl_vector_calloc(N);

    // interior
    for (int i = 1; i < NX; ++i) {
        for (int j = 1; j < NY; ++j) {
            int l = get_l(i, j);
            double temp = DT / (2.0 * DELTA*DELTA);

            gsl_matrix_set(mA, l, l-NX-1, temp);
            gsl_matrix_set(mA, l, l-1, temp);
            gsl_matrix_set(mA, l, l+1, temp);
            gsl_matrix_set(mA, l, l+NX+1, temp);
            gsl_matrix_set(mA, l, l, -4.0 * temp - 1.0);

            gsl_matrix_set(mB, l, l-NX-1, -temp);
            gsl_matrix_set(mB, l, l-1, -temp);
            gsl_matrix_set(mB, l, l+1, -temp);
            gsl_matrix_set(mB, l, l+NX+1, -temp);
            gsl_matrix_set(mB, l, l, 4.0 * temp - 1.0);
        }
    }
    // dirichlet
    for (int j = 0; j <= NY; ++j) {
        int l1 = get_l(0, j);
        int l2 = get_l(NX, j);
        gsl_matrix_set(mA, l1, l1, 1.0);
        gsl_matrix_set(mA, l2, l2, 1.0);
        gsl_matrix_set(mB, l1, l1, 1.0);
        gsl_matrix_set(mB, l2, l2, 1.0);
    }
    // vneumann
    // up
    for (int i = 1; i < NX; ++i) {
        int l = get_l(i, NY);
        gsl_matrix_set(mA, l, l-NX-1, -1.0 / (KB * DELTA));
        gsl_matrix_set(mA, l, l, 1.0 + 1.0 / (KB * DELTA));
        gsl_vector_set(c, l, TB);
        for (int j = 0; j < N; ++j) {
            gsl_matrix_set(mB, l, j, 0.0);
        }
    }
    // down
    for (int i = 1; i < NX; ++i) {
        int l = get_l(i, 0);
        gsl_matrix_set(mA, l, l+NX+1, -1.0 / (KD * DELTA));
        gsl_matrix_set(mA, l, l, 1.0 + 1.0 / (KD * DELTA));
        gsl_vector_set(c, l, TD);
        for (int j = 0; j < N; ++j) {
            gsl_matrix_set(mB, l, j, 0.0);
        }
    }

    // T0
    for (int j = 0; j <= NY; ++j) {
        int l1 = get_l(0, j);
        int l2 = get_l(NX, j);
        gsl_vector_set(T, l1, TA);
        gsl_vector_set(T, l2, TC);
    }

    int s;
    gsl_permutation *p = gsl_permutation_calloc(N);
    gsl_linalg_LU_decomp(mA, p, &s);

    gsl_vector* T_prev = gsl_vector_calloc(N);
    for (int it = 1; it <= IT_MAX; ++it) {
        gsl_vector_memcpy(d, c);
        gsl_blas_dgemv(CblasNoTrans, 1.0, mB, T, 1.0, d);   // d = B * T + c

        gsl_vector_memcpy(T_prev, T);
        gsl_linalg_LU_solve(mA, p, d, T);

        if (it == 100 || it == 200 || it == 500 || it == 1000 || it == 2000) {
            char f_name[24];
            char f_diff_name[24];
            sprintf(f_name, "T_it%d.dat", it);
            sprintf(f_diff_name, "T_diff_it%d.dat", it);

            FILE* file = fopen(f_name, "w");
            FILE* file_diff = fopen(f_diff_name, "w");

            gsl_vector_sub(T_prev, T);
            gsl_vector_scale(T_prev, 1.0 / DT);

            vector_to_file(T, file);
            vector_to_file(T_prev, file_diff);

            fclose(file);
            fclose(file_diff);
        }
    }

    gsl_matrix_free(mA);
    gsl_matrix_free(mB);
    gsl_vector_free(c);
    gsl_vector_free(d);
    gsl_vector_free(T);
    gsl_vector_free(T_prev);
    gsl_permutation_free(p);
}
