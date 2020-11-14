typedef struct
{
    int nx;
    int ny;
    double x_max;
    double y_max;
    double V1;
    double V2;
    double delta;
    double omega;
    double eps;
    double TOL;
} equation_params;

double** alloc_matrix(int x, int y)
{
    double** result = calloc(x, sizeof(double*));
    for (int i = 0; i < x; i++)
    {
        result[i] = calloc(y, sizeof(double));
    }
    return result;
}

void free_matrix(double** matrix, int x)
{
    for (int i = 0; i < x; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}
