#include <mpi.h>
#include <iostream>

double **mk_2D_array(size_t n1, size_t n2, bool zero);

void transpose(double **a, double **b, int M, int N); 

int *mk_1D_array(size_t n, bool zero);

void print(double **a, int M, int N);

int main(int argc, char **argv) {
	MPI_Init(NULL,NULL);
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	

	int M = 4;
	int local_N = M / size;
	double *a = (double *) malloc(8 * sizeof(double *));
	double *b = (double *) malloc(8 * sizeof(double *));
	for (int i = 0; i < 8; i++) {
		a[i] = i * rank;
	}
	if (rank == 0) {
		for (int i = 0; i < 8; i++) {
			std::cout << a[i] << " ";
		}
		std::cout << std::endl;
	}
		

	MPI_Alltoall(a, 2, MPI_DOUBLE, b, 2, MPI_DOUBLE, MPI_COMM_WORLD);
	
	if (rank == 0) {
		for (int i = 0; i < 8; i++) {
			std::cout << b[i] << " ";
		}
		std::cout << std::endl;
	}
	MPI_Finalize();
}

void print(double **a, int M, int N) {
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			if (a[i][j] < 10) {
				std::cout << 0;
			}
			std::cout << a[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

void transpose(double **a, double **b, int M, int N) {
	
}



double **mk_2D_array(size_t n1, size_t n2, bool zero)
{
    // 1
    double **ret = (double **)malloc(n1 * sizeof(double *));

    // 2
    if (zero) {
        ret[0] = (double *)calloc(n1 * n2, sizeof(double));
    }
    else {
        ret[0] = (double *)malloc(n1 * n2 * sizeof(double));
    }
    
    // 3
    for (size_t i = 1; i < n1; i++) {
        ret[i] = ret[i-1] + n2;
    }
    return ret;
}

int *mk_1D_array(size_t n, bool zero)
{
    if (zero) {
        return (int *)calloc(n, sizeof(int));
    }
    return (int *)malloc(n * sizeof(int));
}
