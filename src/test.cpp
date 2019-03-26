#include <mpi.h>
#include <iostream>

double **mk_2D_array(size_t n1, size_t n2, bool zero);

void transpose(double **a, double **b, int M, int N); 

int *mk_1D_array(size_t n, bool zero);


int main(int argc, char **argv) {
	MPI_Init(NULL,NULL);
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	

	int *sendcount = mk_1D_array_int(size, true), *sdispls = mk_1D_array_int(size, true), *recvcount = mk_1D_array_int(size, true), *rdispls = mk_1D_array_int(size, true);
	
	int N = 17;
	int local_M = N / size;
	for (int i = 0; i < size; i++) {
		if (i < N % size)
			sencount
			local_M++;
			
	}


	for (int i = rank; i < N; i += size) {
		local_M++;
	}

	double **a_local = mk_2D_array(local_M, N, false),  **b_local = mk_2D_array(local_M, N, false);
	
	for (int i = 0; i < local_M; i++) {
		for (int j = 0; j < N; j++) {
			a_local[i][j] = i*j*rank;
		}	
	}

	int *sendcount = mk_1D_array_int(size, true), *sdispls = mk_1D_array_int(size, true), *recvcount = mk_1D_array_int(size, true), *rdispls = mk_1D_array_int(size, true);
	for (int i = 0; i < size; i++) {
		
	}
	
//	MPI_Type_vector(local_M * N);
/*//
	// MÅL: send den transponerte:
	int *sdispls = mk_1D_array(local_M * N);

	for (int 

	MPI_Alltoallv(&a_local[0], local_M * N, sdispls
		
	//	for (int i = rank; i < N; i += size) {
	//		for (int j = 0; j < N; j++) {
	//			std::cout << a_local[j][i] << " " ;
	//		} 
	//		std::cout << std::endl;
	//	}
*/
	MPI_Finalize();
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
