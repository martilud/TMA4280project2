#include <iostream>
#include <mpi.h>

double *mk_1D_array(size_t n, bool zero)
{
    if (zero) {
        return (double *)calloc(n, sizeof(double));
    }
    return (double *)malloc(n * sizeof(double));
}

/*
 * The allocation of the two-dimensional array used for storing matrices is done
 * in the following way for a matrix in R^(n1*n2):
 * 1. an array of pointers is allocated, one pointer for each row,
 * 2. a 'flat' array of size n1*n2 is allocated to ensure that the memory space
 *   is contigusous,
 * 3. pointers are set for each row to the address of first element.
 */

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

void printMatrix(double **mat, int N, int M){
    for (int j = 0; j < M; j++) {
		for (int i = 0; i < N; i++) {
			if (mat[i][j] < 10) std::cout << 0;	
			std::cout << mat[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

int main(int argc, char **argv) {
	// Set up MPI
	MPI_Init(NULL, NULL);
	int rank, size;
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	
	// Set up matrix-parameters
	int M = 10;
	int local_N = M/size;
	if (rank < M % size) local_N++;
	
	// Local start colomn, NB, not important!!
	int local_start_N = 0;
	for (int i = 0; i < rank; i++) {
		local_start_N += M/size;
		if (i < M % size) local_start_N++;
		
	}

	// Create send and recieve buffers
	double **b = mk_2D_array(local_N,M, 0);
	double **bt = mk_2D_array(local_N,M, 0);


	// Find send and displacement buffers
	int sendcounts[size], sdispls[size];
	for (int i = 0; i < size; i++) {
		sendcounts[i] = M/size;
		if (i < M % size) sendcounts[i]++;
		sendcounts[i] *= local_N;
		if (rank == 0) std::cout << sendcounts[i] << " ";
	}

	sdispls[0] = 0;
	for (int i = 1; i < size; i++) {
		sdispls[i] = sdispls[i-1] + sendcounts[i-1];
	}

	
	// Fill a with increasing numbers row by row
	for (int i = 0; i < local_N; i++)
		for (int j = 0; j < M; j++)
			b[i][j] = i + j*M + local_start_N; 
	
    if(rank == 0) {
        std::cout << "original matrix" << std::endl;
        printMatrix(b, local_N, M);
    }

	double *temp = mk_1D_array(local_N * M, 0);
	int block_M, current_j = 0, count = 0;	
	for (int block_idx = 0; block_idx < size; block_idx++) {
		block_M = sendcounts[block_idx] / local_N;
		for (int i = 0; i < local_N; i++) {
			for (int j = 0; j < block_M; j++) {
				temp[count++] = b[i][j + current_j];
			}
		}
		current_j += block_M;
	}


	
	double *recieve_temp = mk_1D_array(local_N * M, 0);


	// Let the magic happen
	MPI_Alltoallv(&temp[0], sendcounts, sdispls, MPI_DOUBLE, &recieve_temp[0], sendcounts, sdispls, MPI_DOUBLE, MPI_COMM_WORLD);

	count = 0;
	for (int j= 0; j < M; j++) {
		for (int i = 0; i < local_N; i++) {
			bt[i][j] = recieve_temp[count++];
		}
	}


	// Print b to check result
	if (rank == 0){
        std::cout << "transposed matrix" << std::endl; 
        printMatrix(bt,local_N,M); 
    }

	// Finito!!
	MPI_Finalize();
	return 0;
}

