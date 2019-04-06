#include <iostream>
#include <mpi.h>


int main(int argc, char **argv) {
	// Set up MPI
	MPI_Init(NULL, NULL);
	int rank, size;
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		std::cout << "Excecuting unit test" << std::endl;
	}
	
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
	double **b = (double **)malloc(local_N * sizeof(double *));
	double **bt = (double **)malloc(local_N * sizeof(double *));
	b[0] = (double *)malloc(local_N * M * sizeof(double));
	bt[0] = (double *)malloc(local_N * M * sizeof(double));
    for (size_t i = 1; i < local_N; i++) {
        b[i] = b[i-1] + M;
        bt[i] = bt[i-1] + M;
    }


	// Find send and displacement buffers
	int sendcounts[size], sdispls[size];
	for (int i = 0; i < size; i++) {
		sendcounts[i] = M/size;
		if (i < M % size) sendcounts[i]++;
		sendcounts[i] *= local_N;
	}

	sdispls[0] = 0;
	for (int i = 1; i < size; i++) {
		sdispls[i] = sdispls[i-1] + sendcounts[i-1];
	}

	
	// Fill a with increasing numbers row by row
	for (int i = 0; i < local_N; i++)
		for (int j = 0; j < M; j++)
			b[i][j] = i + j*M + local_start_N; //fill_matrix(i, j, M, local_start_N); // i + j*M + local_start_N; 
	

	// Creating a send and recieve temporary 1D array for wrapping/unwrapping the data
	double *temp = (double *)malloc(local_N * M * sizeof(double));
	double *recieve_temp = (double *)malloc(local_N * M * sizeof(double));


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


	// Let the magic happen
	MPI_Alltoallv(&temp[0], sendcounts, sdispls, MPI_DOUBLE, &recieve_temp[0], sendcounts, sdispls, MPI_DOUBLE, MPI_COMM_WORLD);

	count = 0;
	for (int j= 0; j < M; j++) {
		for (int i = 0; i < local_N; i++) {
			bt[i][j] = recieve_temp[count++];
		}
	}
	
	int correct = 1;
	int recieve_correct;

	for (int i = 0; i < local_N; i++) {
		for (int j = 0; j < M; j++) {
			if (bt[i][j] != j + M * (i + local_start_N)){ //fill_matrix(j, i, M, local_start_N)) {
				correct = 0;
			}
		}
	}

	MPI_Reduce(&correct, &recieve_correct, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		if (recieve_correct == size) {
			std::cout << "Unit test is succesfull!!" << std::endl;
		}
		else {
			std::cout << "Unit test failed!!" << std::endl; 
		}
	}


	// Finito!!
	MPI_Finalize();
	return 0;
}

