#include <iostream>
#include <mpi.h>


int main(int argc, char **argv) {
	// Set up MPI
	MPI_Init(NULL, NULL);
	int rank, size;
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	
	// Set up matrix-parameters
	int M = 9;
	int local_N = M/size;
	if (rank < M % size) local_N++;
	
	// Local start colomn, NB, not important!!
	int local_start_N = 0;
	for (int i = 0; i < rank; i++) {
		local_start_N += M/size;
		if (i < M % size) local_start_N++;
	}

	// Create send and recieve buffers
	double a[M][local_N];
	double b[M][local_N];
	
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
	for (int i = 0; i < M; i++)
		for (int j = 0; j < local_N; j++)
			a[i][j] = i*M + j + local_start_N; 
	
	// Print a at rank i for reference:
	if (rank == 35) {
		std::cout << "a and b from rank: " << rank << std::endl << std::endl << "a:" << std::endl;
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < local_N; j++) {
				std::cout << a[i][j] << " ";
			}
			std::cout << std::endl;
		}
	}
	
	// Let the magic happen
	MPI_Alltoallv(&a[0][0], sendcounts, sdispls, MPI_DOUBLE, &b[0][0], sendcounts, sdispls, MPI_DOUBLE, MPI_COMM_WORLD);
	
	// Transpose the "blocks" of recieved data
	int block_M, current_i = 0;
	for (int block_idx = 0; block_idx < size; block_idx++) {
		block_M = sendcounts[block_idx] / local_N;
		double temp[block_M * local_N];
		for (int i = 0; i < block_M; i++) {
			for (int j = 0; j < local_N; j++) {
				// transpose into temp buffer:
				temp[j + local_N*i] = b[current_i+i][j];
			}
		}
		for (int j = 0; j < local_N; j++) {
			for (int i = 0; i < block_M; i++) {
				// Saving temp buffer into b buffer:
				b[current_i + i][j] = temp[i + j * block_M];
			}
		}
		current_i += block_M;
	}

	// Print b to check result
	if (rank == 0) {
		std::cout << std::endl << "b, which should be transpose: " << std::endl;
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < local_N; j++) {
				std::cout << b[i][j] << " ";
			}
			std::cout << std::endl;
		}
	}

	// Finito!!
	MPI_Finalize();
	return 0;
}

