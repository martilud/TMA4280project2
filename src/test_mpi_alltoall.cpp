#include <iostream>
#include <mpi.h>

int main (int argc, char *argv[])
{
	int NROW = 8;
	int NBLK = NROW/4;

	double a[NROW][NBLK];
	double b[NROW][NBLK];
	
	int rank, size;
	double r0,r1;
	

	MPI_Init (NULL, NULL);
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	
	for (int i = 0; i < NROW; i++)
		for (int j = 0; j < NBLK; j++)
			a[i][j] = i + j + rank*10; 
	

	if (rank == 0) {
		for (int i = 0; i < NROW; i ++) {
			for (int j = 0; j < NBLK; j++) { 
				if (a[i][j] < 10)
					std::cout << " ";
				std::cout << a[i][j] << " "; 
			}
			std::cout << std::endl;
		}
	}
	

  	MPI_Alltoall (&a[0][0],NBLK * NBLK, MPI_DOUBLE, &b[0][0], NBLK * NBLK, MPI_DOUBLE, MPI_COMM_WORLD);
	
	if (rank == 0) {
		std::cout << std::endl << " Finnished Alltoall" << std::endl << std::endl;
		for (int i = 0; i < NROW; i ++) {
			for (int j = 0; j < NBLK; j++) { 
				if (b[i][j] < 10)
					std::cout << " ";
				std::cout << b[i][j] << " "; 
			}
			std::cout << std::endl;
		}
	}

	MPI_Finalize();
}
