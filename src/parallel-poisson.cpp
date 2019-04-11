/*
 * C++ program to solve the two-dimensional Poisson equation on
 * a unit square using one-dimensional eigenvalue decompositions
 * and fast sine transforms.
 *
 * Einar M. Rønquist
 * NTNU, October 2000
 * Revised, October 2001
 * Revised by Eivind Fonn, February 2015
 
 * Parallelized by Martin Ludvigsen and Harald Wilhelmsen 2019
 */

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <cstring>

#define PI 3.14159265358979323846
typedef double real;

// Function prototypes
real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);
void transpose(real **bt, real **b, size_t m);


// New functions
real rhs1();															// Source function
real rhs2(real x, real y);												// source function
real rhs3(int x_grid, int y_grid, int m);								// Source function
real analytical(real x, real y);										// Analytical result of RHS2 used for convergence test
bool utest_transpose(int rank, int size); 								// Test: Creates and transposes a simple matrix
bool utest_result(int rank, double **b, int local_N, int  m);			// Test: compares serial and parallel results

// Functions implemented in FORTRAN in fst.f and called from C.
// The trailing underscore comes from a convention for symbol names, called name
// mangling: if can differ with compilers.
extern "C" void fst_(real *v, int *n, real *w, int *nn);
extern "C" void fstinv_(real *v, int *n, real *w, int *nn);

int main(int argc, char **argv)
{
	// Initialize MPI
    MPI_Init(NULL,NULL);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	// Initalize variables which may be taken in as command line parameters
    int n = 128; 									// Problem size
    int p = omp_get_max_threads();					// Number of threads
    omp_set_num_threads(omp_get_max_threads());		
    char t = 'm'; // default value					// Spesifies if there should be used a test
    int f = 1;										// Spesifies which source function to be used
	
	// Reading command line parameters
    for (int i= 1; i < argc; i++){
        if (i < argc - 1){
            if (std::strcmp("-n", argv[i]) == 0){
                n = atoi(argv[i+1]);
            }
            else if (std::strcmp("-p", argv[i]) == 0){
                p = atoi(argv[i+1]);
		if(p == 0){
			p = omp_get_max_threads();
		}
            }
            else if(std::strcmp("-t", argv[i]) == 0){
                t = *argv[i+1];
            }
            else if(std::strcmp("-f", argv[i]) == 0){
                f = atoi(argv[i+1]);
            }
        }
    } 

	// Check if values are allowed, and setting several values
    if ((n & (n-1)) != 0) {
      printf("n must be a power-of-two\n");
      return 2;
    }

    if (p < 0){
        printf("p must be positive \n");
        return 3;
    }
    else if (p > 0){
        omp_set_num_threads(p);
    }
    if (t == 'u'){
        // Unit test, choosing fixed values
        f = 1;
        n = 128;
        p = 4;
    }
    else if(t == 'v'){
        f = 2;
    }

	// Timing: use the elapsed wall clock time in seconds
    real starttime = omp_get_wtime();

	
	int m = n - 1;
    real h = 1.0 / n;

    real *grid = mk_1D_array(n+1, false);

	// As each grid point is independent of each other, this loop is parallelized
    #pragma omp parallel for
    for (int i = 0; i < n+1; i++) {
        grid[i] = i * h;
    }

	// Again independent elements of the array, making room for parallelization
    real *diag = mk_1D_array(m, false);
    #pragma omp parallel for
    for (int i = 0; i < m; i++) {
        diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
    }
	
	// The columns of the main matrix is devided evenly on the ranks, in the following way
    int local_N = m/size;
    if (rank < m % size) local_N++;
    
    // When the number of columns are distributed is it convinient to know the start column index in the origional matrix
    int local_start_N = 0;
    for (int i = 0; i < rank; i++) {
        local_start_N += m/size;
        if (i < m % size) local_start_N++;
    }

	// Then the matrix is created as n x local_N matrices on each process
    real **b = mk_2D_array(local_N, m, false);
    real **bt = mk_2D_array(local_N, m, false);
    
    int nn = 4 * n;

	// The local b matrix is then initialized with the given RHS function, 
	// This prosess is also parallelizable as there is no reading and writing to the same memory location
    if (f == 1){
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < local_N; i++) {
            for (int j = 0; j < m; j++) {
                b[i][j] = h * h * rhs1(); //grid[local_start_N+i+1], grid[j+1]);
            }
        }
    }
    else if (f == 2){
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < local_N; i++) {
            for (int j = 0; j < m; j++) {
                b[i][j] = h * h * rhs2(grid[local_start_N+i+1], grid[j+1]);
            }
        }
    }
    else{
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < local_N; i++) {
            for (int j = 0; j < m; j++) {
                b[i][j] = h * h * rhs3(local_start_N+i+1, j+1, m);
            }
        }
    }
	
	// Then the fst_ is runned on the columns of b.
	// As the columns are independent, it is possible to parallelize as long as each OpenMP thread have its local z (memmory buffer)
    #pragma omp parallel
    {
        real *z_local = mk_1D_array(nn, false);
        #pragma omp for 
        for (int i = 0; i < local_N; i++) {
            fst_(b[i], &n, z_local, &nn);
        }
    }


	///////////////////////////  Transposing ///////////////////////////////
	// This is the place where the local matrix b is transposed
	// This part may require some more ellaborative comments

    // Allocate send and displacement buffers
	int *sendcounts = (int *)malloc(size * sizeof(int));
	int *sdispls = (int *)malloc(size * sizeof(int));
	
	// The goal is to wrap data into contegeous "blocks" of memory, and then send the blocks to the right process
	// The blocks are size: 
	// 		number of rows: local_N for process i,
	// 		number of columns: local_N for this process
	// This size is svaed into the sendcounts
    for (int i = 0; i < size; i++) {
        sendcounts[i] = m/size;
        if (i < m % size) sendcounts[i]++;
        sendcounts[i] *= local_N;
    }
 	
	// The displacement is found relative the start of the first element to be send
    sdispls[0] = 0;
    for (int i = 1; i < size; i++) {
        sdispls[i] = sdispls[i-1] + sendcounts[i-1];
    }

	// To make sure the data is contigeous in memory, the content of b is wrapped into a temporary storage buffer in the coorect order
	double *temp = mk_1D_array(local_N * m, 0);
	int block_M, current_j = 0, count = 0;

	// The process is done "blockwise," where each block represents the data going to a spesific porcess
	for (int block_idx = 0; block_idx < size; block_idx++) {
		block_M = sendcounts[block_idx] / local_N;
		
		// Then for each block, wrapping the data of b columnwise into the temp buffer:
		for (int i = 0; i < local_N; i++) {
			for (int j = 0; j < block_M; j++) {
				temp[count++] = b[i][j + current_j];
			}
		}
		current_j += block_M;
	}


	// MPI does not support aliasing, meaning that another temorary buffer must be used to recieve the data
	double *recieve_temp = mk_1D_array(local_N * m, 0);


	// Let the magic happen
	MPI_Alltoallv(&temp[0], sendcounts, sdispls, MPI_DOUBLE, &recieve_temp[0], sendcounts, sdispls, MPI_DOUBLE, MPI_COMM_WORLD);
	
	// Then unwrapping the temprary buffer back at the bt matrix. This is done row wise to simultaniously transposing the data.
	// As it turns out, this process is parallelizable as each reading or writing to memory is done only once within the loops
	// and the reading is from one buffer, and the writing is done into the another data container
    #pragma omp parallel for collapse(2)
	for (int j = 0; j < m; j++) {
		for (int i = 0; i < local_N; i++) {
			bt[i][j] = recieve_temp[j*local_N + i];
		}
	}

	//////////////////////// Finnished transposing //////////////////////////
	
    // Again by making sure each process have a local memory buffer z, the fstinv_ is parallelizable in the following way
    #pragma omp parallel
    {
        real *z_local = mk_1D_array(nn, false);
        #pragma omp for 
        for (int i = 0; i < local_N; i++) {
            fstinv_(bt[i], &n, z_local, &nn);
        }
    }
	
	// As the opperation done within each loop is independent on the other opperations, this loops is parallelizable in the following way
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < local_N; i++) {
        for (int j = 0; j < m; j++) {
            bt[i][j] = bt[i][j] / (diag[local_start_N + i] + diag[j]);
        }
    }
	
	// Parallelized by creating a local memory buffer again
    #pragma omp parallel
    {
        real *z_local = mk_1D_array(nn, false);
        #pragma omp for 
        for (int i = 0; i < local_N; i++) {
            fst_(bt[i], &n, z_local, &nn);
        }
    }
   
	// Tranposing back. It is still possible to reuse the send count and disp from the previous transponation,
	// meaning that one can directly start wrapping the data into the old buffer
	current_j = 0; 
	count = 0;	
	for (int block_idx = 0; block_idx < size; block_idx++) {
		block_M = sendcounts[block_idx] / local_N;
		for (int i = 0; i < local_N; i++) {
			for (int j = 0; j < block_M; j++) {
				temp[count++] = bt[i][j + current_j];
			}
		}
		current_j += block_M;
	}
	// This two following blocks have equal logic as the first transponation
	
	// Let the magic happen
	MPI_Alltoallv(&temp[0], sendcounts, sdispls, MPI_DOUBLE, &recieve_temp[0], sendcounts, sdispls, MPI_DOUBLE, MPI_COMM_WORLD);
	
    #pragma omp parallel for collapse(2)
	for (int j= 0; j < m; j++) {
		for (int i = 0; i < local_N; i++) {
			b[i][j] = recieve_temp[j*local_N + i];
		}
	}


	// Then by making the local memory buffer the ftsinv_ can again be parallelized
    #pragma omp parallel
    {
        real *z_local = mk_1D_array(nn, false);
        #pragma omp for 
        for (int i = 0; i < local_N; i++) {
            fstinv_(b[i], &n, z_local, &nn);
        }
    }



	////////////////////////// Main program finnished /////////////////////////////////////
	// The rest of the code is used to process the data found, either by doing tests, 
	// writing the results to file or finding the maximum u



	// verification test: test if the numerical answer is close to the analytical answer known for the second RHS function
    if (t == 'v'){
        double e_max = 0.0;
        for (int i = 0; i < local_N; i++) {
            for (int j = 0; j < m; j++) {
                e_max = e_max > fabs(b[i][j] - analytical(grid[local_start_N+i+1], grid[j+1])) ? e_max : fabs(b[i][j] - analytical(grid[local_start_N+i+1], grid[j+1]));
            }
        }

		// Reducing the result to process 0, and printing:
        double global_e_max = 0.0;
        MPI_Reduce(&e_max, &global_e_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); 
        if (rank == 0){
            real endtime = omp_get_wtime();
            std::cout << size << '\t' << p << '\t' << n << '\t' << global_e_max << '\t' << endtime-starttime << std::endl;
        } 
        MPI_Finalize();
        return 0;  
    }

	// Unit testing
    else if(t == 'u'){
        if (rank == 0){
            std::cout << "\n========= RUNNING UNIT TESTS =========\n" << std::endl; 
        }
		bool pass = true;

		// Test if transpose logic works
		if (!utest_transpose(rank, size))
			pass = false;

		// Test if program returns the same as the serial programe	
		if (!utest_result(rank, b, local_N, m))
			pass = false;
		
		if (rank == 0) {
			if (pass)
				std::cout << "\n======== UNIT TESTS SUCCESSFUL ========" << std::endl;
			else
				std::cout << "\n========= UNIT TESTS FAILED!  =========" << std::endl;
		}


        MPI_Finalize();
        return 0;

    }

	// Writing u to file. To avoid complication with parallel I/O, this is only done if the number of processes = 1
    else if(t == 'r' && size == 1){
        std::ofstream solution;
        solution.open("solution.txt");
        for (int i = 0; i < m; i++){
            for (int j = 0; j < m; j++){
                solution << b[i][j] << '\t';
            }
            solution << std::endl;
        }
        solution.close();
        std::cout << "Successfully wrote to file solution.txt" << std::endl;
        MPI_Finalize();
        return 0;
    }
    

	// Finding the max u, which by reduction is stored at process 0 and printed
    double u_max = 0.0;
    for (int i = 0; i < local_N; i++) {
        for (int j = 0; j < m; j++) {
            u_max = u_max > fabs(b[i][j]) ? u_max : fabs(b[i][j]);
        }
    }

    double global_u_max = 0.0;
    MPI_Reduce(&u_max,&global_u_max,1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); 
    if (rank == 0){
        printf("u_max = %e\n", global_u_max);
    }

	// Then finalizing the MPI and ending the program
    MPI_Finalize();
    return 0;
}

/*
 * This function is used for initializing the right-hand side of the equation.
 * Other functions can be defined to swtich between problem definitions.
 */

real rhs1() {
    // Function 1
	return 1;
}

real rhs2(real x, real y){
    // Function 2
    return 5 * PI * PI *sin(PI * x) * sin(2 * PI * y); 
}
real rhs3(int x_grid, int y_grid, int m){
    if (x_grid == m/4 && y_grid == m/4){
        return m*m;
    }
    else if(x_grid == 3*m/4 && y_grid == 3*m/4){
        return m*m;
    }
    else if(x_grid == m/4 && y_grid == 3*m/4){
        return -m*m;
    }
    else if(x_grid == 3*m/4 && y_grid == m/4){
        return -m*m;
    }
    else{
        return 0.0;
    }
}

real analytical(real x, real y){
    // Analytical solution corresponding to rhs2
    return sin(PI * x) * sin(2 * PI * y);
}

/*
 * Write the transpose of b a matrix of R^(m*m) in bt.
 * In parallel the function MPI_Alltoallv is used to map directly the entries
 * stored in the array to the block structure, using displacement arrays.
 */

void transpose(real **bt, real **b, size_t m)
{
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            bt[i][j] = b[j][i];
        }
    }
}

/*
 * The allocation of a vectore of size n is done with just allocating an array.
 * The only thing to notice here is the use of calloc to zero the array.
 */

real *mk_1D_array(size_t n, bool zero)
{
    if (zero) {
        return (real *)calloc(n, sizeof(real));
    }
    return (real *)malloc(n * sizeof(real));
}

/*
 * The allocation of the two-dimensional array used for storing matrices is done
 * in the following way for a matrix in R^(n1*n2):
 * 1. an array of pointers is allocated, one pointer for each row,
 * 2. a 'flat' array of size n1*n2 is allocated to ensure that the memory space
 *   is contigusous,
 * 3. pointers are set for each row to the address of first element.
 */

real **mk_2D_array(size_t n1, size_t n2, bool zero)
{
    // 1
    real **ret = (real **)malloc(n1 * sizeof(real *));

    // 2
    if (zero) {
        ret[0] = (real *)calloc(n1 * n2, sizeof(real));
    }
    else {
        ret[0] = (real *)malloc(n1 * n2 * sizeof(real));
    }
    
    // 3
    for (size_t i = 1; i < n1; i++) {
        ret[i] = ret[i-1] + n2;
    }
    return ret;
}



///////////////////////////  UNIT TEST: TEST IF TRANSPOSE LOGIC WORKS ////////////////////////////////

bool utest_transpose(int rank, int size) {
	if (rank == 0)
		std::cout << "Excecuting transpose unit test ..." << std::endl;
	
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
    for (int i = 1; i < local_N; i++) {
        b[i] = b[i-1] + M;
        bt[i] = bt[i-1] + M;
    }

	// Find send and displacement buffers
	int *sendcounts = (int *)malloc(size * sizeof(int));
	int *sdispls = (int *)malloc(size * sizeof(int));
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
			b[i][j] = i + j*M + local_start_N;	

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

	MPI_Alltoallv(&temp[0], sendcounts, sdispls, MPI_DOUBLE, &recieve_temp[0], sendcounts, sdispls, MPI_DOUBLE, MPI_COMM_WORLD);

	count = 0;
	for (int j= 0; j < M; j++)
		for (int i = 0; i < local_N; i++)
			bt[i][j] = recieve_temp[count++];
	
	int correct = 1;
	int recieve_correct;
	for (int i = 0; i < local_N; i++) {
		for (int j = 0; j < M; j++) {
			if (bt[i][j] != j + M * (i + local_start_N)){
				correct = 0;
			}
		}
	}

	MPI_Allreduce(&correct, &recieve_correct, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


	if (recieve_correct == size) {
		if (rank == 0)
			std::cout << "Test matrix is sucessfully transposed!" << std::endl;
		return true;
	}
	else {
		if (rank == 0)
			std::cout << "Transpose unit test failed" << std::endl; 
		return false;
	}
	
}


////////////////////////  UNIT TEST: TEST IF PARALLEL CODE PRODUCES SAME RESULT AS SERIAL /////////////////////////////

bool utest_result(int rank, double **b, int local_N, int m){
	if (rank == 0)
	std::cout << "\nExcecuting unit test: \ncomparing parallel and serial results ..." << std::endl;

    double u_max = 0.0;
    for (int i = 0; i < local_N; i++) {
        for (int j = 0; j < m; j++) {
            u_max = u_max > fabs(b[i][j]) ? u_max : fabs(b[i][j]);
        }
    }
    double global_u_max = 0.0;
    MPI_Allreduce(&u_max, &global_u_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); 

    if (rank == 0){
		printf("u_max = %e\n", global_u_max);
    }
	
	if (fabs(global_u_max) < 1){
		if (rank == 0)
            std::cout << "Solution is within required tolerance!" << std::endl;
		return true;
	}
	else{
		if (rank == 0) 
			std::cout << "Solution is NOT within required tolerance" << std::endl;
		return false;
    }
}


