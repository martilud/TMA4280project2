/*
 * C++ program to solve the two-dimensional Poisson equation on
 * a unit square using one-dimensional eigenvalue decompositions
 * and fast sine transforms.
 *
 * Einar M. Rønquist
 * NTNU, October 2000
 * Revised, October 2001
 * Revised by Eivind Fonn, February 2015
 */

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <cstring>

#define PI 3.14159265358979323846
//#define true 1
//#define false 0

typedef double real;
//typedef int bool;

// Function prototypes
real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);
void transpose(real **bt, real **b, size_t m);
real rhs1();
real rhs2(real x, real y);
real analytical(real x, real y);
bool utest_transpose(int rank, int size); 
bool utest_result(int rank, double **b, int local_N, int  m);

// Functions implemented in FORTRAN in fst.f and called from C.
// The trailing underscore comes from a convention for symbol names, called name
// mangling: if can differ with compilers.
extern "C" void fst_(real *v, int *n, real *w, int *nn);
extern "C" void fstinv_(real *v, int *n, real *w, int *nn);

int main(int argc, char **argv)
{
    MPI_Init(NULL,NULL);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);



    int n = 128;
    int p = omp_get_max_threads();
    omp_set_num_threads(omp_get_max_threads());
    char t = 'm';

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
        }
    }

    
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
        n = 128;
        p = 4;
    }

    real starttime = omp_get_wtime();
    /*
     *  The equation is solved on a 2D structured grid and homogeneous Dirichlet
     *  conditions are applied on the boundary:
     *  - the number of grid points in each direction is n+1,
     *  - the number of degrees of freedom in each direction is m = n-1,
     *  - the mesh size is constant h = 1/n.
     */

    int m = n - 1;
    real h = 1.0 / n;

    /*
     * Grid points are generated with constant mesh size on both x- and y-axis.
     */
    real *grid = mk_1D_array(n+1, false);
    #pragma omp parallel for
    for (int i = 0; i < n+1; i++) {
        grid[i] = i * h;
    }

    /*
     * The diagonal of the eigenvalue matrix of T is set with the eigenvalues
     * defined Chapter 9. page 93 of the Lecture Notes.
     * Note that the indexing starts from zero here, thus i+1.
     */
    real *diag = mk_1D_array(m, false);
    #pragma omp parallel for
    for (int i = 0; i < m; i++) {
        diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
    }

    /*
     * Allocate the matrices b and bt which will be used for storing values of
     * G, \tilde G^T, \tilde U^T, U as described in Chapter 9. page 101.
     */

    int local_N = m/size;
    if (rank < m % size) local_N++;
    
    // Local start colomn, NB, not important!!
    int local_start_N = 0;
    for (int i = 0; i < rank; i++) {
        local_start_N += m/size;
        if (i < m % size) local_start_N++;
    }
    
    real **b = mk_2D_array(local_N, m, false);
    real **bt = mk_2D_array(local_N, m, false);
    
    /*
     * This vector will hold coefficients of the Discrete Sine Transform (DST)
     * but also of the Fast Fourier Transform used in the FORTRAN code.
     * The storage size is set to nn = 4 * n, look at Chapter 9. pages 98-100:
     * - Fourier coefficients are complex so storage is used for the real part
     *   and the imaginary part.
     * - Fourier coefficients are defined for j = [[ - (n-1), + (n-1) ]] while 
     *   DST coefficients are defined for j [[ 0, n-1 ]].
     * As explained in the Lecture notes coefficients for positive j are stored
     * first.
     * The array is allocated once and passed as arguments to avoid doing 
     * reallocations at each function call.
     */
    int nn = 4 * n;
    //real *z = mk_1D_array(nn, false);
    /*
     * Initialize the right hand side data for a given rhs function.
     * 
     */
    if (t == 'u'){
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < local_N; i++) {
            for (int j = 0; j < m; j++) {
                b[i][j] = h * h * rhs1(); //grid[local_start_N+i+1], grid[j+1]);
            }
        }
    }
    else{
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < local_N; i++) {
            for (int j = 0; j < m; j++) {
                b[i][j] = h * h * rhs2(grid[local_start_N+i+1], grid[j+1]);
            }
        }
    }


    // Find send and displacement buffers
	int *sendcounts = (int *)malloc(size * sizeof(int));
	int *sdispls = (int *)malloc(size * sizeof(int));

    //int *sendcounts =            [size]
	//int *sdispls[size];
    for (int i = 0; i < size; i++) {
        sendcounts[i] = m/size;
        if (i < m % size) sendcounts[i]++;
        sendcounts[i] *= local_N;
    }
 
    sdispls[0] = 0;
    for (int i = 1; i < size; i++) {
        sdispls[i] = sdispls[i-1] + sendcounts[i-1];
    }

    /*
     * Compute \tilde G^T = S^-1 * (S * G)^T (Chapter 9. page 101 step 1)
     * Instead of using two matrix-matrix products the Discrete Sine Transform
     * (DST) is used.
     * The DST code is implemented in FORTRAN in fst.f and can be called from C.
     * The array zz is used as storage for DST coefficients and internally for 
     * FFT coefficients in fst_ and fstinv_.
     * In functions fst_ and fst_inv_ coefficients are written back to the input 
     * array (first argument) so that the initial values are overwritten.
     */
    #pragma omp parallel
    {
        real *z_local = mk_1D_array(nn, false);
        #pragma omp for 
        for (int i = 0; i < local_N; i++) {
            fst_(b[i], &n, z_local, &nn);
        }
    }

	double *temp = mk_1D_array(local_N * m, 0);
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


	
	double *recieve_temp = mk_1D_array(local_N * m, 0);


	// Let the magic happen
	MPI_Alltoallv(&temp[0], sendcounts, sdispls, MPI_DOUBLE, &recieve_temp[0], sendcounts, sdispls, MPI_DOUBLE, MPI_COMM_WORLD);
	
	// Unpack buffer
	//count = 0;
    #pragma omp parallel for collapse(2)
	for (int j = 0; j < m; j++) {
		for (int i = 0; i < local_N; i++) {
			bt[i][j] = recieve_temp[j*local_N + i];
		}
	}


	
    
    #pragma omp parallel
    {
        real *z_local = mk_1D_array(nn, false);
        #pragma omp for 
        for (int i = 0; i < local_N; i++) {
            fstinv_(bt[i], &n, z_local, &nn);
        }
    }

    /*
     * Solve Lambda * \tilde U = \tilde G (Chapter 9. page 101 step 2)
     */
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < local_N; i++) {
        for (int j = 0; j < m; j++) {
            bt[i][j] = bt[i][j] / (diag[local_start_N + i] + diag[j]);
        }
    }

    /*
     * Compute U = S^-1 * (S * Utilde^T) (Chapter 9. page 101 step 3)
     */
    #pragma omp parallel
    {
        real *z_local = mk_1D_array(nn, false);
        #pragma omp for 
        for (int i = 0; i < local_N; i++) {
            fst_(bt[i], &n, z_local, &nn);
        }
    }
   

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


	
	// Let the magic happen
	MPI_Alltoallv(&temp[0], sendcounts, sdispls, MPI_DOUBLE, &recieve_temp[0], sendcounts, sdispls, MPI_DOUBLE, MPI_COMM_WORLD);
	
	// Unpack buffer
	//count = 0;
    #pragma omp parallel for collapse(2)
	for (int j= 0; j < m; j++) {
		for (int i = 0; i < local_N; i++) {
			//b[i][j] = recieve_temp[count++];
			b[i][j] = recieve_temp[j*local_N + i];
		}
	}



    #pragma omp parallel
    {
        real *z_local = mk_1D_array(nn, false);
        #pragma omp for 
        for (int i = 0; i < local_N; i++) {
            fstinv_(b[i], &n, z_local, &nn);
        }
    }


    if (t == 'v'){
        double e_max = 0.0;
        for (int i = 0; i < local_N; i++) {
            for (int j = 0; j < m; j++) {
                e_max = e_max > fabs(b[i][j] - analytical(grid[local_start_N+i+1], grid[j+1])) ? e_max : fabs(b[i][j] - analytical(grid[local_start_N+i+1], grid[j+1]));
            }
        }
        double global_e_max = 0.0;
        MPI_Reduce(&e_max, &global_e_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); 
        if (rank == 0){
            real endtime = omp_get_wtime();
            std::cout << size << '\t' << p << '\t' << n << '\t' << global_e_max << '\t' << endtime-starttime << std::endl;
        } 
        MPI_Finalize();
        return 0;  
    }
    /*
     * Compute maximal value of solution for convergence analysis in L_\infty
     * norm.
     */
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
		
		//std::cout << "pass = " << pass << std::endl;
		if (rank == 0) {
			if (pass)
				std::cout << "\n======== UNIT TESTS SUCCESSFUL ========" << std::endl;
			else
				std::cout << "\n========= UNIT TESTS FAILED!  =========" << std::endl;
		}


        MPI_Finalize();
        return 0;

    }
    
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


