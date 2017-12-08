//
// Implementation of the iterative Jacobi method.
//
// Given a known, diagonally dominant matrix A and a known vector b, we aim to
// to find the vector x that satisfies the following equation:
//
//		 Ax = b
//
// We first split the matrix A into the diagonal D and the remainder R:
//
//		 (D + R)x = b
//
// We then rearrange to form an iterative solution:
//
//		 x' = (b - Rx) / D
//
// More information:
// -> https://en.wikipedia.org/wiki/Jacobi_method
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <immintrin.h>
#include <omp.h>
#include "mpi.h"

static int N;
static int MAX_ITERATIONS;
static int SEED;
static float CONVERGENCE_THRESHOLD;

#define UNROLL 8

#define SEPARATOR "------------------------------------\n"

// Return the current time in seconds since the Epoch
double get_timestamp();

// Parse command line arguments to set solver parameters
void parse_arguments(int argc, char *argv[]);

// Run the Jacobi solver
// Returns the number of iterations performed
int run(float *A, float *b, float *x, float *xtmp, int rank, int size)
{
	int itr;
	int row;
	float dot;
	float diff;
	float sqdiff;
	float *ptrtmp;

	// Work out local N

	// Loop until converged or maximum iterations reached
	itr = 0;
	do {
		sqdiff = 0.0;
		// Perfom Jacobi iteration
		#pragma omp parallel for schedule(static) private(dot) shared(A,b,x,xtmp) reduction(+:sqdiff)
		for (row = 0; row < N; row += 1) {
			dot = 0.0F;

			#pragma omp simd reduction(+:dot)
			for (int col = 0; col < N; col++) {
				float val = A[(row+0)*N + col] * x[col];
				dot += val;
			}

			dot -= A[row*N + row] * x[row];
			xtmp[row] = (b[row] - dot) / A[(row)*N + row + 0];

			diff = x[row] - xtmp[row];
			diff *= diff;
			sqdiff += diff;
		}

		ptrtmp = x;
		x      = xtmp;
		xtmp   = ptrtmp;
		itr++;
	} while ((itr < MAX_ITERATIONS) && (sqrt(sqdiff) > CONVERGENCE_THRESHOLD));

	return itr;
}

void sendInitialData(float *A, float *b, float *x, int size)
{
	for (int dest = 1; dest < size; dest++) {
		MPI_Send(A, N*N, MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
		MPI_Send(b, N,   MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
		MPI_Send(x, N,   MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
	}
}

void recvInitialData(float *A, float *b, float *x)
{
	MPI_Recv(A, N*N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(b, N,   MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(x, N,   MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

int main(int argc, char *argv[])
{
	parse_arguments(argc, argv);

	// Get into MPI land

	MPI_Init(&argc, &argv);

	// Check if that worked
	int flag;
	MPI_Initialized(&flag);
	if (!flag) {
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}

	// We're here, it worked, lets find out what's going on
	char *hostname = malloc(MPI_MAX_PROCESSOR_NAME); // Process name
	int namelen = 0; // name length
	int size = 0;    // size of group
	int rank = 0;    // which one am I?
	MPI_Get_processor_name(hostname, &namelen);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	printf("Process %d of %d started, called: '%s'\n", rank, size, hostname);

	float *A    = _mm_malloc(N*N*sizeof(float), 64);
	float *b    = _mm_malloc(N*sizeof(float), 64);
	float *x    = _mm_malloc(N*sizeof(float), 64);
	float *xtmp = _mm_malloc(N*sizeof(float), 64);

	// Initialize memory
	#pragma omp parallel for schedule(static)
	for (int row = 0; row < N; row++) {
		for (int col = 0; col < N; col++) {
			A[row*N + col] = 0.0F;
		}
		b[row] = 0.0F;
		x[row] = 0.0F;
	}

	double total_start = 0.0F;
	if (rank == 0) {
		printf(SEPARATOR);
		printf("Matrix size:            %dx%d\n", N, N);
		printf("Maximum iterations:     %d\n", MAX_ITERATIONS);
		printf("Convergence threshold:  %lf\n", CONVERGENCE_THRESHOLD);
		printf(SEPARATOR);

		double total_start = get_timestamp();

		// Initialize random data
		srand(SEED);
		for (int row = 0; row < N; row++) {
			float rowsum = 0.0;
			for (int col = 0; col < N; col++) {
				float value = rand()/(float)RAND_MAX;
				A[row*N + col] = value;
				rowsum += value;
			}
			A[row*N + row] += rowsum; // Add the row sum to the diagonal
			b[row] = rand()/(float)RAND_MAX;
		}
		sendInitialData(A, b, x, size);
	} else { // Rank != 0
		recvInitialData(A, b, x);
	}

	// Ship out initialisation data

	// Run Jacobi solver
	double solve_start = get_timestamp();
	int itr = run(A, b, x, xtmp, rank, size);
	double solve_end = get_timestamp();

	float err = 0.0;
	if (rank == 0) {
		// Check error of final solution
		for (int row = 0; row < N; row++) {
			float tmp = 0.0;
			for (int col = 0; col < N; col++) {
				tmp += A[row*N + col] * x[col];
			}
			tmp = b[row] - tmp;
			err += tmp*tmp;
		}
		err = sqrt(err);
	}

	double total_end = get_timestamp();

	printf("Solution error = %lf\n", err);
	printf("Iterations     = %d\n", itr);
	printf("Total runtime  = %lf seconds\n", (total_end-total_start));
	printf("Solver runtime = %lf seconds\n", (solve_end-solve_start));
	if (itr == MAX_ITERATIONS) {
		printf("WARNING: solution did not converge\n");
	}
	printf(SEPARATOR);

	_mm_free(A);
	_mm_free(b);
	_mm_free(x);
	_mm_free(xtmp);

	// Finish with MPI
	MPI_Finalize();
	return 0;
}

double get_timestamp()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec*1e-6;
}

int parse_int(const char *str)
{
	char *next;
	int value = strtoul(str, &next, 10);
	return strlen(next) ? -1 : value;
}

float parse_float(const char *str)
{
	char *next;
	float value = strtod(str, &next);
	return strlen(next) ? -1 : value;
}

void parse_arguments(int argc, char *argv[])
{
	// Set default values
	N = 1000;
	MAX_ITERATIONS = 20000;
	CONVERGENCE_THRESHOLD = 0.0001;
	SEED = 0;

	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "--convergence") || !strcmp(argv[i], "-c")) {
			if (++i >= argc || (CONVERGENCE_THRESHOLD = parse_float(argv[i])) < 0) {
				printf("Invalid convergence threshold\n");
				exit(1);
			}
		} else if (!strcmp(argv[i], "--iterations") || !strcmp(argv[i], "-i")) {
			if (++i >= argc || (MAX_ITERATIONS = parse_int(argv[i])) < 0) {
				printf("Invalid number of iterations\n");
				exit(1);
			}
		} else if (!strcmp(argv[i], "--norder") || !strcmp(argv[i], "-n")) {
			if (++i >= argc || (N = parse_int(argv[i])) < 0) {
				printf("Invalid matrix order\n");
				exit(1);
			}
		} else if (!strcmp(argv[i], "--seed") || !strcmp(argv[i], "-s")) {
			if (++i >= argc || (SEED = parse_int(argv[i])) < 0) {
				printf("Invalid seed\n");
				exit(1);
			}
		} else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
			printf("\n");
			printf("Usage: ./jacobi [OPTIONS]\n\n");
			printf("Options:\n");
			printf("  -h  --help               Print this message\n");
			printf("  -c  --convergence  C     Set convergence threshold\n");
			printf("  -i  --iterations   I     Set maximum number of iterations\n");
			printf("  -n  --norder       N     Set maxtrix order\n");
			printf("  -s  --seed         S     Set random number seed\n");
			printf("\n");
			exit(0);
		} else {
			printf("Unrecognized argument '%s' (try '--help')\n", argv[i]);
			exit(1);
		}
	}
}
