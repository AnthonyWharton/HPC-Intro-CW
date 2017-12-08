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

void chunker(int *chunk, int *chunkSize, int rank, int size) {
	// Work out how many rows this thread should calculate
	// Share any unequally divided threads amongst last workers
	(*chunkSize) = N / size;
	(*chunk) = (*chunkSize) * rank ;
	int mod = N % size;
	if (mod != 0) {
		if (rank >= size - mod) {
			(*chunkSize)++;
			(*chunk) += rank - (size - mod);
		}
	}
	(*chunkSize)--;
	// printf("(%d) CHUNK: %d -> %d (%d)", rank, chunk, chunk+chunkSize, chunkSize);
}

// void masterInLoop(float *x, float *sqdiff, int *chunks, int *chunkSizes, int size, int itr) {
//
// }
//
// void workerInLoop(float *x, float *sqdiff, int chunk, int chunkSize, int itr) {
//
// }

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
	int chunk;
	int chunkSize;
	int *chunks;     // Only used by master
	int *chunkSizes; // Only used by master

	if (rank == 0) {
		chunks     = _mm_malloc((size-1)*sizeof(int), 64);
		chunkSizes = _mm_malloc((size-1)*sizeof(int), 64);

		for (int i = 0; i < size ; i++) {
			chunker(&chunks[i], &chunkSizes[i], i, size);
			printf("(%d) CHUNK: %d -> %d (%d)\n", i, chunks[i], chunks[i]+chunkSizes[i], chunkSizes[i]);
		}
	}

	chunker(&chunk, &chunkSize, rank, size);

	// Loop until converged or maximum iterations reached
	itr = 0;
	do {
		sqdiff = 0.0F;
		// Perfom Jacobi iteration
		#pragma omp parallel for schedule(static) private(dot) shared(A,b,x,xtmp) reduction(+:sqdiff)
		for (row = chunk; row < (chunk+chunkSize); row += 1) {
			dot = 0.0F;

			#pragma omp simd reduction(+:dot)
			for (int col = 0; col < N; col++) {
				float val = A[row*N + col] * x[col];
				dot += val;
			}

			dot -= A[row*N + row] * x[row];
			xtmp[row] = (b[row] - dot) / A[(row)*N + row];
			diff = x[row] - xtmp[row];
			diff *= diff;
			sqdiff += diff;
		}

		if (rank == 0) {
			// float sqdiffacc = 0.0F;
			for (int remote = 1; remote < size; remote++) {
				// printf("%f, ", x[chunks[remote]]);
				MPI_Recv(&xtmp[chunks[remote]], chunkSizes[remote], MPI_FLOAT, remote, itr, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// MPI_Recv(&sqdiffacc,            1,                  MPI_FLOAT, remote, itr, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// printf("%f, ", sqdiff);
				// sqdiff += sqdiffacc;
				// printf("%f, %f\n", sqdiffacc, sqdiff);
			}

			sqdiff = 0;
			for (int i = 0; i < N; i++) {
				diff = x[row] - xtmp[row];
				sqdiff += diff * diff;
			}

			for (int remote = 1; remote < size; remote++) {
				// printf("%d) R%d, %f %f %f %f %f %f %f %f\n", itr, rank, x[0], x[1], x[498], x[499], x[500], x[501], x[998], x[999]);
				MPI_Send(&xtmp[0], N, MPI_FLOAT, remote, itr, MPI_COMM_WORLD);
				MPI_Send(&sqdiff,  1, MPI_FLOAT, remote, itr, MPI_COMM_WORLD);
			}
		} else {
			// Send back portion of data worked on and sqdiff
			// printf("W %f %f %f %f %f %f\n", x[chunk], x[chunk+1], x[chunk+2], x[chunk+3], x[chunk+4], x[chunk+5]);
			MPI_Send(&xtmp[chunk], chunkSize, MPI_FLOAT, 0, itr, MPI_COMM_WORLD);
			// MPI_Send(&sqdiff,      1,         MPI_FLOAT, 0, itr, MPI_COMM_WORLD);

			// Await updates
			MPI_Recv(&xtmp[0], N, MPI_FLOAT, 0, itr, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			// printf("%d) R%d, %f %f %f %f %f %f %f %f\n", itr, rank, x[0], x[1], x[498], x[499], x[500], x[501], x[998], x[999]);
			MPI_Recv(&sqdiff,  1, MPI_FLOAT, 0, itr, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		ptrtmp = x;
		x      = xtmp;
		xtmp   = ptrtmp;
		itr++;

	} while ((itr < MAX_ITERATIONS) && (sqrt(sqdiff) > CONVERGENCE_THRESHOLD));

	if (rank == 0) {
		_mm_free(chunks);
		_mm_free(chunkSizes);
	}

	return itr;
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
		printf("ABORTED");
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
		b[row]    = 0.0F;
		x[row]    = 0.0F;
		xtmp[row] = 0.0F;
	}

	if (rank == 0) {
		printf(SEPARATOR);
		printf("Matrix size:            %dx%d\n", N, N);
		printf("Maximum iterations:     %d\n", MAX_ITERATIONS);
		printf("Convergence threshold:  %lf\n", CONVERGENCE_THRESHOLD);
		printf(SEPARATOR);
	}

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

	if (rank == 0) {
		printf("Solution error = %lf\n", err);
		printf("Iterations     = %d\n", itr);
		printf("Total runtime  = %lf seconds\n", (total_end-total_start));
		printf("Solver runtime = %lf seconds\n", (solve_end-solve_start));
		if (itr == MAX_ITERATIONS) {
			printf("WARNING: solution did not converge\n");
		}
		printf(SEPARATOR);
	}

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
