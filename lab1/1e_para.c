#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <chrono>
#include <iostream>
#include <fstream>

#define ISIZE 5000
#define JSIZE 5000

#define NULL_EXECUTOR 0
#define ONE_ELEMENT 1
#define TAG_2 2

// mpic++ 1e_para.c -lm
// mpirun -np 2 a.out 

int main(int argc, char** argv)
{
	int i, j;

	MPI_Init(&argc, &argv);

	int commsize, this_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &commsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);

	double* a = (double*) calloc(ISIZE * JSIZE, sizeof(double));
	if (a == NULL) {
		printf("Cannot allocate memory\n");
	}

	FILE *ff;
	for (i = 0; i < ISIZE; i++){
		for (j = 0; j < JSIZE; j++){
			a[i * JSIZE + j] = 10 * i + j;
		}
	}

#define JSIZE_CALC (JSIZE - 8)

    int jsize_per_thread = JSIZE_CALC / commsize;
    int line_start = jsize_per_thread * this_rank;
    int line_end = jsize_per_thread * (this_rank + 1);
	if (this_rank == commsize - 1) line_end = JSIZE_CALC;

    int per_thread_size = line_end - line_start;
    double per_thread_array[JSIZE];

	int* recvcnts = (int*) calloc (commsize, sizeof(int));
	int* displs = (int*) calloc (commsize, sizeof(int));

	MPI_Barrier(MPI_COMM_WORLD);
	for (int k = 0; k < commsize - 1; k++) {
		recvcnts[k] = per_thread_size;
	}
	recvcnts[commsize - 1] = per_thread_size + JSIZE_CALC % commsize;

    for (int k = 0; k < commsize; k++) {
        displs[k] = 8 + k * per_thread_size; 
    }

	auto&& start = std::chrono::high_resolution_clock::now();

	for (i = 0; i < ISIZE - 1; i++) {
        for (j = line_start; j < line_end; j++) {
            per_thread_array[j + 8] = std::sin(5 * a[i * JSIZE + j]);
        }

        MPI_Allgatherv(&per_thread_array[line_start + 8], per_thread_size, MPI_DOUBLE, 
        			&a[(i + 1) * JSIZE], recvcnts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
    }

	auto&& end = std::chrono::high_resolution_clock::now();
	auto&& passed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	MPI_Barrier(MPI_COMM_WORLD);

    std::cerr << "Time in ms: " << passed << std::endl;

	if (this_rank == NULL_EXECUTOR) {
		ff = fopen("result_para.txt","w");
		for(i = 0; i < ISIZE; i++){
			for (j = 0; j < JSIZE; j++){
				fprintf(ff,"%f ",a[i * JSIZE + j]);
			}
			fprintf(ff,"\n");
		}
		fclose(ff);
	}
	free(a);
	free(recvcnts);
	free(displs);

	MPI_Finalize();
}