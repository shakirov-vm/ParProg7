#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

#define ISIZE 15000
#define JSIZE 15000

#define NULL_EXECUTOR 0
#define ONE_ELEMENT 1
#define TAG_2 2

// mpicc 1e_para.c -lm
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
/*	
	int distance = 1 * JSIZE + 8;
	int per_one_distance = distance / commsize + (distance % commsize) / this_rank;
	//начало измерения времени
	int jump_distance = distance - per_one_distance;

	int in_distance = 0;
	for (i = 1; i < ISIZE; i++) {
		for (j = 8; j < JSIZE; j++) {
			a[i * JSIZE + j] = sin(5 * a[(i - 1) * JSIZE - 8]);
			in_distance++;
			if (in_distance >= per_one_distance) {
				in_distance = 0;
				j = (i * JSIZE + j + jump_distance) % JSIZE;
				i = (i * JSIZE + j + jump_distance) / JSIZE;
				MPI_Barrier(MPI_COMM_WORLD);
			}
		}
	}
*/

	clock_t start_time = clock();
	j = this_rank;
	for (i = 0; i < ISIZE - 1; i++) {
		for (; j < JSIZE - 8;) {
			a[(1 + i) * JSIZE + (8 + j)] = sin(5 * a[i * JSIZE + j]);
			j += commsize;
		}
		j %= (JSIZE - 8);
	}
	clock_t end_time = clock();
	printf("Time is %lf seconds\n", ((double) end_time - start_time) / CLOCKS_PER_SEC);
	//окончание измерения времени

	if (this_rank != NULL_EXECUTOR) {
		j = this_rank;
		for (i = 0; i < ISIZE - 1; i++) {
			for (; j < JSIZE - 8;) {
				int index = (1 + i) * JSIZE + (j + 8);
				MPI_Send(&a[index], ONE_ELEMENT, MPI_DOUBLE, NULL_EXECUTOR, TAG_2, MPI_COMM_WORLD);
				j += commsize;
				MPI_Barrier(MPI_COMM_WORLD);
			}
			j %= (JSIZE - 8);
		}	
	}
	else {
		j = 0;
		for (i = 0; i < ISIZE - 1; i++) {
			for (; j < JSIZE - 8;) {
				for (int rank = NULL_EXECUTOR + 1; rank < commsize; rank++) {
					int index = (1 + i + (j + 8 + rank) / JSIZE) * JSIZE + (j + 8 + rank) % JSIZE;
					MPI_Recv(&a[index], ONE_ELEMENT, MPI_DOUBLE, rank, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
				}
				j += commsize;
				MPI_Barrier(MPI_COMM_WORLD);
			}
			j %= (JSIZE - 8);
		}
	}
/*
	for (i = 1; i < ISIZE; i++) {
		for (j = 8; j < JSIZE; j++) {
			for (int in_distance = 0; in_distance < per_one_distance; in_distance++) {
				a[i][j] = sin(5 * a[i - 1][j - 8]);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
	*/
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

	MPI_Finalize();
}