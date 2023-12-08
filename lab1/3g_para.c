#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define ISIZE 15000
#define JSIZE 15000

// gcc 3g_para.c -lm -fopenmp 

int main(int argc, char **argv) {
	
	if (argc != 2) {
		printf("Need 2 arg\n");
		return 1;
	}

	long prog_num_threads = atoi(argv[1]);

	double* a = (double*) calloc(ISIZE * JSIZE, sizeof(double));
	double* b = (double*) calloc(ISIZE * JSIZE, sizeof(double));
	int i, j;
	FILE *ff;
	for (i = 0; i < ISIZE; i++) {
		for (j = 0; j < JSIZE; j++) {
			a[i * JSIZE + j] = 10 * i +j;
			b[i * JSIZE + j] = 0;
		}
	}

	clock_t start_time = clock();
	//начало измерения времени
#pragma omp parallel num_threads(prog_num_threads)
{	
	int i, j;
	int this_rank = omp_get_thread_num();

	int full_size = ISIZE * JSIZE;
	int start_ind = full_size / prog_num_threads * this_rank;
	int finish_ind = full_size / prog_num_threads * (this_rank + 1);

	if (this_rank == prog_num_threads - 1) {
		finish_ind = full_size;
	}

	for (i = start_ind; i < finish_ind; i++) {
		a[i] = sin(0.005 * a[i]);
	}

	#pragma omp barrier

	j = this_rank;
	for (i = 0; i < ISIZE - 5; i++) {
		for (; j < JSIZE - 2;) {
			b[(i + 5) * JSIZE + j] = a[i * JSIZE + j + 2] * 1.5;
			j += prog_num_threads;
		}
		j %= (JSIZE - 2);
	}
}
	//окончание измерения времени
	clock_t end_time = clock();
	printf("Time is %lf seconds\n", ((double) end_time - start_time) / CLOCKS_PER_SEC);

	ff = fopen("result_para.txt","w");
	for(i = 0; i < ISIZE; i++){
		for (j = 0; j < JSIZE; j++){
			fprintf(ff,"%f ", b[i * JSIZE + j]);
		}
		fprintf(ff,"\n");
	}
	fclose(ff);
}