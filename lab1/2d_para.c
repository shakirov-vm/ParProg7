#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <chrono>
#include <iostream>
#include <fstream>
#define ISIZE 15000
#define JSIZE 15000

// gcc 2d_para.c -lm -fopenmp 

int main(int argc, char **argv)
{
	if (argc != 2) {
		printf("Need 2 arg\n");
		return 1;
	}

	long prog_num_threads = atoi(argv[1]);

	double* a = (double*) calloc(ISIZE * JSIZE, sizeof(double));

	if (a == NULL) {
		printf("Cannot allocate memory\n");
	}
	double* stable_a = (double*) calloc(ISIZE * JSIZE, sizeof(double));
	if (stable_a == NULL) {
		printf("Cannot allocate memory\n");
	}

	int k, l;
	FILE *ff;
	for (k = 0; k < ISIZE; k++){
		for (l = 0; l < JSIZE; l++){
			a[k * JSIZE + l] = 10 * k + l;
			stable_a[k * JSIZE + l] = a[k * JSIZE + l];
		}
	}
	auto&& start = std::chrono::high_resolution_clock::now();
//	clock_t start_time = clock();
	//начало измерения времени
#pragma omp parallel num_threads(prog_num_threads)
{	
	int i, j;
	int this_rank = omp_get_thread_num();

	int distance = 4 * (JSIZE - 5);
	int per_one_steps = distance / prog_num_threads + (distance % prog_num_threads) * this_rank;
	int steps = 0;

	j = this_rank;
	for (i = 0; i < ISIZE - 4; i++){
		for (; j < JSIZE - 5;) {
			a[i * JSIZE + j + 5] = sin(0.1 * a[(i + 4) * JSIZE + j]);
			j += prog_num_threads;
			steps++;
		}
		j %= (JSIZE - 5);
		if (steps >= per_one_steps) {
			steps = 0;
			#pragma omp barrier
		}
	}
}
/*
	for (i = 0; i < ISIZE - 4; i++){
		for (j = 0; j < JSIZE - 5; j++){
			a[i][j + 5] = sin(0.1 * stable_a[i + 4][j]);
		}
	}
*/
	//окончание измерения времени
//	clock_t end_time = clock();
	auto&& end = std::chrono::high_resolution_clock::now();
	auto&& passed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//	printf("Time is %lf seconds\n", ((double) end_time - start_time) / CLOCKS_PER_SEC);
    std::cerr << "Time in ms: " << passed << std::endl;
	ff = fopen("result_para.txt","w");
	for(k = 0; k < ISIZE; k++){
		for (l = 0; l < JSIZE; l++){
			fprintf(ff,"%f ",a[k * JSIZE + l]);
		}
		fprintf(ff,"\n");
	}
	fclose(ff);
	free(a);
}