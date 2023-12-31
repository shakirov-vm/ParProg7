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

// g++ 2d_para.c -lm -fopenmp 

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

	int k, l;
	FILE *ff;
	for (k = 0; k < ISIZE; k++){
		for (l = 0; l < JSIZE; l++){
			a[k * JSIZE + l] = 10 * k + l;
		}
	}
	auto&& start = std::chrono::high_resolution_clock::now();
	
	for (int i = 0; i < ISIZE - 4; i++) {        

        #pragma omp parallel for num_threads(prog_num_threads)
        for (int j = 5; j < JSIZE; j++) {
            a[i * JSIZE + j] = sin(0.1 * a[(i + 4) * JSIZE + j - 5]);
        }
        // implicit barrier
    }
    
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