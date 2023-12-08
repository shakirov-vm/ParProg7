#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define ISIZE 15000
#define JSIZE 15000

// gcc 3g_cons.c -lm

int main(int argc, char **argv) {
	double** a = (double**) calloc(ISIZE, sizeof(double));
	double** b = (double**) calloc(ISIZE, sizeof(double));
	if (a == NULL || b == NULL) {
		printf("Cannot allocate memory\n");
		return 1;
	}
	int i, j;
	for (i = 0; i < ISIZE; i++) {
		a[i] = (double*) calloc(JSIZE, sizeof(double));
		b[i] = (double*) calloc(JSIZE, sizeof(double));
	}
	FILE *ff;
	for (i = 0; i < ISIZE; i++) {
		for (j = 0; j < JSIZE; j++) {
			a[i][j] = 10 * i +j;
			b[i][j] = 0;
		}
	}
	clock_t start_time = clock();
	//начало измерения времени
	for (i = 0; i < ISIZE; i++) {
		for (j = 0; j < JSIZE; j++) {
			a[i][j] = sin(0.005 * a[i][j]);
		}
	}
	for (i = 5; i < ISIZE; i++){
		for (j = 0; j < JSIZE - 2; j++){
			b[i][j] = a[i - 5][j + 2] * 1.5;
		}
	}
	//окончание измерения времени
	clock_t end_time = clock();
	printf("Time is %lf seconds\n", ((double) end_time - start_time) / CLOCKS_PER_SEC);

	ff = fopen("result.txt","w");
	for(i = 0; i < ISIZE; i++){
		for (j = 0; j < JSIZE; j++){
			fprintf(ff,"%f ", b[i][j]);
		}
		fprintf(ff,"\n");
	}
	fclose(ff);
}