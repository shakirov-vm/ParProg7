#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define ISIZE 15000
#define JSIZE 15000

// gcc 2d_cons.c -lm

int main(int argc, char **argv)
{
	double** a = (double**) calloc(ISIZE, sizeof(double));
	int i, j;
	if (a == NULL) {
		printf("Cannot allocate memory\n");
	}
	for (i = 0; i < ISIZE; i++) {
		a[i] = (double*) calloc(JSIZE, sizeof(double));
	}
	FILE *ff;
	for (i = 0; i < ISIZE; i++){
		for (j = 0; j < JSIZE; j++){
			a[i][j] = 10 * i + j;
		}
	}
	clock_t start_time = clock();
	//начало измерения времени
	for (i = 0; i < ISIZE - 4; i++){
		for (j = 5; j < JSIZE; j++){
			a[i][j] = sin(0.1 * a[i + 4][j - 5]);
		}
	}
	//окончание измерения времени
	clock_t end_time = clock();
	printf("Time is %lf seconds\n", ((double) end_time - start_time) / CLOCKS_PER_SEC);
	ff = fopen("result.txt","w");
	for(i = 0; i < ISIZE; i++){
		for (j = 0; j < JSIZE; j++){
			fprintf(ff,"%f ",a[i][j]);
		}
		fprintf(ff,"\n");
	}
	fclose(ff);
	free(a);
}