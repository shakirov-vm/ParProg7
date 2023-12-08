#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define ISIZE 15000
#define JSIZE 15000

// gcc 1e_cons.c -lm

int main()
{
	int i, j;

//	double a[ISIZE][JSIZE];
	double** a = (double**) calloc(ISIZE, sizeof(double));
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
	for (i = 1; i < ISIZE; i++){
		for (j = 8; j < JSIZE; j++){
			a[i][j] = sin(5 * a[i - 1][j - 8]);
		}
	}
	clock_t end_time = clock();
	printf("Time is %lf seconds\n", ((double) end_time - start_time) / CLOCKS_PER_SEC);
	//окончание измерения времени
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