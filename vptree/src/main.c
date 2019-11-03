#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "sys/time.h"
#include "vptree.h"

#define DefaultNumPoints 1000000
#define DefaultDim 20

int nop = DefaultNumPoints;
int dim = DefaultDim;
_Bool matlab = 0;

struct timeval startwtime, endwtime;

void help(int argc, char *argv[]);
void export_data(FILE *file, double *X);
void export_struct(FILE *file, vptree *root, const char *root_name, int str_size);

int main(int argc, char *argv[])
{
	help(argc, argv);
	printf("Running with values n=%i and d=%i\n", nop, dim);

	srand((unsigned int)time(NULL));
	FILE *data = NULL;

	// Generate random point set
	printf("Generating random data set. ");
	gettimeofday(&startwtime, NULL);
	double *X = (double *)malloc(nop * dim * sizeof(double));
	for (int i = 0; i < nop; i++) {
		for (int j = 0; j < dim; j++)
			*(X + i * dim + j) = ((double)rand() / (RAND_MAX));
	}
	gettimeofday(&endwtime, NULL);
	double p_time = (double)((endwtime.tv_usec - startwtime.tv_usec) / 1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
	printf("DONE in %fsec!\n", p_time);

	if (matlab) {
		printf("Writting data set to data.m. ");
		gettimeofday(&startwtime, NULL);
		data = fopen("./data.m", "w");
		export_data(data, X);
		fclose(data);
		gettimeofday(&endwtime, NULL);
		p_time = (double)((endwtime.tv_usec - startwtime.tv_usec) / 1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
		printf("DONE in %fsec!\n", p_time);
	}

	// Build search tree
	printf("Building search tree. ");
	gettimeofday(&startwtime, NULL);
	vptree *tree = buildvp(X, nop, dim);
	gettimeofday(&endwtime, NULL);
	p_time = (double)((endwtime.tv_usec - startwtime.tv_usec) / 1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
	printf("DONE in %fsec!\n", p_time);

	if (matlab) {
		printf("Appending tree to data.m ");
		gettimeofday(&startwtime, NULL);
		data = fopen("./data.m", "a");
		export_struct(data, tree, "tree", 5);
		fclose(data);
		gettimeofday(&endwtime, NULL);
		p_time = (double)((endwtime.tv_usec - startwtime.tv_usec) / 1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
		printf("DONE in %fsec!\n", p_time);
	}

	printf("Exiting\n");
	free(X);
	return 0;
}


void help(int argc, char *argv[])
{
	if (argc > 1) {
		for (int i = 1; i < argc; i += 2) {
			if (*argv[i] == '-') {
				if (*(argv[i] + 1) == 'n')
					nop = atoi(argv[i + 1]);
				else if (*(argv[i] + 1) == 'd')
					dim = atoi(argv[i + 1]);
				else if (*(argv[i] + 1) == 'm') {
					matlab = 1;
					i--;
				}
				else {
					help(1, argv);
					return;
				}
			}
			else {
				help(1, argv);
				return;
			}
		}
		return;
	}
	printf("Flags to use:\n");
	printf("-n [Number] :Number of points (default:%i)\n", DefaultNumPoints);
	printf("-d [Dimension] :Dimension of the space (default: %i)\n", DefaultDim);
	printf("-m :Print results into data.m file to evaluate in MATLAB\n");
}

void export_data(FILE *file, double *X)
{
	fprintf(file, "n = %i;\n", nop);
	fprintf(file, "dim = %i;\n", dim);
	fprintf(file, "X=[");
	for (int i = 0; i < nop; i++) {
		fprintf(file, "[");
		for (int j = 0; j < dim; j++)
			fprintf(file, "%lf ", *(X + i * dim + j));
		fprintf(file, "]; ");
	}
	fprintf(file, "];\n");
}

void export_struct(FILE *file, vptree *root, const char *root_name, int str_size)
{
	if (root == NULL) {
		fprintf(file, "%s = [];\n", root_name);
		return;
	}
	else {
		fprintf(file, "%s.vp = [", root_name);
		for (int j = 0; j < dim; j++)
			fprintf(file, "%lf ", *(getVP(root) + j));
		fprintf(file, "];\n");
		fprintf(file, "%s.md = %lf;\n", root_name, getMD(root));
		fprintf(file, "%s.idx = %i;\n", root_name, getIDX(root) + 1);
		char *tmpi = (char*)malloc((str_size + 6) * sizeof(char));
		memcpy(tmpi, root_name, str_size - 1);
		memcpy(tmpi + str_size - 1, ".inner", 7);
		char *tmpo = (char*)malloc((str_size + 6) * sizeof(char));
		memcpy(tmpo, root_name, str_size - 1);
		memcpy(tmpo + str_size - 1, ".outer", 7);
		export_struct(file, getInner(root), tmpi, str_size + 6);
		export_struct(file, getOuter(root), tmpo, str_size + 6);
	}
}
