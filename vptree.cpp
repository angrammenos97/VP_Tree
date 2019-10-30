#include "pch.h"

#include "vptree.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/////////////////////////////////
double *distance_from_last(double *X, int n, int dim)
{
	if (n == 1)
		exit(-1);
	double *d = (double*)malloc((n - 1) * sizeof(double));
	for (int i = 0; i < n - 1; i++) {
		*(d + i) = 0.0;
		for (int j = 0; j < dim; j++)
			*(d + i) += pow(*(X + i * dim + j) - *(X + (n - 1) * dim + j), 2);
		*(d + i) = sqrt(*(d + i));
	}
	return d;
}

void SWAPX(double *X, int dim, int a, int b)
{
	double *tempX = (double*)malloc(dim * sizeof(double));
	for (int j = 0; j < dim; j++)
		*(tempX + j) = *(X + a * dim + j);
	for (int j = 0; j < dim; j++)
		*(X + a * dim + j) = *(X + b * dim + j);
	for (int j = 0; j < dim; j++)
		*(X + b * dim + j) = *(tempX + j);
}

double quick_select(double *d, double *X, int *idx, int len, int k,int dim)
{
#define SWAP(a, b) { tmpd = d[a]; d[a] = d[b]; d[b] = tmpd; tmpi = idx[a]; idx[a] = idx[b]; idx[b] = tmpi;}
	int i, st, tmpi;
	double tmpd;
	for (st = i = 0; i < len - 1; i++) {
		if (d[i] > d[len - 1]) continue;
		SWAP(i, st);
		SWAPX(X, dim, i, st);
		st++;
	}
	SWAP(len - 1, st);
	SWAPX(X, dim, len - 1, st);
	return k == st ? d[st]
		: st > k ? quick_select(d, X, idx, st, k, dim)
		: quick_select(d + st, X + st * dim, idx + st, len - st, k - st, dim);
}

double median(double *X, int *idx, int n, int dim)
{
	if (n == 1)
		return 0.0;
	double *d = distance_from_last(X, n, dim);		
	return quick_select(d, X, idx, (n - 1), (n - 2) / 2, dim);
}

vptree *vpt(double *X, int *idx, int n, int dim)
{
	vptree *tree = (vptree*)malloc(sizeof(vptree));
	if (n == 0) {
		tree = NULL;
		return tree;
	}
	tree->vp = (double *)malloc(1 * dim * sizeof(double));
	tree->vp = (X + (n - 1) * dim);
	tree->md = median(X, idx, n, dim);
	tree->idx = *(idx + n - 1);
	// split and recurse
	if ((n - 1) % 2 == 0) {
		tree->inner = vpt(X, idx, (n - 1) / 2, dim);
		tree->outer = vpt((X + ((n - 1) / 2)*dim), (idx + (n - 1) / 2), (n - 1) / 2, dim);
	}
	else {
		tree->inner = vpt(X, idx, (n - 1) / 2 + 1, dim);
		tree->outer = vpt((X + ((n - 1) / 2 + 1)*dim), (idx + (n - 1) / 2 + 1), (n - 1) / 2, dim);
	}
	return tree;
}
/////////////////////////////////

vptree * buildvp(double *X, int n, int d)
{
	double *X_copy = (double*)malloc(n * d * sizeof(double));
	int *idx = (int*)malloc(n * sizeof(int));
	for (int i = 0; i < n; i++) {
		*(idx + i) = i;
		for (int j = 0; j < d; j++)
			*(X_copy + i * d + j) = *(X + i * d + j);
	}
	return vpt(X_copy, idx, n, d);
}

vptree * getInner(vptree * T)
{
	return T->inner;
}
vptree * getOuter(vptree * T)
{
	return T->outer;
}
double getMD(vptree * T) 
{
	return T->md;
}
double * getVP(vptree * T)
{
	return T->vp;
}
int getIDX(vptree * T) 
{
	return T->idx;
}
