#include "pch.h"

#include "vptree.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/////////////////////////////////
double *distance_from_last(double *X, int *idx, int n, int dim)
{
	if (n == 1)
		exit(-1);

	double *d = (double*)malloc((n - 1) * sizeof(double));
	memset(d, 0, (n - 1) * sizeof(double));	// set hole array with zeroes
	for (int i = 0; i < n - 1; i++) {
		for (int j = 0; j < dim; j++)
			*(d + i) += pow(*(X + *(idx + i) * dim + j) - *(X + *(idx + n - 1) * dim + j), 2);
		*(d + i) = sqrt(*(d + i));
	}
	return d;
}

double quick_select(double *v, int *idx, int len, int k)
{
#define SWAP(a, b) { tmpd = v[a]; v[a] = v[b]; v[b] = tmpd; tmpi = idx[a]; idx[a] = idx[b]; idx[b] = tmpi;}
	int i, st, tmpi;
	double tmpd;
	for (st = i = 0; i < len - 1; i++) {
		if (v[i] > v[len - 1]) continue;
		SWAP(i, st);
		st++;
	}
	SWAP(len - 1, st);
	return k == st ? v[st]
		: st > k ? quick_select(v, idx, st, k)
		: quick_select(v + st, idx, len - st, k - st);
}

double median(double *d, int *idx, int n)
{
	if (n == 0)
		exit(-1);
	else
	{
		if ((n % 2) == 1)
			return quick_select(d, idx, n, n / 2);
		else
			return ((quick_select(d, idx, n, n / 2 - 1) + quick_select(d, idx, n, n / 2)) / 2.0);
	}
}

vptree *vpt(double *X, int *idx, int n, int dim)
{
	vptree *tree = (vptree*)malloc(sizeof(vptree));
	tree->vp = (double *)malloc(1 * dim * sizeof(double));
	if (n == 0) {
		tree = NULL;
		return tree;
	}
	else if (n == 1) {
		tree->vp = (X + *(idx + n - 1) * dim);
		tree->md = 0.0;
		tree->idx = *(idx + n - 1);
		tree->inner = NULL;
		tree->outer = NULL;
	}
	else {
		tree->vp = (X + *(idx + n - 1) * dim);
		double *d = distance_from_last(X, idx, n, dim);
		tree->md = median(d, idx, n - 1);
		tree->idx = *(idx + n - 1);
		// split and recurse	
		if ((n - 1) % 2 == 0) {
			tree->inner = vpt(X, idx, (n - 1) / 2, dim);
			tree->outer = vpt(X, (idx + (n - 1) / 2), (n - 1) / 2, dim);
		}
		else {
			tree->inner = vpt(X, idx, (n - 1) / 2 + 1, dim);
			tree->outer = vpt(X, (idx + (n - 1) / 2 + 1), (n - 1) / 2, dim);
		}
	}
	return tree;
}
/////////////////////////////

vptree * buildvp(double *X, int nop, int d)
{
	int *idx = (int*)malloc(nop * sizeof(int));
	for (int i = 0; i < nop; i++)
		*(idx + i) = i;

	return vpt(X, idx, nop, d);
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
