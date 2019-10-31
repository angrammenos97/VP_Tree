#include "vptree.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SWAP(x, y) { double temp = *x; *x = *y; *y = temp; }

void swapPoints(double* array1, double* array2, int dim) {
	double tmp;
	int i;

	for (i = 0; i < dim; i++) {
		tmp = array1[i];
		array1[i] = array2[i];
		array2[i] = tmp;
	}
}

//calculates the distance between elements of the given part of the array assuming the last element is the vp and using the
//index array as reference
//validated
double * calculateDist(int dim, int size, double arrayList[size][dim]) {


	double* dist = (double*)malloc((size - 1) * sizeof(double));
	int i;

	//using i to dereference index table and then access the original holder table of nxd
	for (i = 0; i < size - 1; ++i) {
		double distSqrd = 0.0;
		int y;

		//for loop to calculate distance for all dimensions
		for (y = 0; y < dim; ++y) {
			distSqrd += pow(arrayList[i][y] - arrayList[size - 1][y], 2);
		}

		dist[i] = sqrt(distSqrd);

	}

	return dist;
}



// Partition using Lomdto partition scheme and parallel update of the initial index table
//validated
double* partition(int size, int dim, int* a, double list[size][dim], double* inner, double* outer, double* pivotIndex)
{
	// Pick pivotIndex as pivot from the array
	double pivot = *pivotIndex;
	int* pivotA = a + (pivotIndex - inner);
	double(*pivotArray)[dim] = list + (pivotIndex - inner);

	// Move pivot to end
	SWAP(pivotIndex, outer);
	SWAP(pivotA, (a + size - 1));
	swapPoints(*pivotArray, *(list + size - 1), dim);

	// elements less than pivot will be pushed to the inner of pIndex
	// elements more than pivot will be pushed to the outer of pIndex
	// equal elements can go either way
	pivotIndex = inner;
	pivotA = a;
	pivotArray = list;

	int i;

	// each time we finds an element less than or equal to pivot, pIndex
	// is incremented and that element would be placed before the pivot and index table is also updated.
	for (i = 0; i < size - 1; i++)
	{
		if (inner[i] <= pivot)
		{
			SWAP((inner + i), pivotIndex);
			SWAP((a + i), pivotA);
			swapPoints(*(list + i), *pivotArray, dim);
			pivotIndex++;
			pivotA++;
			pivotArray++;
		}
	}

	// Move pivot to its final place
	SWAP(pivotIndex, outer);
	SWAP(pivotA, (a + size - 1));
	swapPoints(*(pivotArray), *(list + size - 1), dim);

	// return pIndex (index of pivot element)
	return pivotIndex;
}

// Returns the k-th smallest element of list within inner..outer
// (i.e. inner <= k <= outer) while updating initial index table using lomdto partition
//validated
double* quickselect(int size, int dim, int* a, double list[size][dim], double* inner, double* outer, int k)
{
	// If the array contains only one element, return that element
	if (inner == outer)
		return inner;

	// select a pivotIndex between inner and outer
	double*  pivotIndex = inner + (rand() % (outer - inner + 1));

	pivotIndex = partition(outer - inner + 1, dim, a, list, inner, outer, pivotIndex);

	// The pivot is in its final sorted position
	if ((inner + k - 1) == pivotIndex)
		return pivotIndex;

	// if k is less than the pivot index
	else if ((inner + k - 1) < pivotIndex)
		return quickselect((pivotIndex - inner), dim, a, list, inner, pivotIndex - 1, k);

	// if k is more than the pivot index
	else
		return quickselect(outer - pivotIndex, dim, a + (pivotIndex - inner) + 1, list + (pivotIndex - inner) + 1, pivotIndex + 1, outer, k - (pivotIndex - inner) - 1);
}


//finds the median value - median calculated as the first of the middle elements in case of even
//# of elements - and returns that value while having rearranged the index table properly for further use in recursion
//validated
double findMedian(int dim, int size, double arrayList[size][dim], int index[size]) {

	//checks whether table has only one element
	if (size == 1)
		return 0.0;

	//calculates median and updates index table
	double* dist = calculateDist(dim, size, arrayList);
	double md = *quickselect(size - 1, dim, index, arrayList, dist, dist + (size - 2), (size - 1) / 2);

	return md;
}


//creates Vptree assuming vp is the last element in index and calls recursively until
//there are no points inner
//validated
vptree * createVpTree(int dim, int size, int index[size], double list[size][dim]) {

	if (size == 0)
		return NULL;

	vptree* node = (vptree*)malloc(sizeof(vptree));
	node->idx = index[size - 1];
	node->vp = list[size - 1];
	node->md = findMedian(dim, size, list, index);

	//calls recursively taking into consideration whether size is
	//odd or even number
	if (size % 2 != 0) {
		node->inner = createVpTree(dim, (size - 1) / 2, index, list);
		node->outer = createVpTree(dim, (size - 1) / 2, index + (size - 1) / 2, list + (size - 1) / 2);
	}
	else {
		node->inner = createVpTree(dim, (size - 1) / 2 + 1, index, list);
		node->outer = createVpTree(dim, (size - 1) / 2, index + (size - 1) / 2 + 1, list + (size - 1) / 2 + 1);
	}

	return node;
}

//prints VpTree with reference to its nodes
void printTree(vptree* node, int* counter) {

	if (node == NULL) {
		printf("\n");
		return;
	}
	(*counter)++;
	int temp = *counter;
	printf("printing inner->%d \n", temp);
	printTree(node->inner, counter);

	printf("%lf %lf %lf ->%d\n", (node->vp)[0], (node->vp)[1], (node->vp)[2], temp);

	printf("printing outer->%d\n", temp);
	printTree(node->outer, counter);
}

//function written only to meet the header criteria of the project
//thus simply making use of createVpTree plus creating the index table
//validated
vptree* buildvp(double* X, int n, int d) {

	double* list = (double*)malloc(sizeof(double)*n*d);
	memcpy(list, X, sizeof(double)*n*d);


	int* index = (int*)malloc(sizeof(int)*n);
	int i;
	for (i = 0; i < n; ++i) {
		index[i] = i;
	}

	vptree* root = createVpTree(d, n, index, (double**)list);
	return root;
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
