#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#define MAX_ITER 100

// Maximum value of the matrix element
#define MAX 100
#define TOL 0.000001




// Generate a random float number with the maximum value of max
float rand_float(int max) {
	return ((float)rand() / (float)(RAND_MAX)) * max;
}




// Calculates how many rows are going to be send to each node
int rows_per_node(int num_nodes, int n) {
	return ((int)floor(n/nodes) + 2);
}




// Gets the index from which each node should get rows
// SUPPOSITION 1: STARTS WITH NODE 0
// SUPPOSITION 2: THE RETURNED INDEX IS INCLUDED
int get_lower_index(int node_id, int rows_per_node) {
	return = node_id * (rows_per_node-2);
}




// Gets the index until which each node should get rows
// SUPPOSITION 1: STARTS WITH NODE 0
// SUPPOSITION 2: THE RETURNED INDEX IS INCLUDED
int get_upper_index(int node_id, int rows_per_node, int n) {

	int index = (node_id+1) * (rows_per_node-2) + 1;
	if (index >= n) {
		return n-1;
	}
	else {
		return index;
	}
}




// Allocate 2D matrix
void allocate_init_2Dmatrix(float ***mat, int n, int m){

	int i, j;
	*mat = (float **) malloc(n * sizeof(float *));

	for (i = 0; i < n; i++) {
		(*mat)[i] = (float *)malloc(m * sizeof(float));

		for (j = 0; j < m; j++) {
			(*mat)[i][j] = rand_float(MAX);
		}
	}
}




// Solver
void solver(float ***mat, int n, int m) {

	float diff = 0, temp;
	int done = 0, cnt_iter = 0, i, j, myrank;

  	while (!done && (cnt_iter < MAX_ITER)) {
  		diff = 0;

  		for (i = 1; i < n - 1; i++) {
  			for (j = 1; j < m - 1; j++) {
  				temp = (*mat)[i][j];
				(*mat)[i][j] = 0.2 * ((*mat)[i][j] + (*mat)[i][j - 1] + (*mat)[i - 1][j] + (*mat)[i][j + 1] + (*mat)[i + 1][j]);
				diff += abs((*mat)[i][j] - temp);
  			}
      	}

		if (diff/n/n < TOL){
			done = 1;
		}
		cnt_iter ++;
	}


	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if (myrank == 0) {
		if (done) {
			printf("Solver converged after %d iterations\n", cnt_iter);
		}
		else {
			printf("Solver not converged after %d iterations\n", cnt_iter);
		}
	}
}




int main(int argc, char *argv[]) {

	int np, myrank, n, communication;
	float **a;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	if (argc < 3) {
		if (myrank == 0) {
			printf("Call this program with two parameters: matrix_size communication \n");
			printf("\t matrix_size: Add 2 to a power of 2 (e.g. : 18, 1026)\n");
			printf("\t communication:\n");
			printf("\t\t 0: initial and final using point-to-point communication\n");
			printf("\t\t 1: initial and final using collective communication\n");
		}

		MPI_Finalize();
		exit(1);
	}

	n = atoi(argv[1]);
	communication =	atoi(argv[2]);
	printf("Matrix size = %d communication = %d\n", n, communication);


	// Calculate how many rows are going to be sent to each node
	int num_rows = rows_per_node(np, n);


	switch(communication) {
		case 0: {
			if (myrank == 0) {

				// Allocating memory for the whole matrix
				allocate_init_2Dmatrix(&a, n, n);
				int i;

				// Master sends chuncks to every other node
				for (i = 1; i < np; i++) {
					int lower_index = get_lower_index(i, num_rows);
					int upper_index = get_upper_index(i, num_rows, n);
					int num_elems = (upper_index-lower_index+1) * n;

					MPI_Send(&a[lower_index], num_elems, MPI_FLOAT, i, MPI_ANY_TAG, MPI_COMM_WORLD);
				}
			}
			else {

				// Allocating the exact memory to the rows receiving
				allocate_init_2Dmatrix(&a, num_rows, n);
				int status;

				int lower_index = get_lower_index(myrank, num_rows);
				int upper_index = get_upper_index(myrank, num_rows, n);
				int num_elems = (upper_index-lower_index+1) * n;

				// Receiving the data from the master node
				MPI_Recv(&a, num_elems, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}
			break;
		}

		/*
		case 1: {
			if (myrank == 0) {
				allocate_init_2Dmatrix(&a, n, n);
			}
			else {
				allocate_init_2Dmatrix(&a, num_rows, n);
			}

			// collective communication for scattering the matrix
			break;
		}
		*/
	}


	// --------- SOLVER ---------
	if (myrank == 0) {
		//TODO (Aqui va a haber problemas por el código del solver)
	}
	else {
		solver(&a, num_rows, n);
	}


	switch(communication) {
		case 0: {
			// p2p communication for scattering the matrix
			break;
		}
		case 1: {
			// collective communication for gathering the matrix
			break;
		}
	}

	MPI_Finalize();
	return 0;
}
