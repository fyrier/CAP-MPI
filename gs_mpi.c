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




// Calculates how many rows are given, as maximum, to each node
int get_max_rows(int num_nodes, int n) {
	return ((int) floor(n / num_nodes) + 2);
}




// Calculates how many rows are going to a given node
int get_node_rows(int lower_bound, int upper_bound) {
	return upper_bound - lower_bound + 1;
}




// Calculates how many elements are going to a given node
int get_node_elems(int node_rows, int n) {
	return node_rows * n;
}




// Gets the index from which each node should get rows
// SUPPOSITION 1: STARTS WITH NODE 0
// SUPPOSITION 2: THE RETURNED INDEX IS INCLUDED
int get_lower_index(int node_id, int max_rows) {
	return node_id * (max_rows-2);
}




// Gets the index until which each node should get rows
// SUPPOSITION 1: STARTS WITH NODE 0
// SUPPOSITION 2: THE RETURNED INDEX IS INCLUDED
int get_upper_index(int node_id, int max_rows, int n) {

	int index = (node_id+1) * (max_rows-2) + 1;
	if (index >= n) {
		return n-1;
	}
	else {
		return index;
	}
}




// Allocate 2D matrix in the master node
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




// Allocate 2D matrix in the slaves nodes
void allocate_nodes_2Dmatrix(float ***mat, int n, int m) {

	int i;
	*mat = (float **) malloc(n * sizeof(float *));

	for (i = 0; i < n; i++) {
		(*mat)[i] = (float *)malloc(m * sizeof(float));
	}
}




// Solves as many rows as specified at the argument "n"
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
	if (done) {
		printf("Node %d: Solver converged after %d iterations\n", myrank, cnt_iter);
	}
	else {
		printf("Node %d: Solver not converged after %d iterations\n", myrank, cnt_iter);
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


	int max_rows = get_max_rows(np, n);
	int lower_index = get_lower_index(myrank, max_rows);
	int upper_index = get_upper_index(myrank, max_rows, n);

	// Calculate the proper number of rows and elements for each node
	int num_rows = get_node_rows(lower_index, upper_index);
	int num_elems = get_node_elems(num_rows, n);


	switch(communication) {
		case 0: {
			if (myrank == 0) {

				// Allocating memory for the whole matrix
				allocate_init_2Dmatrix(&a, n, n);
				int i;

				// Master sends chuncks to every other node
				for (i = 1; i < np; i++) {
					int i_lower_index = get_lower_index(i, max_rows);
					int i_upper_index = get_upper_index(i, max_rows, n);
					int i_num_rows = get_node_rows(i_lower_index, i_upper_index);
					int i_num_elems = get_node_elems(i_num_rows, n);

					MPI_Send(&a[i_lower_index], i_num_elems, MPI_FLOAT, i, MPI_ANY_TAG, MPI_COMM_WORLD);
				}
			}
			else {

				// Allocating the exact memory to the rows receiving
				allocate_nodes_2Dmatrix(&a, num_rows, n);
				int status;

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
	solver(&a, num_rows, n);


	switch(communication) {
		case 0: {
			if (myrank == 0) {
				int i, status;

				// Master sends chuncks to every other node
				for (i = 1; i < np; i++) {
					int i_lower_index = get_lower_index(i, max_rows) + 1;		// +1 to skip cortex values
					int i_upper_index = get_upper_index(i, max_rows, n) - 1;	// -1 to skip cortex values
					int i_num_rows = get_node_rows(i_lower_index, i_upper_index);
					int i_num_elems = get_node_elems(i_num_rows, n);

					// Receiving the data from the slave nodes
					MPI_Recv(&a[i_lower_index], i_num_elems, MPI_FLOAT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				}
			}
			else {
				int solved_lower_index = 1;						// Start at 1 to skip cortex values
				int solved_upper_index = num_rows-2;			// Reach num_rows-2 to skip cortex values
				int solved_rows = get_node_rows(solved_lower_index, solved_upper_index);
				int solved_elems = get_node_elems(solved_rows, n);

				// Compute num_elems sin la corteza
				MPI_Send(&a[solved_lower_index], solved_elems, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
			}

			break;
		}
		/*
		case 1: {
			// collective communication for gathering the matrix
			break;
		}
		*/
	}

	MPI_Finalize();
	return 0;
}
