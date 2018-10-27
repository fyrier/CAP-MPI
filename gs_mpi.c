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




// Gets the position from which elements are gonna be sent / received
int get_node_offset(int node_id, int n, int max_rows) {
	return node_id * n * (max_rows-2);
}




// Calculates how many elements are going to a given node
int get_node_elems(int node_offset, int n, int max_rows) {
	int node_elems = max_rows*n

	// Case in which the node receive the full set of elements
	if (node_offset + node_elems < (n*n)) {
		return node_elems
	}

	// Case of the last node, which could get less elements
	else {
		return (node_offset+node_elems) - (n*n)
	}
}




// Allocate 2D matrix in the master node
void allocate_root_matrix(float **mat, int n, int m){

	*mat = (float *) malloc(n * m * sizeof(float));
	for (int i = 0; i < (n*m); i++) {
		mat[i] = rand_float(MAX);
	}
}




// Allocate 2D matrix in the slaves nodes
void allocate_node_matrix(float **mat, int num_elems) {
	*mat = (float *) malloc(num_elems * sizeof(float));
}




// Free all memory from the allocated 2D matrices
void free_2Dmatrix(float **mat, int n) {
	free(*mat);
}




// Solves as many rows as specified at the argument "n"
void solver(float *mat, int n, int m) {

	float diff = 0, temp;
	int done = 0, cnt_iter = 0, myrank;

  	while (!done && (cnt_iter < MAX_ITER)) {
  		diff = 0;

  		for (int i = 1; i < n - 1; i++) {
  			for (int j = 1; j < m - 1; j++) {

  				int pos = i * m + j;
  				int pos_up = (i - 1) * m + j;
  				int pos_do = (i + 1) * m + j;
  				int pos_le = i * m + (j-1);
  				int pos_ri = i * m + (j+1);

  				temp = mat[pos];
				mat[pos] = 0.2 * (mat[pos] + mat[pos_le] + mat[pos_up] + mat[pos_ri] + mat[pos_do]);
				diff += abs(mat[pos] - temp);
  			}
      	}

		if (diff/n/n < TOL) {
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

	int np, myrank, n, communication, i;
	float *a;

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


	// Calculate common relevant values for each node
	int max_rows = get_max_rows(np, n);
	int offset = get_node_offset(myrank, n, max_rows);
	int elems = get_node_elems(myrank, n, max_rows);

	double tscom1 = MPI_Wtime();


	switch(communication) {
		case 0: {
			if (myrank == 0) {

				// Allocating memory for the whole matrix
				allocate_root_matrix(&a, n, n);

				// Master sends chuncks to every other node
				for (i = 1; i < np; i++) {
					int i_offset = get_node_offset(i, n, max_rows);
					int i_elems = get_node_elems(i, n, max_rows);

					MPI_Send(&a[i_offset], i_elems, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
				}
			}
			else {

				// Allocating the exact memory to the rows receiving
				allocate_node_matrix(&a, elems);
				MPI_Status status;

				// Receiving the data from the master node
				MPI_Recv(&a, elems, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}
			break;
		}

		case 1: {
			if (myrank == 0) {
				// Allocating memory for the whole matrix
				allocate_root_matrix(&a, n, n);		
			}
			else {
				// Allocating the exact memory to the rows receiving
				allocate_node_matrix(&a, elems);
			}

			// Collective communication for scattering the matrix
			// Info: https://www.mpich.org/static/docs/v3.1/www3/MPI_Scatterv.html
			//MPI_Scatterv();
			break;
		}
	}


	double tfcom1 = MPI_Wtime();
	double tsop = MPI_Wtime();

	// --------- SOLVER ---------
	solver(&a, num_rows, n);

	double tfop = MPI_Wtime();
	double tscom2 = MPI_Wtime();


	switch(communication) {
		case 0: {
			if (myrank == 0) {

				MPI_Status status;

				// Master sends chuncks to every other node
				for (i = 1; i < np; i++) {
					int i_offset = get_node_offset(i, n, max_rows) + n; 	// +n to skip cortex values
					int i_elems = get_node_elems(i, n, max_rows) - (2*n);	// -2n to skip cortex values

					// Receiving the data from the slave nodes
					MPI_Recv(&a[i_offset], i_elems, MPI_FLOAT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				}
			}
			else {
				int solved_offset = n;										// Start at n to skip cortex values
				int solved_elems = num_elems - (2*n);						// Reach num_elems-2n to skip cortex values

				// Compute num_elems sin la corteza
				MPI_Send(&a[solved_offset], solved_elems, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
			}

			break;
		}
		
		case 1: {

			// Collective communication for gathering the matrix
			// Info: http://www.mpich.org/static/docs/v3.2.1/www/www3/MPI_Gatherv.html
			//MPI_Gatherv();
			break;
		}
	}

	double tfcom2 = MPI_Wtime();

	printf("Tiempo de comunicacion: %f", (tfcom1-tscom1) + (tfcom2-tscom2));
	printf("Tiempo de operacion: %f", tfop - tsop);
	printf("Tiempo total: %f", (tfcom1-tscom1) + (tfcom2-tscom2) + (tfop-tsop));


	// Finally, free the 2D matrix memory allocated
	free_2Dmatrix(&a);
	MPI_Finalize();
	return 0;
}
