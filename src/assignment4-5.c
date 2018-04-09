/***************************************************************************/
/* Template for Asssignment 4/5 ********************************************/
/* Team Names Here              **(*****************************************/
/***************************************************************************/

/***************************************************************************/
/* Includes ****************************************************************/
/***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>

#include<clcg4.h>

#include<mpi.h>


/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

#define ALIVE 1
#define DEAD  0
typedef unsigned short int celltype;

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

// You define these


/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// You define these

void my_swap(celltype** a, celltype** b);

/** Play one tick of a game of life serially with given array of known number
 *  of rows and columns. Assumed that rows wrap around, but first and last
 *  rows are ghosted.
 */
void play_gol(celltype** univ, size_t nrows, size_t ncols);

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{
//    int i = 0;
  celltype *a, *b;
  a = (celltype *) malloc(3 * sizeof(celltype));
  b = (celltype *) malloc(3 * sizeof(celltype));
  a[0] = 1; a[1] = 2; a[2] = 3;
  b[0] = 4; b[1] = 5; b[2] = 6;
    int mpi_myrank;
    int mpi_commsize;
// Example MPI startup and using CLCG4 RNG
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    
// Init 16,384 RNG streams - each rank has an independent stream
    InitDefault();
    
// Note, used the mpi_myrank to select which RNG stream to use.
// You must replace mpi_myrank with the right row being used.
// This just show you how to call the RNG.    
    printf("Rank %d of %d has been started and a first Random Value of %lf\n", 
	   mpi_myrank, mpi_commsize, GenVal(mpi_myrank));
    
    MPI_Barrier( MPI_COMM_WORLD );
    
// Insert your code
    if (0 == mpi_myrank) {
      printf("preswap, a is {%d, %d, %d}\n", a[0], a[1], a[2]);
      printf("preswap, b is {%d, %d, %d}\n", b[0], b[1], b[2]);
      my_swap(&a, &b);
      printf("swapped, a is {%d, %d, %d}\n", a[0], a[1], a[2]);
      printf("swapped, b is {%d, %d, %d}\n", b[0], b[1], b[2]);
    } 

// END -Perform a barrier and then leave MPI
    MPI_Barrier( MPI_COMM_WORLD );
    free(a);
    free(b);
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/

void my_swap(celltype** a, celltype** b) {
  celltype* temp = *a;
  *a = *b;
  *b = temp;
}

void play_gol(celltype** univ, size_t nrows, size_t ncols)
{
  celltype *rowbackup, *lastrow, sum;
  rowbackup = (celltype *) malloc(ncols * sizeof(celltype));
  lastrow = (celltype *) malloc(ncols * sizeof(celltype));
  size_t i, j;
  for (j = 0; j < ncols; ++j) {
    lastrow[j] = univ[0][j];
  }

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; ++j) {
      rowbackup[j] = univ[i+1][j];
    }

    for (j = 0; j < ncols; ++j) {
      sum = (lastrow[(j-1)%ncols] + lastrow[j] + lastrow[(j+1)%ncols] +
	     rowbackup[(j-1)%ncols] + rowbackup[(j+1)%ncols] +
	     univ[i+2][(j-1)%ncols] + univ[i+2][j] + univ[i+2][(j+1)%ncols]);
      switch (sum) {
      case 2:
        /* Do nothing to the state */
	break;
      case 3:
	/* Alive, either by reproduction or by survival */
	univ[i+1][j] = ALIVE;
	break;
      default:
	/* Dead, either by underpopulation (0--1 nbs) or overpopulation 
	   (4--8 nbs) */
	univ[i+1][j] = DEAD;
	break;
      }
    }

    my_swap(&lastrow, &rowbackup);
  }

  free(rowbackup);
  free(lastrow);
}
