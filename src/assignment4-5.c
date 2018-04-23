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
#define CELLS_PER_PIXEL 4
typedef unsigned short int celltype;
typedef unsigned short int pixel_t;

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

// You define these


/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// You define these

void my_swap(celltype** a, celltype** b);

void draw_img(celltype** univ, pixel_t** img, size_t nrows, size_t ncols);

void fprint_self(FILE* out_file, celltype** univ, size_t nrows, size_t ncols);

void fprint_univ(MPI_File out_file, celltype** univ, size_t nrows, size_t ncols);

void fprint_pixels(MPI_File out_file, pixel_t** img, size_t nrows, size_t ncols);

/** Play one tick of a game of life serially with given array of known number
 *  of rows and columns. Assumed that rows wrap around, but first and last
 *  rows are ghosted.
 */
void play_gol(celltype** univ, size_t nrows, size_t ncols);

/** Play one tick of the game of life using MPI */
void play_gol_mpi(celltype**univ, size_t nrows, size_t ncols, size_t tick);

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{
//    int i = 0;
  size_t nrows = 8, ncols = 16, i, j;
  size_t nirows = nrows / CELLS_PER_PIXEL;
  size_t nicols = ncols / CELLS_PER_PIXEL;
  celltype **univ;
  pixel_t **img;
  FILE *out_file = NULL;
  MPI_File mpi_out_fh;
  MPI_File mpi_img_fh;
  char out_filename[16];
    int mpi_myrank;
    int mpi_commsize;

// Example MPI startup and using CLCG4 RNG
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);

    MPI_File_open(MPI_COMM_WORLD, "out_all",
		  MPI_MODE_CREATE | MPI_MODE_WRONLY,
		  MPI_INFO_NULL, &mpi_out_fh);

    MPI_File_open(MPI_COMM_WORLD, "out_all_img",
		  MPI_MODE_CREATE | MPI_MODE_WRONLY,
		  MPI_INFO_NULL, &mpi_img_fh);

// Init 16,384 RNG streams - each rank has an independent stream
    InitDefault();

// Note, used the mpi_myrank to select which RNG stream to use.
// You must replace mpi_myrank with the right row being used.
// This just show you how to call the RNG.
    printf("Rank %d of %d has been started and a first Random Value of %lf\n",
	   mpi_myrank, mpi_commsize, GenVal(mpi_myrank));

    MPI_Barrier( MPI_COMM_WORLD );

// Insert your code
    if (1) {
      sprintf(out_filename, "out_rank_%d", mpi_myrank);
      if ((out_file = fopen(out_filename, "w")) == NULL) {
	printf("Failed to open output file: %s.\n", out_filename);
      }

      /* To reiterate: 2 ghosted rows */
      univ  = (celltype **)malloc(sizeof(celltype *) *(nrows+2));
      univ[0] = (celltype *)malloc(sizeof(celltype) * ncols *(nrows+2));

      /* Defining initial state */
      /* TODO: Randomize */
      for(i = 0; i < (nrows+2); i++) {
	univ[i] = (*univ + ncols * i);
	for (j = 0; j < ncols; ++j) {
	  univ[i][j] = DEAD;
	}
      }

      img = (pixel_t **) malloc(sizeof(pixel_t *) * (nirows+2));
      img[0] = (pixel_t *) malloc(sizeof(pixel_t) * nicols * nirows);

      for(i = 0; i < nirows; i++) {
	img[i] = (*img + nicols * i);
      }

      if (0 == mpi_myrank) {
	univ[2][2] = ALIVE;
	univ[3][3] = ALIVE;
	univ[4][1] = ALIVE;
	univ[4][2] = ALIVE;
	univ[4][3] = ALIVE;
      }

      size_t tick;
      for (tick = 0; tick < 32; ++tick) {
	fprintf(out_file, "Rank %d, tick %ld:\n", mpi_myrank, tick);
	fprint_self(out_file, univ, nrows, ncols);
	/* TODO: play one tick */
	play_gol_mpi(univ, nrows, ncols, tick);
      }

      fprint_univ(mpi_out_fh, univ, nrows, ncols);
      draw_img(univ, img, nrows, ncols);
      fprint_pixels(mpi_img_fh, img, nirows, nicols);
      fprintf(out_file, "Rank %d, tick %ld:\n", mpi_myrank, tick);
      fprint_self(out_file, univ, nrows, ncols);
    }

// END -Perform a barrier and then leave MPI
    MPI_Barrier( MPI_COMM_WORLD );
    free(univ[0]);
    free(univ);
    fclose(out_file);
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

void draw_img(celltype** univ, pixel_t** img, size_t nrows, size_t ncols)
{
  /* TODO: erase old image? */
  size_t nirows = nrows/CELLS_PER_PIXEL;
  size_t nicols = ncols/CELLS_PER_PIXEL;
  size_t ci, cj, ii, ij;
  for (ii = 0; ii < nirows; ++ii) {
    for (ij = 0; ij < nicols; ++ij) {
      img[ii][ij] = 0;
      for (ci = 0; ci < CELLS_PER_PIXEL; ++ci) {
	for (cj = 0; cj < CELLS_PER_PIXEL; ++cj) {
	  /* Here we ghost, one more time */
	  img[ii][ij] += univ[ii*CELLS_PER_PIXEL + ci + 1][ij*CELLS_PER_PIXEL + cj];
	}
      }
    }
  }

}

void fprint_self(FILE* out_file, celltype** univ, size_t nrows, size_t ncols)
{
  size_t i, j;
  fprintf(out_file, "================================\n");
  for (i = 1; i <= nrows; ++i) {
    for (j = 0; j < ncols; ++j) {
      fprintf(out_file, "%d", univ[i][j]);
    }
    fprintf(out_file, "\n");
  }
  fprintf(out_file, "================================\n");
}

void fprint_univ(MPI_File out_file, celltype** univ, size_t nrows, size_t ncols)
{
  /* size_t i, j; */
    int mpi_myrank;
    int mpi_commsize;
    MPI_Status status;
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    /* Note: univ[1] because of the ghost row */
    MPI_File_write_at(out_file, mpi_myrank * nrows * ncols * sizeof(celltype), univ[1], nrows * ncols, MPI_UNSIGNED_SHORT, &status);
}

void fprint_pixels(MPI_File out_file, pixel_t** img, size_t nrows, size_t ncols)
{
    int mpi_myrank;
    int mpi_commsize;
    MPI_Status status;
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    /* Note: univ[1] because of the ghost row */
    MPI_File_write_at(out_file, mpi_myrank * nrows * ncols * sizeof(pixel_t), img[0], nrows * ncols, MPI_UNSIGNED_SHORT, &status);
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

void play_gol_mpi(celltype**univ, size_t nrows, size_t ncols, size_t tick)
{
  int mpi_myrank, mpi_commsize;
  /* up/down refers to direction of message */
  MPI_Request recv_req_up, recv_req_down,
    send_req_up, send_req_down;
  MPI_Status recv_status_up, recv_status_down,
    send_status_up, send_status_down;
  int up_tag = 2*tick, down_tag = 2*tick+1;
  int rank_above, rank_below;

  MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
  MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);

  /* Simply (mpi_myrank-1)%mpi_commsize gives -1 for rank 0.
     Even GDB can't help! */
  rank_above = (mpi_myrank+mpi_commsize-1)%mpi_commsize;
  rank_below = (mpi_myrank+1)%mpi_commsize;
  /* TODO: Request ghosts */
  MPI_Irecv(univ[nrows+1], ncols, MPI_UNSIGNED_SHORT,
	    rank_below, up_tag, MPI_COMM_WORLD,
	    &recv_req_up);
  MPI_Irecv(univ[0], ncols, MPI_UNSIGNED_SHORT,
	    rank_above, down_tag, MPI_COMM_WORLD,
	    &recv_req_down);

  /* TODO: Send shosts */
  MPI_Isend(univ[1], ncols, MPI_UNSIGNED_SHORT,
	    rank_above, up_tag, MPI_COMM_WORLD,
	    &send_req_up);
  MPI_Isend(univ[nrows], ncols, MPI_UNSIGNED_SHORT,
	    rank_below, down_tag, MPI_COMM_WORLD,
	    &send_req_down);

  /* TODO: Play the game */
  MPI_Wait(&recv_req_up, &recv_status_up);
  MPI_Wait(&recv_req_down, &recv_status_down);
  MPI_Wait(&send_req_up, &send_status_up);
  MPI_Wait(&send_req_down, &send_status_down);
  play_gol(univ, nrows, ncols);
}
