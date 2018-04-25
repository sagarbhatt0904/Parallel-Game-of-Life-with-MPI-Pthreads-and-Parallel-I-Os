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
#include<pthread.h>


/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

#define ALIVE 1
#define DEAD  0
#define NTICKS 32
#define UNIV_EDGE_SIZE 32
#define CELLS_PER_PIXEL 4
#define RAND_THRESH 0.25
#define NUM_THREADS 2
typedef unsigned short int celltype;
typedef unsigned short int pixel_t;

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

struct datapack {
  celltype** univ;
  celltype** univ_new;
  pthread_barrier_t barrier;
  size_t tid;
  size_t nrows;
  size_t ncols;
  int mpi_myrank;
};
// You define these


/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// You define these

void my_swap(celltype** a, celltype** b);

void my_swap_2d(celltype*** a, celltype*** b);

void draw_img(celltype** univ, pixel_t** img, size_t nrows, size_t ncols);

void fprint_self(FILE* out_file, celltype** univ, size_t nrows, size_t ncols);

void fprint_univ(MPI_File out_file, celltype** univ, size_t nrows, size_t ncols);

void fprint_pixels(MPI_File out_file, pixel_t** img, size_t nrows, size_t ncols);

void do_ghosting(celltype**univ, size_t nrows, size_t ncols, size_t tick);

void* play_gol_pthreads(void *data);

/** Play one tick of a game of life serially with given array of known number
 *  of rows and columns. Assumed that rows wrap around, but first and last
 *  rows are ghosted.
 */
void play_gol(celltype** univ, size_t nrows, size_t ncols);
void play_gol_2(celltype** univ, celltype** univ_new, size_t nrows, size_t ncols, size_t my_start_row);

/** Play one tick of the game of life using MPI */
void play_gol_mpi(celltype**univ, celltype** univ_new, size_t nrows, size_t ncols, size_t tick);

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{
//    int i = 0;
  size_t nrows, ncols, i, j;
  size_t nirows;
  size_t nicols;
  size_t tick;
  celltype **univ, **univ_new;
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

    nrows = UNIV_EDGE_SIZE / mpi_commsize;
    ncols = UNIV_EDGE_SIZE;
    nirows = nrows / CELLS_PER_PIXEL;
    nicols = ncols / CELLS_PER_PIXEL;

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

      for(i = 0; i < (nrows+2); i++) {
	univ[i] = (*univ + ncols * i);
      }

      univ_new  = (celltype **)malloc(sizeof(celltype *) *(nrows+2));
      univ_new[0] = (celltype *)malloc(sizeof(celltype) * ncols *(nrows+2));

      for(i = 0; i < (nrows+2); i++) {
	univ_new[i] = (*univ_new + ncols * i);
      }

      img = (pixel_t **) malloc(sizeof(pixel_t *) * (nirows+2));
      img[0] = (pixel_t *) malloc(sizeof(pixel_t) * nicols * nirows);

      for(i = 0; i < nirows; i++) {
	img[i] = (*img + nicols * i);
      }

      /* Defining initial state */
      /* TODO: Randomize */
      for (i = 0; i < nrows; ++i) {
	for (j = 0; j < ncols; ++j) {
	  univ[i+1][j] = (GenVal(i) < 0.5) ? DEAD : ALIVE;
	}
      }

      for (tick = 0; tick < NTICKS; ++tick) {
	fprintf(out_file, "Rank %d, tick %ld:\n", mpi_myrank, tick);
	fprint_self(out_file, univ, nrows, ncols);
	/* TODO: play one tick */
	do_ghosting(univ, nrows, ncols, tick);
	play_gol_2(univ, univ_new, nrows, ncols, mpi_myrank * nrows);
	my_swap_2d(&univ, &univ_new);
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
    MPI_File_close(&mpi_out_fh);
    MPI_File_close(&mpi_img_fh);
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

void my_swap_2d(celltype*** a, celltype*** b) {
  celltype** temp = *a;
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

void play_gol_2(celltype** univ, celltype** univ_new, size_t nrows, size_t ncols, size_t my_start_row)
{
  size_t i, j;
  celltype sum;
  for (i = 0; i < nrows; ++i) {
    for (j = 0; j < ncols; ++j) {
      if (GenVal(my_start_row + i) < RAND_THRESH) {
	univ_new[i+1][j] = (GenVal(my_start_row + i) < 0.5) ? DEAD : ALIVE;
      }
      sum = (univ[i][(j+ncols-1)%ncols] + univ[i][j] + univ[i][(j+1)%ncols] +
	     univ[i + 1][(j+ncols-1) % ncols] + univ[i + 1][(j+1) % ncols] +
	     univ[i + 2][(j+ncols-1) % ncols] + univ[i + 2][j] + univ[i + 2][(j+1) % ncols]);
      switch (sum) {
      case 2:
        /* Do nothing to the state */
	univ_new[i+1][j] = univ[i+1][j];
	break;
      case 3:
	/* Alive, either by reproduction or by survival */
	univ_new[i+1][j] = ALIVE;
	break;
      default:
	/* Dead, either by underpopulation (0--1 nbs) or overpopulation
	   (4--8 nbs) */
	univ_new[i+1][j] = DEAD;
	break;
      }
    }
  }
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

void* play_gol_pthreads(void *data)
{
  struct datapack* data1 = (struct datapack *) data;
  size_t i;
  size_t ntrows = data1->nrows / NUM_THREADS;
  size_t start_row = ntrows * data1->tid;
  size_t true_start_row = data1->nrows * data1->mpi_myrank + start_row;
  /* Thread is started. Go into ticks */
  for (i = 0; i < NTICKS; ++i) {
    if (0 == data1->tid) {
      /* TODO: Do all the ghosting broughaha */
      do_ghosting(data1->univ, data1->nrows, data1->ncols, i);
    }
    pthread_barrier_wait(&(data1->barrier));
    play_gol_2(&(data1->univ[start_row]), &(data1->univ_new[start_row]), ntrows, data1->ncols, true_start_row);
    pthread_barrier_wait(&(data1->barrier));
    if (0 == data1->tid) {
      my_swap_2d(&(data1->univ), &(data1->univ_new));
    }
  }

  return NULL;
}

void do_ghosting(celltype **univ, size_t nrows, size_t ncols, size_t tick)
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
}

void play_gol_mpi(celltype**univ, celltype** univ_new, size_t nrows, size_t ncols, size_t tick)
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
  play_gol_2(univ, univ_new, nrows, ncols, mpi_myrank*nrows);
}
