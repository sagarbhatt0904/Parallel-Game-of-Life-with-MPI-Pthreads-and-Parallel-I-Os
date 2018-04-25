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
#define NTICKS 128
#define UNIV_EDGE_SIZE 16
#define CELLS_PER_PIXEL 4
#define RAND_THRESH 0.25
#define NUM_THREADS 1
// #define NROWS 32
// #define NCOLS 32
typedef unsigned short int celltype;
typedef unsigned short int pixel_t;
pthread_barrier_t barrier;
celltype** univ;
celltype** univ_new;
int mpi_myrank;
int my_mpichunk;
int globalCommSize;
/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/


// You define these


/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// You define these
void contour(FILE * file1, celltype** univ, int Nx, int Ny, size_t tick, char* str);
void my_swap(celltype** a, celltype** b);
void *play_gol_thread(void *ptr);
void my_swap_2d(celltype*** a, celltype*** b);

void draw_img(celltype** univ, pixel_t** img, size_t nrows, size_t ncols);

void fprint_self(FILE* out_file, celltype** univ, size_t nrows, size_t ncols);

void fprint_univ(MPI_File out_file, celltype** univ, size_t nrows, size_t ncols);

void fprint_pixels(MPI_File out_file, pixel_t** img, size_t nrows, size_t ncols);

void do_ghosting(celltype** univ, size_t nrows, size_t ncols, size_t tick);

void* play_gol_pthreads(void *data);

/** Play one tick of a game of life serially with given array of known number
 *  of rows and columns. Assumed that rows wrap around, but first and last
 *  rows are ghosted.
 */
void play_gol(celltype** univ, size_t nrows, size_t ncols);
void play_gol_2(celltype** univ, celltype** univ_new, size_t nrows, size_t ncols, size_t my_start_row);

/** Play one tick of the game of life using MPI */
void play_gol_mpi(celltype** univ, celltype** univ_new, size_t nrows, size_t ncols, size_t tick);

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{
  MPI_Init( &argc, &argv);
//    int i = 0;
  size_t nrows, ncols, i, j;
  size_t nirows;
  size_t nicols;
  size_t tick;
  // celltype ** univ, ** univ_new;
  // struct datapack data1;
  pixel_t **img;
  FILE *out_file = NULL;
  MPI_File mpi_out_fh;
  MPI_File mpi_img_fh;
  char out_filename[16];
  
  

// Example MPI startup and using CLCG4 RNG
  
  MPI_Comm_size( MPI_COMM_WORLD, &globalCommSize);
  MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);

  nrows = UNIV_EDGE_SIZE / globalCommSize;
  ncols = UNIV_EDGE_SIZE;
  nirows = nrows / CELLS_PER_PIXEL;
  nicols = ncols / CELLS_PER_PIXEL;
  my_mpichunk = UNIV_EDGE_SIZE / globalCommSize;
  MPI_File_open(MPI_COMM_WORLD, "out_all",
                MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &mpi_out_fh);

  MPI_File_open(MPI_COMM_WORLD, "out_all_img",
                MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &mpi_img_fh);
  int rc = pthread_barrier_init(&barrier, NULL, NUM_THREADS);
  if (rc != 0)
  {
    printf("Failed to init barrier\n");
    exit(-1);
  }
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
    univ  = (celltype **)malloc(sizeof(celltype *) * (nrows + 2));
    univ[0] = (celltype *)malloc(sizeof(celltype) * ncols * (nrows + 2));

    for (i = 0; i < (nrows + 2); i++) {
      univ[i] = (* univ + ncols * i);
    }

    univ_new  = (celltype **)malloc(sizeof(celltype *) * (nrows + 2));
    univ_new[0] = (celltype *)malloc(sizeof(celltype) * ncols * (nrows + 2));

    for (i = 0; i < (nrows + 2); i++) {
      univ_new[i] = (* univ_new + ncols * i);
    }

    img = (pixel_t **) malloc(sizeof(pixel_t *) * (nirows + 2));
    img[0] = (pixel_t *) malloc(sizeof(pixel_t) * nicols * nirows);

    for (i = 0; i < nirows; i++) {
      img[i] = (*img + nicols * i);
    }

    /* Defining initial state */
    /* TODO: Randomize */
    for (i = 0; i < nrows; ++i) {
      for (j = 0; j < ncols; ++j) {
        univ[i + 1][j] = DEAD; //(GenVal(i) < 0.5) ? DEAD : ALIVE;
      }
    }
    univ[2][2] = ALIVE;
    univ[3][3] = ALIVE;
    univ[4][1] = ALIVE;
    univ[4][2] = ALIVE;
    univ[4][3] = ALIVE;
    
    int tid[NUM_THREADS];
    pthread_t thread[NUM_THREADS];
    for (int ithreads = 0; ithreads < NUM_THREADS; ++ithreads) {
      

      /* TODO: play one tick */

      // do_ghosting( univ, nrows, ncols, tick);
      // play_gol_2( univ,  univ_new, nrows, ncols, mpi_myrank * nrows);
      tid[ithreads] = ithreads;
      printf("In main: creating thread %d\n", ithreads);

      rc = pthread_create(&thread[ithreads], NULL, &play_gol_thread, (void*) &tid[ithreads]);
      if (rc)
      {
        printf("ERROR; return code from pthread_create() is %d\n", rc);
        exit(-1);
      }


      // contour(file1,  univ, nrows, ncols, tick, str);
    }
    // my_swap_2d(& univ, & univ_new);




//
    // fprint_univ(mpi_out_fh,  univ, nrows, ncols);
    // draw_img( univ, img, nrows, ncols);
    // fprint_pixels(mpi_img_fh, img, nirows, nicols);
    // fprintf(out_file, "Rank %d, tick %ld:\n", mpi_myrank, tick);
    // fprint_self(out_file,  univ, nrows, ncols);
  }

// END -Perform a barrier and then leave MPI
  MPI_Barrier( MPI_COMM_WORLD );
  free( univ[0]);
  free( univ);
  fclose(out_file);
  MPI_File_close(&mpi_out_fh);
  MPI_File_close(&mpi_img_fh);
  MPI_Barrier( MPI_COMM_WORLD );
  printf("just before final\n");
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
  size_t nirows = nrows / CELLS_PER_PIXEL;
  size_t nicols = ncols / CELLS_PER_PIXEL;
  size_t ci, cj, ii, ij;
  for (ii = 0; ii < nirows; ++ii) {
    for (ij = 0; ij < nicols; ++ij) {
      img[ii][ij] = 0;
      for (ci = 0; ci < CELLS_PER_PIXEL; ++ci) {
        for (cj = 0; cj < CELLS_PER_PIXEL; ++cj) {
          /* Here we ghost, one more time */
          img[ii][ij] +=  univ[ii * CELLS_PER_PIXEL + ci + 1][ij * CELLS_PER_PIXEL + cj];
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
      fprintf(out_file, "%d",  univ[i][j]);
    }
    fprintf(out_file, "\n");
  }
  fprintf(out_file, "================================\n");
}

void fprint_univ(MPI_File out_file, celltype** univ, size_t nrows, size_t ncols)
{
  /* size_t i, j; */
  
  
  MPI_Status status;
  //MPI_Comm_size( MPI_COMM_WORLD, &globalCommSize);
  
  /* Note:  univ[1] because of the ghost row */
  MPI_File_write_at(out_file, mpi_myrank * nrows * ncols * sizeof(celltype),  univ[1], nrows * ncols, MPI_UNSIGNED_SHORT, &status);
}

void fprint_pixels(MPI_File out_file, pixel_t** img, size_t nrows, size_t ncols)
{
  
  
  MPI_Status status;
  //MPI_Comm_size( MPI_COMM_WORLD, &globalCommSize);
  
  /* Note:  univ[1] because of the ghost row */
  MPI_File_write_at(out_file, mpi_myrank * nrows * ncols * sizeof(pixel_t), img[0], nrows * ncols, MPI_UNSIGNED_SHORT, &status);
}

void play_gol_2(celltype** univ, celltype** univ_new, size_t nrows, size_t ncols, size_t my_start_row)
{
  size_t i, j;
  celltype sum;
  for (i = my_start_row; i < nrows; ++i) {
    for (j = 0; j < ncols; ++j) {
      if (GenVal(my_start_row + i) < RAND_THRESH) {
        univ_new[i + 1][j] = (GenVal(my_start_row + i) < 0.5) ? DEAD : ALIVE;
      }
      sum = ( univ[i][(j + ncols - 1) % ncols] +  univ[i][j] +  univ[i][(j + 1) % ncols] +
              univ[i + 1][(j + ncols - 1) % ncols] +  univ[i + 1][(j + 1) % ncols] +
              univ[i + 2][(j + ncols - 1) % ncols] +  univ[i + 2][j] +  univ[i + 2][(j + 1) % ncols]);
      switch (sum) {
      case 2:
        /* Do nothing to the state */
        univ_new[i + 1][j] =  univ[i + 1][j];
        break;
      case 3:
        /* Alive, either by reproduction or by survival */
        univ_new[i + 1][j] = ALIVE;
        break;
      default:
        /* Dead, either by underpopulation (0--1 nbs) or overpopulation
           (4--8 nbs) */
        univ_new[i + 1][j] = DEAD;
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
    lastrow[j] =  univ[0][j];
  }

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; ++j) {
      rowbackup[j] =  univ[i + 1][j];
    }

    for (j = 0; j < ncols; ++j) {
      sum = (lastrow[(j - 1) % ncols] + lastrow[j] + lastrow[(j + 1) % ncols] +
             rowbackup[(j - 1) % ncols] + rowbackup[(j + 1) % ncols] +
             univ[i + 2][(j - 1) % ncols] +  univ[i + 2][j] +  univ[i + 2][(j + 1) % ncols]);
      switch (sum) {
      case 2:
        /* Do nothing to the state */
        break;
      case 3:
        /* Alive, either by reproduction or by survival */
        univ[i + 1][j] = ALIVE;
        break;
      default:
        /* Dead, either by underpopulation (0--1 nbs) or overpopulation
           (4--8 nbs) */
        univ[i + 1][j] = DEAD;
        break;
      }
    }

    my_swap(&lastrow, &rowbackup);
  }

  free(rowbackup);
  free(lastrow);
}



void do_ghosting(celltype ** univ, size_t nrows, size_t ncols, size_t tick)
{
  
  /* up/down refers to direction of message */
  MPI_Request recv_req_up, recv_req_down,
              send_req_up, send_req_down;
  MPI_Status recv_status_up, recv_status_down,
             send_status_up, send_status_down;
  int up_tag = 2 * tick, down_tag = 2 * tick + 1;
  int rank_above, rank_below;

  //MPI_Comm_size( MPI_COMM_WORLD, &globalCommSize);
  

  /* Simply (mpi_myrank-1)%globalCommSize gives -1 for rank 0.
     Even GDB can't help! */
  
  rank_above = (mpi_myrank + globalCommSize - 1) % globalCommSize;
  rank_below = (mpi_myrank + 1) % globalCommSize;
  /* TODO: Request ghosts */
  MPI_Irecv( univ[nrows + 1], ncols, MPI_UNSIGNED_SHORT,
             rank_below, up_tag, MPI_COMM_WORLD,
             &recv_req_up);
  MPI_Irecv( univ[0], ncols, MPI_UNSIGNED_SHORT,
             rank_above, down_tag, MPI_COMM_WORLD,
             &recv_req_down);

  /* TODO: Send shosts */
  MPI_Isend( univ[1], ncols, MPI_UNSIGNED_SHORT,
             rank_above, up_tag, MPI_COMM_WORLD,
             &send_req_up);
  MPI_Isend( univ[nrows], ncols, MPI_UNSIGNED_SHORT,
             rank_below, down_tag, MPI_COMM_WORLD,
             &send_req_down);

  /* TODO: Play the game */
  MPI_Wait(&recv_req_up, &recv_status_up);
  MPI_Wait(&recv_req_down, &recv_status_down);
  MPI_Wait(&send_req_up, &send_status_up);
  MPI_Wait(&send_req_down, &send_status_down);
}

void play_gol_mpi(celltype** univ, celltype** univ_new, size_t nrows, size_t ncols, size_t tick)
{
  
  /* up/down refers to direction of message */
  MPI_Request recv_req_up, recv_req_down,
              send_req_up, send_req_down;
  MPI_Status recv_status_up, recv_status_down,
             send_status_up, send_status_down;
  int up_tag = 2 * tick, down_tag = 2 * tick + 1;
  int rank_above, rank_below;

  //MPI_Comm_size( MPI_COMM_WORLD, &globalCommSize);
  

  /* Simply (mpi_myrank-1)%globalCommSize gives -1 for rank 0.
     Even GDB can't help! */
  rank_above = (mpi_myrank + globalCommSize - 1) % globalCommSize;
  rank_below = (mpi_myrank + 1) % globalCommSize;
  /* TODO: Request ghosts */
  MPI_Irecv( univ[nrows + 1], ncols, MPI_UNSIGNED_SHORT,
             rank_below, up_tag, MPI_COMM_WORLD,
             &recv_req_up);
  MPI_Irecv( univ[0], ncols, MPI_UNSIGNED_SHORT,
             rank_above, down_tag, MPI_COMM_WORLD,
             &recv_req_down);

  /* TODO: Send shosts */
  MPI_Isend( univ[1], ncols, MPI_UNSIGNED_SHORT,
             rank_above, up_tag, MPI_COMM_WORLD,
             &send_req_up);
  MPI_Isend( univ[nrows], ncols, MPI_UNSIGNED_SHORT,
             rank_below, down_tag, MPI_COMM_WORLD,
             &send_req_down);

  /* TODO: Play the game */
  MPI_Wait(&recv_req_up, &recv_status_up);
  MPI_Wait(&recv_req_down, &recv_status_down);
  MPI_Wait(&send_req_up, &send_status_up);
  MPI_Wait(&send_req_down, &send_status_down);
  play_gol_2( univ,  univ_new, nrows, ncols, mpi_myrank * nrows);
}

void *play_gol_thread(void *ptr)
{

  int tid = *((int*)ptr);
  FILE* file1 = NULL; char str[20];
  int rc;
  int thread_chunk = my_mpichunk / NUM_THREADS;
  int start = tid * thread_chunk;
  int end = (tid + 1) * thread_chunk;
  size_t ncols = UNIV_EDGE_SIZE;

  if (tid == 0)
  {
    start++;
  }
  if (tid == NUM_THREADS - 1)
  {
    end = my_mpichunk;
  }

  celltype tmp;
  size_t tick;
  for (tick = 0; tick < 128; ++tick)
  {
    rc = pthread_barrier_wait (&barrier);
    if (rc == PTHREAD_BARRIER_SERIAL_THREAD)
    {

      do_ghosting( univ, my_mpichunk, ncols, tick);
    }
    rc = pthread_barrier_wait (&barrier);

    // play_gol_2(univ, univ_new, thread_chunk, ncols, start);
play_gol_mpi( univ, univ_new, thread_chunk, ncols, tick);

    rc = pthread_barrier_wait (&barrier);

    if (rc == PTHREAD_BARRIER_SERIAL_THREAD)
    {
printf("ran GOL %d\n", tick);
      // tmp =  univ;
      //  univ =   univ_new;
      //   univ_new = tmp;
      for (int i = 0; i < my_mpichunk; ++i)
      {
        for (int j = 0; j < ncols; ++j)
        {
          
          tmp =  univ[i][j];
          univ[i][j] =   univ_new[i][j];
          univ_new[i][j] = tmp;
        }
      }
      // contour(file1,  univ, my_mpichunk, ncols, tick + 1, str);
    }
    else if (rc != 0)
    {
      printf("%d\n", tid);
      printf("Failed to create barrier\n");
      exit(-1);
    }
   
    rc = pthread_barrier_wait(&barrier);
    if (rc != 0 && rc  != PTHREAD_BARRIER_SERIAL_THREAD)
    {
      printf("%d\n", tid);
      printf("Failed to create barrier\n");
      exit(-1);
    }

  }


  return 0;
}


void contour(FILE * file1, celltype** univ, int Nx, int Ny, size_t tick, char* str)
{
  sprintf(str, "out_%ld.dat", tick);
  file1 = fopen(str, "w");
  fprintf(file1, "TITLE= GOL\n");
  fprintf(file1, "VARIABLES=X,Y,A\n");
  fprintf(file1, "ZONE I=%d, J=%d, F=POINT\n", Nx, Ny);

  for (int j = 0; j < Ny; j++)
  {
    for (int i = 0; i < Nx; i++)
    {

      fprintf(file1, "%d %d %d \n", i, j,  univ[i][j]);
    }
  }
  fclose(file1);
}