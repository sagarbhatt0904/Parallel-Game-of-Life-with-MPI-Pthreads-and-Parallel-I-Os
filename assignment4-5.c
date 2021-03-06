/***************************************************************************/
/* Asssignment 4/5 *********************************************************/
/* Ajinkya Dahale, Sagar Bhatt *********************************************/
/*CSCI 6360: Parallel Computing*********************************************/
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

// #define BGQ 1 // when running BG/Q, comment out when running on kratos
#ifdef BGQ
#include<hwi/include/bqc/A2_inlines.h>
double processor_frequency = 1600000000.0;
unsigned long long start_cycles = 0;
unsigned long long end_cycles = 0;
#else
#define GetTimeBase MPI_Wtime
double start_cycles = 0.0, end_cycles = 0.0;
#endif
double time_in_secs1 = 0.0;

/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

#define ALIVE 1
#define DEAD  0
#define NTICKS 128
#define UNIV_EDGE_SIZE 16384
#define CELLS_PER_PIXEL 16
#define RAND_THRESH 0.25
#define NUM_THREADS 32
typedef unsigned short int celltype;
typedef unsigned short int pixel_t;

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

celltype** univ;
celltype** univ_new;
int mpi_myrank;
int mpi_commsize;
size_t nrows;
size_t ncols;
pthread_barrier_t barrier;
FILE *out_file = NULL;

struct datapack {
  size_t tid;
};


/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

void contour(FILE * file1, celltype** univ, int Nx, int Ny, size_t tick, char* str);
void my_swap_2d(celltype*** a, celltype*** b);
void draw_img(celltype** univ, pixel_t** img, size_t nrows, size_t ncols);
void fprint_self(FILE* out_file, celltype** univ, size_t nrows, size_t ncols);
void fprint_univ(MPI_File out_file, celltype** univ, size_t nrows, size_t ncols);
void fprint_pixels(MPI_File out_file, pixel_t** img, size_t nrows, size_t ncols);
void do_ghosting(celltype**univ, size_t nrows, size_t ncols, size_t tick);
void* play_gol_pthreads(void *data);
void play_gol(celltype** univ, celltype** univ_new, size_t nrows, size_t ncols, size_t my_start_row);

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{
  size_t i, j, thr;
  size_t nirows;
  size_t nicols;
  pthread_t my_threads[NUM_THREADS];
  struct datapack datapacks[NUM_THREADS];
  pixel_t **img;
  // MPI_File mpi_out_fh;
  MPI_File mpi_img_fh;
  // char out_filename[16];

// Example MPI startup and using CLCG4 RNG
  MPI_Init( &argc, &argv);
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
  MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);

  nrows = UNIV_EDGE_SIZE / mpi_commsize;
  ncols = UNIV_EDGE_SIZE;
  nirows = nrows / CELLS_PER_PIXEL;
  nicols = ncols / CELLS_PER_PIXEL;

  // MPI_File_open(MPI_COMM_WORLD, "out_all",
  //               MPI_MODE_CREATE | MPI_MODE_WRONLY,
  //               MPI_INFO_NULL, &mpi_out_fh);

  MPI_File_open(MPI_COMM_WORLD, "out_all_img",
                MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &mpi_img_fh);

  pthread_barrier_init(&barrier, NULL, NUM_THREADS);

// Init 16,384 RNG streams - each rank has an independent stream
  InitDefault();

  MPI_Barrier( MPI_COMM_WORLD );

  /* To reiterate: 2 ghosted rows */
  univ  = (celltype **)malloc(sizeof(celltype *) * (nrows + 2));
  univ[0] = (celltype *)malloc(sizeof(celltype) * ncols * (nrows + 2));

  for (i = 0; i < (nrows + 2); i++) {
    univ[i] = (*univ + ncols * i);
  }

  univ_new  = (celltype **)malloc(sizeof(celltype *) * (nrows + 2));
  univ_new[0] = (celltype *)malloc(sizeof(celltype) * ncols * (nrows + 2));

  for (i = 0; i < (nrows + 2); i++) {
    univ_new[i] = (*univ_new + ncols * i);
  }

  img = (pixel_t **) malloc(sizeof(pixel_t *) * (nirows + 2));
  img[0] = (pixel_t *) malloc(sizeof(pixel_t) * nicols * nirows);

  for (i = 0; i < nirows; i++) {
    img[i] = (*img + nicols * i);
  }

  /* Defining initial state */
  /*  Randomize */
  for (i = 0; i < nrows; ++i) {
    for (j = 0; j < ncols; ++j) {
       univ[i+1][j] = (GenVal(i) < 0.5) ? DEAD : ALIVE; 
    }
  }

  
  // Calculate Start time
  start_cycles = GetTimeBase();   // P2P time start

  for (thr = 0; thr < NUM_THREADS; ++thr) {
    datapacks[thr].tid = thr;
    pthread_create(&my_threads[thr], NULL,play_gol_pthreads, (void *) &datapacks[thr]);
  }

  for (thr = 0; thr < NUM_THREADS; ++thr) {
    pthread_join(my_threads[thr], NULL);
  }

  // fprint_univ(mpi_out_fh, univ, nrows, ncols);
  draw_img(univ, img, nrows, ncols);
  fprint_pixels(mpi_img_fh, img, nirows, nicols);
  
  //Calculate end time
  end_cycles = GetTimeBase();
#ifdef BGQ

  time_in_secs1 = ((double)(end_cycles - start_cycles)) / processor_frequency;
#else

  time_in_secs1 = ((double)(end_cycles - start_cycles));
#endif
  // fprint_self(out_file, univ, nrows, ncols);

// END -Perform a barrier and then leave MPI
  MPI_Barrier( MPI_COMM_WORLD );
  free(univ[0]);
  free(univ);
  // fclose(out_file);
  // MPI_File_close(&mpi_out_fh);
  MPI_File_close(&mpi_img_fh);

  if (mpi_myrank==0)
  {
    /**********************************************************
      Writing time to file for Plotting
    **********************************************************/
    FILE *ptrFileOut;
    ptrFileOut = fopen("time.dat", "a");
    fprintf( ptrFileOut, "%d %d %f \n", mpi_commsize, NUM_THREADS, time_in_secs1);
    fclose(ptrFileOut);
  }
  MPI_Finalize();
  return 0;
}

/***************************************************************************/
/* Other Functions *********************************************************/
/***************************************************************************/

void my_swap_2d(celltype*** a, celltype*** b) {
  celltype** temp = *a;
  *a = *b;
  *b = temp;
}

// Function To reduce the size of the grid to visualize in order to save space
void draw_img(celltype** univ, pixel_t** img, size_t nrows, size_t ncols)
{
  size_t nirows = nrows / CELLS_PER_PIXEL;
  size_t nicols = ncols / CELLS_PER_PIXEL;
  size_t ci, cj, ii, ij;
  for (ii = 0; ii < nirows; ++ii) {
    for (ij = 0; ij < nicols; ++ij) {
      img[ii][ij] = 0;
      for (ci = 0; ci < CELLS_PER_PIXEL; ++ci) {
        for (cj = 0; cj < CELLS_PER_PIXEL; ++cj) {
          /* Here we ghost, one more time */
          img[ii][ij] += univ[ii * CELLS_PER_PIXEL + ci + 1][ij * CELLS_PER_PIXEL + cj];
        }
      }
    }
  }

}
// Print function for debugging
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

void play_gol(celltype** univ, celltype** univ_new, size_t nrows, size_t ncols, size_t my_start_row)
{
  size_t i, j;
  celltype sum;
  for (i = 0; i < nrows; ++i) {
    for (j = 0; j < ncols; ++j) {
       if (GenVal(my_start_row + i) < RAND_THRESH) { 
        univ_new[i+1][j] = (GenVal(my_start_row + i) < 0.5) ? DEAD : ALIVE; 
        continue; 
       } 
      sum = (univ[i][(j + ncols - 1) % ncols] + univ[i][j] + univ[i][(j + 1) % ncols] +
             univ[i + 1][(j + ncols - 1) % ncols] + univ[i + 1][(j + 1) % ncols] +
             univ[i + 2][(j + ncols - 1) % ncols] + univ[i + 2][j] + univ[i + 2][(j + 1) % ncols]);
      switch (sum) {
      case 2:
        /* Do nothing to the state */
        univ_new[i + 1][j] = univ[i + 1][j];
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

void* play_gol_pthreads(void *data)
{
  struct datapack* data1 = (struct datapack *) data;
  int rc;
  size_t i;
  size_t ntrows = nrows / NUM_THREADS;
  size_t start_row = ntrows * data1->tid;
  size_t true_start_row = nrows * mpi_myrank + start_row;
  // FILE* file1 = NULL; char str[20];
  /* Thread is started. Go into ticks */
  for (i = 0; i < NTICKS; ++i) {
    rc = pthread_barrier_wait(&barrier);
    if (PTHREAD_BARRIER_SERIAL_THREAD == rc) {
      
      do_ghosting(univ, nrows, ncols, i);
    }
    // printf("rank %d, thread %ld tick %ld\n", mpi_myrank, data1->tid, i);
    rc = pthread_barrier_wait(&(barrier));

    play_gol(&(univ[start_row]), &(univ_new[start_row]), ntrows, ncols, true_start_row);

    rc = pthread_barrier_wait(&(barrier));
    if (PTHREAD_BARRIER_SERIAL_THREAD == rc) {
      my_swap_2d(&(univ), &(univ_new));
      // contour(file1,  univ, nrows, ncols, i, str);
    }
  }

  pthread_exit(NULL);
}

void do_ghosting(celltype **univ, size_t nrows, size_t ncols, size_t tick)
{
  /* int mpi_myrank, mpi_commsize; */
  /* up/down refers to direction of message */
  MPI_Request recv_req_up, recv_req_down,
              send_req_up, send_req_down;
  MPI_Status recv_status_up, recv_status_down,
             send_status_up, send_status_down;
  int up_tag = 2 * tick, down_tag = 2 * tick + 1;
  int rank_above, rank_below;
  /* Simply (mpi_myrank-1)%mpi_commsize gives -1 for rank 0.
     Even GDB can't help! */
  rank_above = (mpi_myrank + mpi_commsize - 1) % mpi_commsize;
  rank_below = (mpi_myrank + 1) % mpi_commsize;
  /*  Request ghosts */
  MPI_Irecv(univ[nrows + 1], ncols, MPI_UNSIGNED_SHORT,
            rank_below, up_tag, MPI_COMM_WORLD,
            &recv_req_up);
  MPI_Irecv(univ[0], ncols, MPI_UNSIGNED_SHORT,
            rank_above, down_tag, MPI_COMM_WORLD,
            &recv_req_down);

  /*  Send shosts */
  MPI_Isend(univ[1], ncols, MPI_UNSIGNED_SHORT,
            rank_above, up_tag, MPI_COMM_WORLD,
            &send_req_up);
  MPI_Isend(univ[nrows], ncols, MPI_UNSIGNED_SHORT,
            rank_below, down_tag, MPI_COMM_WORLD,
            &send_req_down);

  /*  Play the game */
  MPI_Wait(&recv_req_up, &recv_status_up);
  MPI_Wait(&recv_req_down, &recv_status_down);
  MPI_Wait(&send_req_up, &send_status_up);
  MPI_Wait(&send_req_down, &send_status_down);
}

// Function for visualization. (Note, this runs on 1 rank (for debugging purposes))
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
