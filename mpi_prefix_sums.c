#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>

int main(void)
{
  int i, num, temp, *arr, *input;
  int psum = 0;

  MPI_Comm comm;
  int comm_sz, my_rank;

  comm = MPI_COMM_WORLD;

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(comm, &my_rank);
  MPI_Comm_size(comm, &comm_sz);

  arr = (int*)malloc(comm_sz*sizeof(int));

  if (my_rank == 0) {
      srand(getpid());
      input = (int*)malloc(comm_sz*sizeof(int));
      printf("Input is:");
      for (i = 0; i < comm_sz; i++) {
          input[i] = rand()%10;
          printf(" %d", input[i]);
      }
      printf("\n");
  }
  
  MPI_Scatter(input, 1, MPI_INT, &num, 1, MPI_INT, 0, comm);

  MPI_Gather(&psum, 1, MPI_INT, arr, 1, MPI_INT, 0, comm);

  if (my_rank == 0)
    psum = num+num;
  
  MPI_Reduce(&num, &temp, 1, MPI_INT, MPI_SUM, 0, comm);

  MPI_Allreduce(&num, &temp, 1, MPI_INT, MPI_SUM, comm);

  MPI_Gather(&psum, 1, MPI_INT, arr, 1, MPI_INT, 0, comm);

  if (my_rank == 0) {
      for (i = 0; i < comm_sz; i++)
          printf(" %d", arr[i]);
      printf("\n");
  }

  free(arr);

  MPI_Finalize();

}