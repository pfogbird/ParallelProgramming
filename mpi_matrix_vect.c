#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

int my_rank, comm_sz;
MPI_Comm comm;

void Get_size(int* m_p, int* local_m_p, int* n_p, int* local_n_p);
void Read_matrix(char title[], double local_matrix[], int m, int local_m, int n, int local_n);
void Read_vector(char title[], double local_vec[], int n, int local_n);
void Print_matrix(char title[], double local_matrix[], int m, int local_m, int n);
void Print_vector(char title[], double local_vec[], int n, int local_n);
void Mat_vect_mult(double local_matrix[], double local_vector[], double local_result[], double temp[], int m, int local_m, int n, int local_n);

int main(void) {
   double* local_matrix;
   double* local_vector;
   double* local_result;
   double* temp;
   int m, local_m, n, local_n;

   MPI_Init(NULL, NULL);
   comm = MPI_COMM_WORLD;
   MPI_Comm_size(comm, &comm_sz);
   MPI_Comm_rank(comm, &my_rank);

   Get_size(&m, &local_m, &n, &local_n);
   local_matrix = malloc(local_m*n*sizeof(double));
   local_vector = malloc(local_n*sizeof(double));
   local_result = malloc(local_m*sizeof(double));
   Read_matrix("Input Matrix", local_matrix, m, local_m, n, local_n);
   
   Print_matrix("Input Matrix", local_matrix, m, local_m, n);

   Read_vector("Input vector", local_vector, n, local_n);
   
   Print_vector("Input vector", local_vector, n, local_n);

   temp = malloc(n*sizeof(double));
   MPI_Barrier(comm);
   
   Mat_vect_mult(local_matrix, local_vector, local_result, temp, m, local_m, n, local_n);

   Print_vector("Result", local_result, m, local_m);
   
   free(local_matrix);
   free(local_vector);
   free(local_result);
   free(temp);
   MPI_Finalize();
   return 0;
}  

void Get_size(int* m_p, int* local_m_p, int* n_p, int* local_n_p) {
  
   if (my_rank == 0){

        printf("What is the order of the vector? \n");
        scanf("%d", n_p);
    }

    if (my_rank == 0){

        printf("What is the order of the vector? \n");
        scanf("%d", m_p);
    }

   MPI_Bcast(m_p, 1, MPI_INT, 0, comm);
   MPI_Bcast(n_p, 1, MPI_INT, 0, comm);

   *local_m_p = *m_p/comm_sz;
   *local_n_p = *n_p/comm_sz;
}  

void Read_matrix(char title[], double local_matrix[], int m, int local_m, int n, int local_n) {
   double* matrix = NULL;
   int i, j;

   if (my_rank == 0) {
      matrix = malloc(m*n*sizeof(double));
      
      printf("Enter the matrix %s\n", title);
      for (i = 0; i < m; i++)
         for (j = 0; j < n; j++)
            scanf("%lf", &matrix[i*n+j]);
      MPI_Scatter(matrix, local_m*n, MPI_DOUBLE, local_matrix, local_m*n, MPI_DOUBLE, 0, comm);
      free(matrix);
   } else {
      
      MPI_Scatter(matrix, local_m*n, MPI_DOUBLE, local_matrix, local_m*n, MPI_DOUBLE, 0, comm);

   }
}  

void Read_vector(char title[], double local_vec[], int n, int local_n) {
   double* vec = NULL;
   int i;

   if (my_rank == 0) {
      vec = malloc(n*sizeof(double));
      
      printf("Enter the vector %s\n", title);
      for (i = 0; i < n; i++)
         scanf("%lf", &vec[i]);
      MPI_Scatter(vec, local_n, MPI_DOUBLE, local_vec, local_n, MPI_DOUBLE, 0, comm);
      free(vec);
   } else {
      
      MPI_Scatter(vec, local_n, MPI_DOUBLE, local_vec, local_n, MPI_DOUBLE, 0, comm);
   }
}  

void Print_matrix(char title[], double local_matrix[], int m, int local_m, int n) {

   double* matrix = NULL;
   int i, j;

   if (my_rank == 0) {
      matrix = malloc(m*n*sizeof(double));
      
      MPI_Gather(local_matrix, local_m*n, MPI_DOUBLE, matrix, local_m*n, MPI_DOUBLE, 0, comm);

      printf("\nThe matrix %s\n", title);
      for (i = 0; i < m; i++) {
         for (j = 0; j < n; j++)
            printf("%f ", matrix[i*n+j]);
         printf("\n");
      }
      printf("\n");
      free(matrix);
   } else {
      
      MPI_Gather(local_matrix, local_m*n, MPI_DOUBLE, matrix, local_m*n, MPI_DOUBLE, 0, comm);
   }
}  

void Print_vector(char title[], double local_vec[], int n, int local_n) {
   double* vec = NULL;
   int i;

   if (my_rank == 0) {
      vec = malloc(n*sizeof(double));
      
      MPI_Gather(local_vec, local_n, MPI_DOUBLE, vec, local_n, MPI_DOUBLE, 0, comm);
      printf("\nThe vector %s\n", title);
      for (i = 0; i < n; i++)
         printf("%f ", vec[i]);
      printf("\n");
      free(vec);
   }  else {
      
      MPI_Gather(local_vec, local_n, MPI_DOUBLE, vec, local_n, MPI_DOUBLE, 0, comm);
   }
}  

void Mat_vect_mult(double local_matrix[], double local_vector[], double  local_result[], double temp[], int m, int local_m, int n, int local_n) {
   int i, j;

   MPI_Allgather(local_vector, local_n, MPI_DOUBLE, temp, local_n, MPI_DOUBLE, comm);

   for (i = 0; i < local_m; i++) {
      local_result[i] = 0.0;
      for (j = 0; j < n; j++)
         local_result[i] += local_matrix[i*n+j]*temp[j];
   }
}  
