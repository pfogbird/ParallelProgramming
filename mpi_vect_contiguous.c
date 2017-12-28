#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

void Read_n(int* n_p, int* local_n_p, int my_rank, int comm_sz,
            MPI_Comm comm);

void Read_vector(double local_vec1[], double local_vec2[],
               int local_n, MPI_Datatype type, int my_rank, int comm_sz, MPI_Comm comm);


void Print_vector(double local_vec[], int local_n, int n, char title[],
                  int my_rank, MPI_Comm comm);

void Add_vector(double local_vec1[], double local_vec2[], double local_result[], int local_n);

int main(void) {

    int n, local_n;
    double *local_vec1, *local_vec2;
    double *local_add_vec;
    int comm_sz, my_rank;
    MPI_Comm comm;

    MPI_Init(NULL, NULL);

    comm = MPI_COMM_WORLD;

    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);

    MPI_Datatype type;

    Read_n(&n, &local_n, my_rank, comm_sz, comm);

    local_vec1 = malloc(local_n*sizeof(double));
    local_vec2 = malloc(local_n*sizeof(double));
    local_add_vec = malloc(local_n*sizeof(double));

    MPI_Type_contiguous(local_n, MPI_INT, &type);
    MPI_Type_commit(&type);

    Read_vector(local_vec1, local_vec2, local_n, type, my_rank, comm_sz, comm);

    
    if(my_rank == 0)
        printf("\n\n input data \n");
    
    Print_vector(local_vec1, local_n, n, "first vector is", my_rank, comm);
    Print_vector(local_vec2, local_n, n, "second vector is", my_rank, comm);

    Add_vector(local_vec1, local_vec2, local_add_vec, local_n);

    Print_vector(local_add_vec, local_n, n, "The added vector  of x + y is", my_rank, comm);

    free(local_add_vec);
    free(local_vec2);
    free(local_vec1);
    MPI_Finalize();
    return 0;

}

void Read_n(int* n_p, int* local_n_p, int my_rank, int comm_sz, MPI_Comm comm) {

    if (my_rank == 0){

        printf("What is the order of the vector? \n");
        scanf("%d", n_p);
    }

    MPI_Bcast(n_p, 1, MPI_INT, 0, MPI_COMM_WORLD);
    *local_n_p = *n_p / comm_sz;
}  

void Read_vector(double local_vec1[], double local_vec2[], int local_n, MPI_Datatype type, int my_rank, int comm_sz, MPI_Comm comm) {
    double* a = NULL;
    int i;
    
    if (my_rank == 0){
        a = malloc(local_n * comm_sz * sizeof(double));
        printf("Enter the first vector\n");
        for (i = 0; i < local_n * comm_sz; i++)
            scanf("%lf", &a[i]);

        MPI_Scatter(a, local_n, MPI_DOUBLE, local_vec1, local_n, MPI_DOUBLE, 0, comm);


        printf("Enter the second vector\n");
        for (i = 0; i < local_n * comm_sz; i++)
            scanf("%lf", &a[i]);

        MPI_Scatter(a, local_n, MPI_DOUBLE, local_vec2, local_n, MPI_DOUBLE, 0, comm);

        free(a);

    } else {
        MPI_Scatter(a, local_n, MPI_DOUBLE, local_vec1, local_n, MPI_DOUBLE, 0, comm);
        MPI_Scatter(a, local_n, MPI_DOUBLE, local_vec2, local_n, MPI_DOUBLE, 0, comm);

    }
}  

void Print_vector(double local_vec[], int local_n, int n, char title[],
                  int my_rank, MPI_Comm comm) {
    double* a = NULL;
    int i;

    if (my_rank == 0) {
        a = malloc(n * sizeof(double));
        MPI_Gather(local_vec, local_n, MPI_DOUBLE, a, local_n,
                   MPI_DOUBLE, 0, comm);
        printf("%s \n", title);
        for (i = 0; i < n; i++)
            printf("%.2f ", a[i]);
        printf("\n");
        free(a);
    } else {
        MPI_Gather(local_vec, local_n, MPI_DOUBLE, a, local_n,
                   MPI_DOUBLE, 0, comm);
    }

} 

void Add_vector(double local_vec1[], double local_vec2[], double local_result[], int local_n) {

    for (int i = 0; i < local_n; i++) {

        local_result[i] = local_vec1[i] + local_vec2[i];

    }

}
