#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
	

void getInput(int* n, double* a, double* b, int my_rank, int comm_sz);
double Trap(double left_endpt, double right_endpt, int trap_count, double base_length);
double f(double x);
	
int main(void){
	int my_rank, comm_sz, n, local_n;
	double a, b, h, local_a, local_b;
	double local_int, total_int;
	int source;
	
	MPI_Init(NULL, NULL); 
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	
	getInput(&n, &a, &b, my_rank, comm_sz); 
	
	h = (b-a)/n; 
	local_n = n/comm_sz; 
	
	local_a = a + my_rank*local_n*h; 
	local_b = local_a + local_n * h; 
	local_int = Trap(local_a, local_b, local_n, h); 
	
	
	if(my_rank == 0){

		total_int = local_int;
		for (source = 1; source < comm_sz; source ++){
			MPI_Recv(&local_int, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
			total_int += local_int;
		}
	}
	else{
		MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); 
	}
	
	if(my_rank == 0){ 
		printf("With n = %d trapezoids, out estimate \n", n);
		printf("of the integral from %f to %f = %.15e\n",a,b,total_int);
	}
	
	MPI_Finalize(); 
	return 0;
}
	
void getInput(int* n, double* a, double* b, int my_rank, int comm_sz){
	
	int dest;
	
	if(my_rank == 0){ 
		printf("Print n, a, and b\n");
		scanf("%d %lf %lf", n, a, b);
	
		for(dest = 1; dest < comm_sz; dest++){ 
			MPI_Send(n, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
			MPI_Send(a, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
			MPI_Send(b, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
		}
	}
	else{
		MPI_Recv(n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
		MPI_Recv(a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(b, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}
	
	
	
double Trap(double left_endpt, double right_endpt, int trap_count, double base_length){
	int i;
	double estimate,x;
	
	estimate = (f(left_endpt) + f(right_endpt))/2.0; 
	for (i = 1; i <= trap_count-1; i ++){
		x = left_endpt + i*base_length;
		estimate += f(x);
	}
	
	estimate = estimate*base_length;
	
	return estimate;
}
	
double f(double x){
	return 2 * x + 1;
}
