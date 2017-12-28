#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
 
int main(void){

    int local_number_in_circle = 10000000;
    int my_rank;                       
    double x,y;                    
    int i;
    int count=0;                
    double z;                       
    double pi;                      
    int number_of_tosses;                   
    int number_in_circle;                   
    srand(time(0));
 
    MPI_Init(NULL, NULL);                 
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);           
    if(my_rank != 0){
    
        for (i=0; i<local_number_in_circle; ++i){                  
        
            x = (double)rand()/RAND_MAX;          
            y = (double)rand()/RAND_MAX;          
            z = sqrt((x*x)+(y*y));              
            if (z<=1){
            
                ++count;                
            }
        }
    }
 
    MPI_Reduce(&count, &number_of_tosses, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&local_number_in_circle, &number_in_circle, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    number_in_circle -= local_number_in_circle;                  
 
    if (my_rank == 0){                      
    
        pi = ((double)number_of_tosses/(double)number_in_circle)*4.0;               
        printf("Pi: %f\n%i\n%d\n", pi, number_of_tosses, number_in_circle);
        
    }
 
    MPI_Finalize();                     
    return 0;
}
