#include<stdio.h>
#include<stdint.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define MULTIPLIER 1664525
#define INCREMENT 1013904223
#define MODULUS 4294967296
#define BILLION 1000000000L
const double pi = 3.141592;

 

/*
*                      Φ = cos(2πy)
*       (-1,1) +--------------------+ (1,1)
*              |                    |
*              |         y          |
*              |         ^          |
* Φ = cos(2πy) |         |          | Φ = cos(2πy)
*              |         ---> x     |
*              |                    |
*              |                    |
*      (-1,-1) +--------------------+ (1,-1)
*                      Φ = cos(2πy)
*/

void saveToCSV(int n, double (*data)[n]) {
      FILE *file = fopen("data2.csv", "w");
    
    if (file == NULL) {
        printf("Error opening file.\n");
        return;
    }

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            fprintf(file, "%lf", data[i][j]);
            if (j != n - 1) {
                fprintf(file, ",");
            }
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

double f(double x, double y){
    return 0;
    return -1*((x*x)+(y*y));
}

double g(double x, double y){
    if(x==-1.0){
        return cos(2*pi*y);
        // return 0.0;
    }
    if(x==1.0){
        return cos(2*pi*y); 
        // return 0.0;
    }
    if(y==-1.0){
        return cos(2*pi*x);
        // return 0.0;
    }
    return cos(2*pi*x); 
        // return 0.0;

}

void boundary(int n, double (*a)[n], double delta){
    
    for (int i = 0; i < n; i++)
    { 
        a[0][i] = g(-1, (i*delta)-1);
        a[n-1][i] = g(1, (i*delta)-1);
    }
    for (int i = 0; i < n; i++)
    { 
        a[i][0] = g((i*delta)-1, -1);
        a[i][n-1] = g((i*delta)-1, 1);
    }
      
}    
double generateRandom(unsigned int *seed) {
    *seed = (*seed * MULTIPLIER + INCREMENT) % MODULUS;
    return (double)(*seed) / MODULUS; // Convert to double and scale to [0,1)
}

void random_dance(int n, double (*a)[n], double delta, int N, unsigned int seed){
    #pragma acc parallel loop num_gangs(10) copy(a[0:n ][0:n])
    for (int i = 1; i < (n-1); i++)
    {
        for (int j = 1; j < (n-1); j++)
        { 
            a[i][j] = 0.0; 
            seed += (i*68)+(j*46); 
            for (int k = 0; k < N; k++) 
            {
                int i1 = i, j1 = j; 
                while(1){

                    a[i][j] += f((i1*delta)-1, (j1*delta)-1)*delta*delta/4;

                    double rand_num = generateRandom(&seed);
                    if(rand_num <=0.25){i1--;} 
                    else if(rand_num <= 0.5){j1--;}
                    else if(rand_num <= 0.75){i1++;}
                    else{j1++;}
                    

                    if((i1 == 0) || (j1 == 0) || (i1 == n-1) || (j1 == n-1)){
                        a[i][j] += a[i1][j1] ;
                        break;
                    } 

                }

            }
            a[i][j] /= N;             
            
        }
        
    }
    
}

int main(){
    unsigned int seed = (unsigned int)time(NULL);
    double delta;
    printf("delta : ");
    scanf("%lf", &delta);
    int n = (2.0/delta)+1 ;  
    int N; 
    printf("N : ");   
    scanf("%d", &N);
    double phi[n][n];

    struct timespec start_t,end_t;
    uint64_t diff;

    clock_gettime(CLOCK_MONOTONIC,&start_t);

    boundary(n, phi, delta);

    random_dance(n, phi, delta, N, seed);
 
    saveToCSV(n, phi); 

    clock_gettime(CLOCK_MONOTONIC,&end_t);

    diff = BILLION*(end_t.tv_sec - start_t.tv_sec) + end_t.tv_nsec - start_t.tv_nsec;
    printf("\n elapsed time = %lf seconds \n",(double)diff/BILLION);
 
    return 0;
}  