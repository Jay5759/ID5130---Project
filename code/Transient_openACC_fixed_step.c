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
    return -1*((x*x)+(y*y));
}

double g(double x, double y, double t){
    if(x==-1.0){
        return 100.0;
        // return 0.0;
    }
    if(x==1.0){
        return 100.0;
        // return 0.0;
    }
    if(y==-1.0){
        return 100.0;
        // return 0.0;
    }
    if(y==1.0){
        return 100.0; 
        // return 0.0;
    }
    if(t==0.0){
        return 1000.0;
    }
    return 0.0; 
}

void boundary(int n, double (*a)[n], double delta, double t){
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a[i][j] = g((i*delta)-1, (j*delta)-1, t);
        }   
    }
     
}    

double generateRandom(unsigned int *seed) {
    *seed = (*seed * MULTIPLIER + INCREMENT) % MODULUS;
    return (double)(*seed) / MODULUS; // Convert to double and scale to [0,1)
}

 
void random_dance(int n, double (*a)[n], double delta, int N, double t_delta, double time_t, double alpha, unsigned int seed){
    #pragma acc parallel loop num_gangs(10) copy(a[0:n ][0:n])
    for (int i = 0; i < (n); i++)
    {     
        for (int j = 0; j < (n); j++)
        {

            if(i==0 || i==(n-1) || j==0 || j==(n-1)){
                a[i][j] = g((i*delta)-1, (j*delta)-1, time_t); 
            }
            else{

            a[i][j] = 0.0;  
            seed += (i*68)+(j*46); 
            for (int k = 0; k < N; k++) 
            { 
                int i1 = i, j1 = j; double t1  = time_t;  
                        
                while(1){ 

                    // a[i][j] += f((i1*delta)-1, (j1*delta)-1)*delta*delta/4;

                    // unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                    // std::mt19937 generator(seed);
                    // std::uniform_real_distribution<double> distribution(0.0, 1.0);
                    // double rand_num = distribution(generator); 
                    double rand_num = generateRandom(&seed);
                    if(rand_num <= (alpha*t_delta)/(delta*delta)){i1--;}
                    else if(rand_num <= (2*alpha*t_delta)/(delta*delta)){j1--;}
                    else if(rand_num <= (3*alpha*t_delta)/(delta*delta)){i1++;}
                    else if(rand_num <= (4*alpha*t_delta)/(delta*delta)){j1++;}
                    t1 -= t_delta;

                    if((i1 == 0) || (j1 == 0) || (i1 == n-1) || (j1 == n-1) || t1<=0.0){
                        if(t1<0.0)t1 = 0.0;
                        a[i][j] += g((i1*delta)-1, (j1*delta)-1, t1) ; 
                        break;
                    }  

                }  
 
            }
            a[i][j] /= N;             
            
            }
              
        }
        
    }
    
}

int main(){
    unsigned int seed = (unsigned int)time(NULL);
    double delta; 
    printf("delta : ");
    scanf("%lf", &delta);

    double t_delta; 
    printf("t_delta : ");
    scanf("%lf", &t_delta);
    double time_t;
    printf("time_t : ");
    scanf("%lf", &time_t);
    double alpha;
    printf("alpha : ");
    scanf("%lf", &alpha);
    int N; 
    printf("N : ");
    scanf("%d", &N); 
    int n = (2.0/delta)+1 ;  
    double phi_0[n][n];
    double phi_t[n][n]; 

    struct timespec start_t,end_t;
    uint64_t diff;

    clock_gettime(CLOCK_MONOTONIC,&start_t);

    // vector<vector<double>> phi_0(n,vector<double>(n,0));
    // vector<vector<double>> phi_t(n,vector<double>(n,0));


    boundary(n, phi_0, delta, 0.0);
    random_dance(n, phi_t, delta, N, t_delta, time_t, alpha, seed);
 
    saveToCSV(n, phi_t); 

    clock_gettime(CLOCK_MONOTONIC,&end_t);

    diff = BILLION*(end_t.tv_sec - start_t.tv_sec) + end_t.tv_nsec - start_t.tv_nsec;
    printf("\n elapsed time = %lf seconds \n",(double)diff/BILLION);

    return 0; 
} 