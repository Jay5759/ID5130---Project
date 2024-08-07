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
double min(double a, double b) {
    return (a < b) ? a : b;
}

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

double closest_boundary_dist(double x, double y, double delta, int n){

    double dist = 0.0;
    dist = min(1.0+x, min(1.0-x,min(1.0+y,1.0-y)));
    return dist; 
}

double generateRandom(unsigned int *seed) {
    *seed = (*seed * MULTIPLIER + INCREMENT) % MODULUS;
    return (double)(*seed) / MODULUS; // Convert to double and scale to [0,1)
}

double H(double h){
    if(h<0.1){return ((0.01312)+(3.3082*h)+(-91.011*h*h)+(1348.1*h*h*h)+(-9524.2*h*h*h*h)+(25594*h*h*h*h*h));}
    if(h<0.3){return ((0.052654)+(0.36498*h)+(-0.45109*h*h)+(0.66164*h*h*h));}
    if(h<0.6){return ((0.051155)+(0.35391*h)+(-0.33104*h*h)+(0.44125*h*h*h));}
    return -0.17292*log(0.62423*(1-h));
}

 
void random_dance(int n, double (*a)[n], double delta, int N,double r_min, double t_delta, double time_t, double alpha, unsigned int seed){
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
                double xi = (i*delta)-1, yi = (j*delta)-1;
                double ti = time_t; 
                         
                while(1){ 

                    double min_dist = closest_boundary_dist(xi, yi, delta, n);

                    if(min_dist<r_min || ti<=0.0){
                        if(ti<=0.0){
                            a[i][j] += g(xi, yi, 0.0);
                        }
                        else if(min_dist == 1.0-xi){
                            a[i][j] += g(1.0, yi, ti);
                        }
                        else if(min_dist == 1.0+xi){
                            a[i][j] += g(-1.0, yi, ti);
                        }
                        else if(min_dist == 1.0-yi){
                            a[i][j] += g(xi, 1.0, ti);
                        }
                        else if(min_dist == 1.0+yi){
                            a[i][j] += g(xi, -1.0, ti);
                        }
                        break;
                    }
                    double rand_num = generateRandom(&seed);
                    xi += (min_dist * cos(rand_num*pi*2));
                    yi += (min_dist * sin(rand_num*pi*2));

                    double rand_num2 = generateRandom(&seed);
                    ti -= (min_dist*min_dist/alpha)*H(rand_num2);
 

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
    double r_min = delta/4;

    struct timespec start_t,end_t;
    uint64_t diff;

    clock_gettime(CLOCK_MONOTONIC,&start_t);

    // vector<vector<double>> phi_0(n,vector<double>(n,0));
    // vector<vector<double>> phi_t(n,vector<double>(n,0));


    // boundary(n, phi_0, delta, 0.0);
    random_dance(n, phi_t, delta, N, r_min, t_delta, time_t, alpha, seed);
 
    saveToCSV(n, phi_t); 

    clock_gettime(CLOCK_MONOTONIC,&end_t);

    diff = BILLION*(end_t.tv_sec - start_t.tv_sec) + end_t.tv_nsec - start_t.tv_nsec;
    printf("\n elapsed time = %lf seconds \n",(double)diff/BILLION);

    return 0; 
}  