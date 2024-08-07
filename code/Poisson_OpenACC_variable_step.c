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


double closest_boundary_dist(double x, double y, double delta, int n){

    double dist = 0.0;
    dist = min(1.0+x, min(1.0-x,min(1.0+y,1.0-y)));
    return dist; 
} 

double generateRandom(unsigned int *seed) {
    *seed = (*seed * MULTIPLIER + INCREMENT) % MODULUS;
    return (double)(*seed) / MODULUS; // Convert to double and scale to [0,1)
}

void random_dance(int n, double (*a)[n], double delta, int N, double r_min, unsigned int seed){
    #pragma acc parallel loop num_gangs(10) copy(a[0:n ][0:n])
    for (int i = 1; i < (n-1); i++)
    {
        for (int j = 1; j < (n-1); j++)
        {
            a[i][j] = 0.0; 
            seed += (i*68)+(j*46);
            for (int k = 0; k < N; k++) 
            {
                double xi = (i*delta)-1, yi = (j*delta)-1;
                while(1){

                    double min_dist = closest_boundary_dist(xi, yi, delta, n);
                    
                    a[i][j] += f(xi,yi)*min_dist*min_dist/4;
                    if(min_dist<r_min){

                        if(min_dist == 1.0-xi){
                            a[i][j] += g(1.0, yi);
                        }
                        else if(min_dist == 1.0+xi){
                            a[i][j] += g(-1.0, yi);
                        }
                        else if(min_dist == 1.0-yi){
                            a[i][j] += g(xi, 1.0);
                        }
                        else if(min_dist == 1.0+yi){
                            a[i][j] += g(xi, -1.0);
                        }
                        break;

                    }
                    // unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                    // std::mt19937 generator(seed);
                    // std::uniform_real_distribution<double> distribution(0.0, 1.0);
                    // double rand_num = distribution(generator);
                    double rand_num = generateRandom(&seed); 

                    xi += (min_dist * cos(rand_num*pi*2));
                    yi += (min_dist * sin(rand_num*pi*2));
                }

            }
            a[i][j] /= N;             
            
        }
        
    }
    
}

int main(){
    unsigned int seed = (unsigned int)time(NULL);
    double delta; 
    // cin>>delta;
    printf("delta : ");
    scanf("%lf", &delta);
    int n = (2.0/delta)+1 ; 
    double r_min = delta/4;
    int N; 
    printf("N : "); 
    scanf("%d", &N);
    // vector<vector<double>> phi(n,vector<double>(n,0));
    double phi[n][n];

    struct timespec start_t,end_t;
    uint64_t diff;

    clock_gettime(CLOCK_MONOTONIC,&start_t);

    boundary(n,phi,delta);

    random_dance(n, phi, delta, N, r_min, seed);
 
    saveToCSV(n, phi);

    clock_gettime(CLOCK_MONOTONIC,&end_t);

    diff = BILLION*(end_t.tv_sec - start_t.tv_sec) + end_t.tv_nsec - start_t.tv_nsec;
    printf("\n elapsed time = %lf seconds \n",(double)diff/BILLION);

    return 0;
} 