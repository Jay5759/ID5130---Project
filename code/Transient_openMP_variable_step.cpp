#include<bits/stdc++.h>
#include<omp.h>
#define BILLION 1000000000L
const double pi = 3.141592;

 
using namespace std;

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

void saveToCSV(vector<vector<double>>& data) {
    ofstream file("data2.csv");
    
 
    for (const auto& row : data) { 
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i != row.size() - 1) {
                file << ",";
            }
        }
        file << endl;
    }

    file.close();
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

void boundary(vector<vector<double>> &a, double delta, double t){
    for (int i = 0; i < a.size(); i++)
    {
        for (int j = 0; j < a[0].size(); j++)
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

double H(double h){
    if(h<0.1){return ((0.01312)+(3.3082*h)+(-91.011*h*h)+(1348.1*h*h*h)+(-9524.2*h*h*h*h)+(25594*h*h*h*h*h));}
    if(h<0.3){return ((0.052654)+(0.36498*h)+(-0.45109*h*h)+(0.66164*h*h*h));}
    if(h<0.6){return ((0.051155)+(0.35391*h)+(-0.33104*h*h)+(0.44125*h*h*h));}
    return -0.17292*log(0.62423*(1-h));
}

void random_dance(vector<vector<double>> &a, double delta, int N,double r_min, int n, double t_delta, double time, double alpha){
    #pragma omp parallel for num_threads(8) collapse(2)

    for (int i = 0; i < (a.size()); i++)
    {   
        for (int j = 0; j < (a[0].size()); j++)
        {
            if(i==0 || i==(a.size()-1) || j==0 || j==(a[0].size()-1)){
                a[i][j] = g((i*delta)-1, (j*delta)-1, time);
            }
            else{
            for (int k = 0; k < N; k++) 
            {
                double xi = (i*delta)-1, yi = (j*delta)-1;
                double ti = time; 
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
                    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                    std::mt19937 generator(seed);
                    std::uniform_real_distribution<double> distribution(0.0, 1.0);
                    double rand_num = distribution(generator); 
                    xi += (min_dist * cos(rand_num*pi*2));
                    yi += (min_dist * sin(rand_num*pi*2));

                    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
                    std::mt19937 generator1(seed1);
                    std::uniform_real_distribution<double> distribution1(0.0, 1.0);
                    double rand_num2 = distribution1(generator1); 
                    ti -= (min_dist*min_dist/alpha)*H(rand_num2);

                } 

            }
            a[i][j] /= N;             
            }
             
        }
        
    }
    
}

int32_t main(){

    double delta; cout<<"delta : "; cin>>delta;
    double t_delta; cout<<"t_delta : "; cin>>t_delta;
    double time; cout<<"time : "; cin>>time;
    double alpha; cout<<"alpha : "; cin>>alpha;
    int N; cout<<"N : "; cin>>N;  
    int n = (2.0/delta)+1 ;  
    vector<vector<double>> phi_0(n,vector<double>(n,0));
    vector<vector<double>> phi_t(n,vector<double>(n,0));
    double r_min = delta/4;

    struct timespec start_t,end_t;
    uint64_t diff;

    clock_gettime(CLOCK_MONOTONIC,&start_t);

    // boundary(phi_0, delta, 0.0);

    random_dance(phi_t, delta, N, r_min, n, t_delta, time, alpha);
 
    saveToCSV(phi_t); 

    clock_gettime(CLOCK_MONOTONIC,&end_t);

    diff = BILLION*(end_t.tv_sec - start_t.tv_sec) + end_t.tv_nsec - start_t.tv_nsec;
    printf("\n elapsed time = %lf seconds \n",(double)diff/BILLION);

    return 0;  
} 