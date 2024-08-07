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

void boundary(vector<vector<double>> &a, double delta){
    
    for (int i = 0; i < a[0].size(); i++)
    { 
        a[0][i] = g(-1, (i*delta)-1);
        a[a.size()-1][i] = g(1, (i*delta)-1);
    }
    for (int i = 0; i < a.size(); i++)
    { 
        a[i][0] = g((i*delta)-1, -1);
        a[i][a[0].size()-1] = g((i*delta)-1, 1);
    }
      
}    

void random_dance(vector<vector<double>> &a, double delta, int N){
    #pragma omp parallel for num_threads(8) collapse(2)  
    for (int i = 1; i < (a.size()-1); i++)
    {
        for (int j = 1; j < (a[0].size()-1); j++)
        {
            for (int k = 0; k < N; k++) 
            { 
                int i1 = i, j1 = j;  
                while(1){

                    a[i][j] += f((i1*delta)-1, (j1*delta)-1)*delta*delta/4;

                    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                    std::mt19937 generator(seed);
                    std::uniform_real_distribution<double> distribution(0.0, 1.0);
                    double rand_num = distribution(generator); 
                    if(rand_num <=0.25){i1--;}
                    else if(rand_num <= 0.5){j1--;}
                    else if(rand_num <= 0.75){i1++;}
                    else{j1++;}
                    

                    if((i1 == 0) || (j1 == 0) || (i1 == a.size()-1) || (j1 == a[0].size()-1)){
                        a[i][j] += a[i1][j1] ;
                        break; 
                    } 

                }

            }
            a[i][j] /= N;             
            
        }
        
    }
    
}

int32_t main(){

    double delta;
    cin>>delta;
    int n = (2.0/delta)+2;
    int N; cin>>N; 
    vector<vector<double>> phi(n,vector<double>(n,0));

    struct timespec start_t,end_t;
    uint64_t diff;

    clock_gettime(CLOCK_MONOTONIC,&start_t);

    boundary(phi,delta);

    random_dance(phi, delta, N);
 
    saveToCSV(phi);

    clock_gettime(CLOCK_MONOTONIC,&end_t);

    diff = BILLION*(end_t.tv_sec - start_t.tv_sec) + end_t.tv_nsec - start_t.tv_nsec;
    printf("\n elapsed time = %lf seconds \n",(double)diff/BILLION);

    return 0;
}   