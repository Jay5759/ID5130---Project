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

void random_dance(vector<vector<double>> &a, double delta, int N, double t_delta, double time, double alpha){
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
                int i1 = i, j1 = j; double t1 = time; 
                while(1){

                    // a[i][j] += f((i1*delta)-1, (j1*delta)-1)*delta*delta/4;

                    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                    std::mt19937 generator(seed);
                    std::uniform_real_distribution<double> distribution(0.0, 1.0);
                    double rand_num = distribution(generator); 
                    if(rand_num <= (alpha*t_delta)/(delta*delta)){i1--;}
                    else if(rand_num <= (2*alpha*t_delta)/(delta*delta)){j1--;}
                    else if(rand_num <= (3*alpha*t_delta)/(delta*delta)){i1++;}
                    else if(rand_num <= (4*alpha*t_delta)/(delta*delta)){j1++;}
                    t1 -= t_delta;

                    if((i1 == 0) || (j1 == 0) || (i1 == a.size()-1) || (j1 == a[0].size()-1) || t1<=0.0){
                        if(t1<0.0)t1=0.0;
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

int32_t main(){

    double delta; cout<<"delta : "; cin>>delta;
    double t_delta; cout<<"t_delta : "; cin>>t_delta;
    double time; cout<<"time : "; cin>>time;
    double alpha; cout<<"alpha : "; cin>>alpha;
    int N; cout<<"N : "; cin>>N;  
    int n = (2.0/delta)+1 ;  
    vector<vector<double>> phi_0(n,vector<double>(n,0));
    vector<vector<double>> phi_t(n,vector<double>(n,0));

    struct timespec start_t,end_t;
    uint64_t diff;

    clock_gettime(CLOCK_MONOTONIC,&start_t);

    boundary(phi_0, delta, 0.0);

    random_dance(phi_t, delta, N, t_delta, time, alpha);
 
    saveToCSV(phi_t); 

    clock_gettime(CLOCK_MONOTONIC,&end_t);

    diff = BILLION*(end_t.tv_sec - start_t.tv_sec) + end_t.tv_nsec - start_t.tv_nsec;
    printf("\n elapsed time = %lf seconds \n",(double)diff/BILLION);

    return 0; 
}  