#include<bits/stdc++.h>
#include<mpi.h>
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
    ofstream file("data.csv");
    

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
    return 1*((x*x)+(y*y));
    // return 0.0; 
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
 
void random_dance(vector<vector<double>> &a, double delta, int N, int nprocs, int myid, int n, double t_delta, double time, double alpha){
    
    int i_myid = 0; 
    for (int i = 0; i < n; i++)
    {
        if((i%nprocs) == myid){
        for (int j = 0; j < (a[0].size()); j++)
        {  
            if(i==0 || i==n-1 || j==0 || j==(a[0].size()-1)){
                a[i_myid][j] = g((i*delta)-1, (j*delta)-1, time);
            }
            else{
            for (int k = 0; k < N; k++) 
            {
                int i1 = i, j1 = j; double t1 = time;  
                while(1){

                    // a[i_myid][j] += f((i1*delta)-1, (j1*delta)-1)*delta*delta/4;

                    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                    std::mt19937 generator(seed);
                    std::uniform_real_distribution<double> distribution(0.0, 1.0);
                    double rand_num = distribution(generator); 
                    if(rand_num <= (alpha*t_delta)/(delta*delta)){i1--;}
                    else if(rand_num <= (2*alpha*t_delta)/(delta*delta)){j1--;}
                    else if(rand_num <= (3*alpha*t_delta)/(delta*delta)){i1++;}
                    else if(rand_num <= (4*alpha*t_delta)/(delta*delta)){j1++;}
                    t1 -= t_delta;
                    

                    if((i1 == 0) || (j1 == 0) || (i1 == n-1) || (j1 == a[0].size()-1) || (t1<=0.0)){
                        if(t1<0.0)t1=0.0;
                        a[i_myid][j] += g((i1*delta)-1, (j1*delta)-1, t1); 
                        break;
                    } 

                }
 
            }
            a[i_myid][j] /= N;  
            }             
            
        } 
        i_myid++;
        }
  

        
    }
    
}

int32_t main(int argc, char* argv[]){

    int myid, nprocs, proc;

    struct timespec start_t,end_t;
    uint64_t diff;

    clock_gettime(CLOCK_MONOTONIC,&start_t);
    
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    double delta = 0.0, t_delta = 0.0, time = 0.0, alpha = 0.0;
    

    int N = 1; 
    if(myid == 0){
    cout<<"delta : "; cin>>delta;
    cout<<"t_delta : "; cin>>t_delta;
    cout<<"time : "; cin>>time;
    cout<<"alpha : "; cin>>alpha;
    cout<<"N : "; cin>>N; 
    }
    MPI_Bcast(&delta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&t_delta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);


    int n = (2.0/delta)+1 ;  
    
    
    vector<vector<double>> phi((n/nprocs)+nprocs,vector<double>(n,0));


    random_dance(phi, delta, N, nprocs, myid, n, t_delta, time, alpha);
 
      
    if(myid!=0){
        int limit = (n/nprocs);
        if(myid<(n%nprocs)){limit++;}
        for (int i = 0; i < limit; i++)
        {
            MPI_Send(&phi[i][0], n, MPI_DOUBLE, 0, (i*nprocs)+myid, MPI_COMM_WORLD);
        }
         
    }
    else{  
        vector<vector<double>> final_phi(n,vector<double>(n,0));
        int i_myid_0 = 0;
        for(int i=0; i<n; i++){
 
            if((i%nprocs) == myid){
                for (int j = 0; j < n; j++)
                { 
                    final_phi[i][j] = phi[i_myid_0][j];
                }
                i_myid_0++; 
                
            }  
            else{
                MPI_Recv(&final_phi[i][0], n, MPI_DOUBLE, (i%nprocs), i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }  
    saveToCSV(final_phi);  
    }

    if(myid==0){
        clock_gettime(CLOCK_MONOTONIC,&end_t);

        diff = BILLION*(end_t.tv_sec - start_t.tv_sec) + end_t.tv_nsec - start_t.tv_nsec;
        printf("\n elapsed time = %lf seconds \n",(double)diff/BILLION);
    }

    MPI_Finalize(); 
    return 0;
}   