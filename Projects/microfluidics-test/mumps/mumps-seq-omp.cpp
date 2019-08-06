#include <cstdio>
#include <cmath>
#include <chrono>
#include <utility>
#include <vector>
#include <string>
#include <mpi.h>
#include <dmumps_c.h>
#include <mkl.h>
#include <omp.h>
#define USE_COMM_WORLD -987654
using namespace std;

typedef double T;

#define ICNTL(I) icntl[(I)-1]

void solve(int dim,int nnz,vector<T>& entries,vector<int>& row,vector<int>& col,vector<T>& b)
{
    const int JOB_INIT=-1,JOB_END=-2,JOB_SOLVE=6;
    MUMPS_INT myid;
    MPI_Init(0,0);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);

    DMUMPS_STRUC_C id;
    id.comm_fortran=USE_COMM_WORLD;
    id.par=1;
    id.sym=2;
    id.job=JOB_INIT;
    dmumps_c(&id);

    id.n=dim;
    id.nnz=nnz;
    id.irn=&row[0];
    id.jcn=&col[0];
    id.a=&entries[0];

    // No debug outputs
    id.ICNTL(1)=-1;
    id.ICNTL(2)=-1;
    id.ICNTL(3)=-1;
    id.ICNTL(4)=0;

    id.rhs=&b[0];
    id.job=JOB_SOLVE;
    dmumps_c(&id);
    if(id.infog[0]<0)
        printf("MUMPS ERROR %d %d\n",id.infog[0],id.infog[1]);
    id.job=JOB_END;
    dmumps_c(&id);
    MPI_Finalize();
}

void check(const vector<T>& x,const vector<T>& ref)
{
    T err_linf=0;
    for(int i=0;i<x.size();i++)
        err_linf=max(err_linf,fabs(x[i]-ref[i]));
    printf("error %e\n",err_linf);
}

// ./prog <output_dir> <threads>
int main(int argc,char* argv[])
{
    vector<pair<chrono::steady_clock::time_point,const char*> > tm;
    tm.reserve(32);
    auto timer=[&tm](const char* name){tm.push_back({chrono::steady_clock::now(),name});};

    string dir;
    int threads;
    dir=argv[1];
    sscanf(argv[2],"%d",&threads);
    omp_set_num_threads(threads);
    mkl_set_num_threads(threads);

    vector<int> col,row;
    vector<T> rhs,entries,ref_sol;
    int dim,nnz;
    auto read_bin=[&dir](auto& arr,const char* name)
    {
        FILE* file=fopen((dir+name).c_str(),"rb");
        int s;
        fread(&s,sizeof(int),1,file);
        arr.resize(s);
        fread(&arr[0],sizeof(arr[0]),s,file);
        fclose(file);
    };
    timer("init");

    read_bin(col,"/col.bin");
    read_bin(row,"/row.bin");
    for(auto& e:col) e+=1;
    for(auto& e:row) e+=1;
    read_bin(rhs,"/rhs.bin");
    read_bin(ref_sol,"/x.bin");
    FILE* file=fopen((dir+"/entries.bin").c_str(),"rb");
    fread(&dim,sizeof(int),1,file);
    fread(&nnz,sizeof(int),1,file);
    entries.resize(nnz);
    fread(&entries[0],sizeof(T),nnz,file);
    fclose(file);
    printf("name %s dim %d nnz %d threads %d\n",dir.c_str(),dim,nnz,threads);
    timer("import");

    solve(dim,nnz,entries,row,col,rhs);
    timer("solve");

    check(rhs,ref_sol);
    timer("check");

    for(int i=1;i<tm.size();i++)
        printf("%20s %5.0f ms\n",tm[i].second,
            chrono::duration_cast<chrono::duration<double> >(tm[i].first-tm[i-1].first).count()*1000);
    
    return 0;
}

