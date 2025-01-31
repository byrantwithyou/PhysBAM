#include <chrono>
#include <cmath>
#include <cstdio>
#include <dmumps_c.h>
#include <mkl.h>
#include <mpi.h>
#include <omp.h>
#include <string>
#include <utility>
#include <vector>
#define USE_COMM_WORLD -987654
using namespace std;

typedef double T;

#define ICNTL(I) icntl[(I)-1]

void solve(MUMPS_INT myid,int dim,int nnz,vector<T>& entries,vector<int>& row,vector<int>& col,vector<T>& b)
{
    const int JOB_INIT=-1,JOB_END=-2,JOB_SOLVE=6;
    DMUMPS_STRUC_C id;
    id.comm_fortran=USE_COMM_WORLD;
    id.par=1;
    id.sym=2;
    id.job=JOB_INIT;
    dmumps_c(&id);

    if(myid==0)
    {
        id.n=dim;
        id.nnz=nnz;
        id.irn=&row[0];
        id.jcn=&col[0];
        id.a=&entries[0];
        id.rhs=&b[0];
    }

    // No debug outputs
    id.ICNTL(1)=-1;
    id.ICNTL(2)=-1;
    id.ICNTL(3)=-1;
    id.ICNTL(4)=0;

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
    MUMPS_INT myid;
    MPI_Init(0,0);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);

    vector<pair<chrono::steady_clock::time_point,const char*> > tm;
    tm.reserve(32);
    auto timer=[&tm](const char* name){tm.push_back({chrono::steady_clock::now(),name});};

    string dir;
    int threads;
    dir=argv[1];
    sscanf(argv[2],"%d",&threads);
    omp_set_num_threads(1);
    mkl_set_num_threads(1);

    int dim=0,nnz=0;
    vector<int> half_row,half_col;
    vector<T> half_entries,rhs,ref_sol;
    if(myid==0)
    {
        vector<int> col,offsets;
        vector<T> entries;
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
        read_bin(offsets,"/offsets.bin");
        read_bin(rhs,"/rhs.bin");
        read_bin(ref_sol,"/x.bin");
        FILE* file=fopen((dir+"/entries.bin").c_str(),"rb");
        fread(&dim,sizeof(int),1,file);
        int s;
        fread(&s,sizeof(int),1,file);
        entries.resize(s);
        fread(&entries[0],sizeof(T),s,file);
        fclose(file);

        for(int i=0;i<dim;i++)
            for(int j=offsets[i];j<offsets[i+1];j++)
            {
                int r=i+1;
                int c=col[j]+1;
                if(c>=r)
                {
                    half_row.push_back(r);
                    half_col.push_back(c);
                    half_entries.push_back(entries[j]);
                    ++nnz;
                }
            }
        printf("name %s dim %d nnz %d threads %d\n",dir.c_str(),dim,nnz,threads);
        timer("import");
    }

    solve(myid,dim,nnz,half_entries,half_row,half_col,rhs);
    timer("solve");

    if(myid==0)
        check(rhs,ref_sol);
    timer("check");

    if(myid==0)
    {
        for(int i=1;i<tm.size();i++)
            printf("%20s %5.0f ms\n",tm[i].second,
                chrono::duration_cast<chrono::duration<double> >(tm[i].first-tm[i-1].first).count()*1000);
    }
    
    return 0;
}

