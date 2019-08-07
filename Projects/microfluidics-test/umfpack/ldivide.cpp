#include <cstdio>
#include <cmath>
#include <chrono>
#include <utility>
#include <vector>
#include <string>
#include <suitesparse/umfpack.h>
#include <omp.h>
#include <mkl.h>
using namespace std;

typedef double T;

void check(const vector<T>& x,const vector<T>& ref)
{
    T err_linf=0;
    for(int i=0;i<x.size();i++)
        err_linf=max(err_linf,fabs(x[i]-ref[i]));
    printf("error %e\n",err_linf);
}

// ./prog <output_dir>
int main(int argc,char* argv[])
{
    vector<pair<chrono::steady_clock::time_point,const char*> > tm;
    tm.reserve(32);
    auto timer=[&tm](const char* name){tm.push_back({chrono::steady_clock::now(),name});};

    string dir;
    dir=argv[1];
    int threads;
    sscanf(argv[2],"%d",&threads);
    omp_set_num_threads(threads);
    mkl_set_num_threads(threads);
    timer("init");

    vector<int> col,offsets;
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

    read_bin(col,"/col.bin");
    read_bin(offsets,"/offsets.bin");
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

    void* symbolic;
    umfpack_di_symbolic(dim,dim,&offsets[0],&col[0],&entries[0],&symbolic,0,0);
    void* numeric;
    umfpack_di_numeric(&offsets[0],&col[0],&entries[0],symbolic,&numeric,0,0);
    vector<T> x(dim);
    umfpack_di_solve(UMFPACK_A,&offsets[0],&col[0],&entries[0],&x[0],&rhs[0],numeric,0,0);
    timer("solve");

    check(x,ref_sol);
    timer("check");

    umfpack_di_free_symbolic(&symbolic);
    umfpack_di_free_numeric(&numeric);

    for(int i=1;i<tm.size();i++)
        printf("%20s %5.0f ms\n",tm[i].second,
            chrono::duration_cast<chrono::duration<double> >(tm[i].first-tm[i-1].first).count()*1000);
    
    return 0;
}

