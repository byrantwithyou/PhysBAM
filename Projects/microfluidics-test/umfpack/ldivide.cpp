#include <cassert>
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

vector<T> residue(int dim,const vector<int>& offsets,const vector<int>& col,const vector<T>& entries,
    const vector<T>& x,const vector<T>& b)
{
    vector<T> result(dim);
    int index=offsets[0];
    for(int i=0;i<dim;i++){
        int end=offsets[i+1];
        T sum=(T)0;
        for(;index<end;index++)
            sum+=entries[index]*x[col[index]];
        result[i]=sum-b[i];
    }
    return result;
}

T l2(const vector<T>& x)
{
    T s=0;
    for(const auto& a:x)
    {
        s+=a*a;
    }
    s/=x.size();
    return sqrt(s);
}

T linf(const vector<T>& x)
{
    T s=abs(x[0]);
    for(const auto& a:x)
    {
        s=max(s,abs(a));
    }
    return s;
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

    vector<long int> col_l,offsets_l;
    col_l.resize(col.size());
    offsets_l.resize(offsets.size());
    for(int i=0;i<col.size();i++) col_l[i]=(long int)col[i];
    for(int i=0;i<offsets.size();i++) offsets_l[i]=(long int)offsets[i];
    timer("import");

    double ctrl[UMFPACK_CONTROL];
    double info[UMFPACK_INFO];
    umfpack_dl_defaults(ctrl);
    //ctrl[UMFPACK_STRATEGY]=UMFPACK_STRATEGY_SYMMETRIC;
    //ctrl[UMFPACK_PIVOT_TOLERANCE]=0.1;
    void* symbolic;
    long int ret;
    ret=umfpack_dl_symbolic((long)dim,(long)dim,&offsets_l[0],&col_l[0],&entries[0],&symbolic,ctrl,info);
    //printf("symbolic %ld\n",ret);
    assert(ret==UMFPACK_OK);
    void* numeric;
    ret=umfpack_dl_numeric(&offsets_l[0],&col_l[0],&entries[0],symbolic,&numeric,ctrl,info);
    //printf("numeric %ld\n",ret);
    assert(ret==UMFPACK_OK);
    vector<T> x(dim);
    ret=umfpack_dl_solve(UMFPACK_A,&offsets_l[0],&col_l[0],&entries[0],&x[0],&rhs[0],numeric,ctrl,info);
    //printf("ret %ld\n",ret);
    assert(ret==UMFPACK_OK);
    timer("solve");

    printf("INFO NZ %.16f\n",info[UMFPACK_NZ]);
    printf("INFO RCOND %.16f\n",info[UMFPACK_RCOND]);
    printf("INFO RSMIN %.16f\n",info[UMFPACK_RSMIN]);
    printf("INFO RSMAX %.16f\n",info[UMFPACK_RSMAX]);
    printf("INFO NOFF_DIAG %.16f\n",info[UMFPACK_NOFF_DIAG]);
    printf("INFO NZDROPPED %.16f\n",info[UMFPACK_NZDROPPED]);
    timer("info");

    check(x,ref_sol);
    auto res=residue(dim,offsets,col,entries,x,rhs);
    printf("residue linf %e l2 %e\n",linf(res),l2(res));
    timer("check");

    umfpack_dl_free_symbolic(&symbolic);
    umfpack_dl_free_numeric(&numeric);

    for(int i=1;i<tm.size();i++)
        printf("%20s %5.0f ms\n",tm[i].second,
            chrono::duration_cast<chrono::duration<double> >(tm[i].first-tm[i-1].first).count()*1000);
    
    return 0;
}

