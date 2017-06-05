//#####################################################################
// Copyright 2002-2007, Ron Fedkiw, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Arrays_Nd/ARRAYS_ND_VIEW.h>
#include <Core/Log/LOG.h>
#include <Grid_Tools/Fourier_Transforms/FFT.h>
namespace PhysBAM{

#ifdef USE_FFTW
#include <fftw3.h>
template<class T> struct FFT_POLICY;
template<> struct FFT_POLICY<float>{typedef fftwf_plan PLAN;};
template<> struct FFT_POLICY<double>{typedef fftw_plan PLAN;};

//fftwf_plan
//fftw_plan

inline fftw_plan plan_dft(int rank,const int *n,fftw_complex *in,fftw_complex *out,int sign,unsigned flags){return fftw_plan_dft(rank,n,in,out,sign,flags);}
inline fftw_plan plan_dft(int rank,const int *n,fftw_complex *in,double *out,int sign,unsigned flags){return fftw_plan_dft_c2r(rank,n,in,out,flags);}
inline fftw_plan plan_dft(int rank,const int *n,double *in,fftw_complex *out,int sign,unsigned flags){return fftw_plan_dft_r2c(rank,n,in,out,flags);}
inline void destroy_plan(fftw_plan plan){fftw_destroy_plan(plan);}
inline void execute_dft(const fftw_plan p,fftw_complex *in,fftw_complex *out){fftw_execute_dft(p,in,out);}
inline void execute_dft(const fftw_plan p,fftw_complex *in,double *out){fftw_execute_dft_c2r(p,in,out);}
inline void execute_dft(const fftw_plan p,double *in,fftw_complex *out){fftw_execute_dft_r2c(p,in,out);}

inline fftwf_plan plan_dft(int rank,const int *n,fftwf_complex *in,fftwf_complex *out,int sign,unsigned flags){return fftwf_plan_dft(rank,n,in,out,sign,flags);}
inline fftwf_plan plan_dft(int rank,const int *n,fftwf_complex *in,float *out,int sign,unsigned flags){return fftwf_plan_dft_c2r(rank,n,in,out,flags);}
inline fftwf_plan plan_dft(int rank,const int *n,float *in,fftwf_complex *out,int sign,unsigned flags){return fftwf_plan_dft_r2c(rank,n,in,out,flags);}
inline void destroy_plan(fftwf_plan plan){fftwf_destroy_plan(plan);}
inline void execute_dft(const fftwf_plan p,fftwf_complex *in,fftwf_complex *out){fftwf_execute_dft(p,in,out);}
inline void execute_dft(const fftwf_plan p,fftwf_complex *in,float *out){fftwf_execute_dft_c2r(p,in,out);}
inline void execute_dft(const fftwf_plan p,float *in,fftwf_complex *out){fftwf_execute_dft_r2c(p,in,out);}

inline fftw_complex* conv(const std::complex<double>* x){return (fftw_complex*)x;}
inline fftwf_complex* conv(const std::complex<float>* x){return (fftwf_complex*)x;}
inline float* conv(const float* x){return (float*)x;}
inline double* conv(const double* x){return (double*)x;}

//#####################################################################
// Destructor
//#####################################################################
template<class TV> FFT<TV>::
~FFT()
{
    for(int c=0;c<2;c++)
        for(int i=0;i<2;i++)
            if(plan[c][i])
                destroy_plan((typename FFT_POLICY<T>::PLAN)plan[c][i]);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FFT<TV>::
FFT(FFT&& p)
{
    for(int c=0;c<2;c++)
        for(int i=0;i<2;i++){
            plan[c][i]=p.plan[c][i];
            counts[c][i]=p.counts[c][i];
            p.plan[c][i]=0;}
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FFT<TV>& FFT<TV>::
operator=(FFT&& p)
{
    for(int c=0;c<2;c++)
        for(int i=0;i<2;i++){
            if(plan[c][i])
                destroy_plan((typename FFT_POLICY<T>::PLAN)plan[c][i]);
            plan[c][i]=p.plan[c][i];
            counts[c][i]=p.counts[c][i];
            p.plan[c][i]=0;}
    return *this;
}
//#####################################################################
// Function Pack_Copy
//#####################################################################
template<int n,class T> static typename std::enable_if<n==1>::type
Pack_Copy(const std::complex<T>* s,std::complex<T>* d,int* pack_sizes,int* ds,int* dd)
{
    memmove(d,s,*pack_sizes*sizeof(*s));
}
//#####################################################################
// Function Pack_Copy
//#####################################################################
template<int n,class T> static typename std::enable_if<n!=1>::type
Pack_Copy(const std::complex<T>* s,std::complex<T>* d,int* pack_sizes,int* ds,int* dd)
{
    for(int i=0;i<*pack_sizes;i++)
        Pack_Copy<n-1>(s+*ds*i,d+*dd*i,pack_sizes+1,ds+1,dd+1);
}
//#####################################################################
// Function Pack
//#####################################################################
template<class T,int d> static void
Pack(ARRAY<std::complex<T>,VECTOR<int,d> >& a)
{
    std::complex<T>* p=a.Get_Array_Pointer();
    int pack_size[d],pack_stride[d],full_stride[d];
    pack_size[d-1]=a.domain.max_corner(d-1)/2+1;
    pack_stride[d-1]=1;
    full_stride[d-1]=1;
    for(int i=d-2;i>=0;i--){
        pack_size[i]=a.domain.max_corner(i);
        pack_stride[i]=pack_stride[i+1]*pack_size[i+1];
        full_stride[i]=full_stride[i+1]*a.domain.max_corner(i+1);}
    Pack_Copy<d>(p,p,pack_size,full_stride,pack_stride);
}
//#####################################################################
// Function Pack
//#####################################################################
template<class T> static void
Pack(ARRAY<std::complex<T> >& a)
{
    std::complex<T>* p=a.Get_Array_Pointer();
    int pack_size[1]={a.m/2+1},stride[1]={1};
    Pack_Copy<1>(p,p,pack_size,stride,stride);
}
//#####################################################################
// Function Unpack_Copy
//#####################################################################
template<int n,class T> static typename std::enable_if<n==1>::type
Unpack_Copy(const std::complex<T>* s,std::complex<T>* d,int* pack_sizes,int* ds,int* dd)
{
    memmove(d,s,*pack_sizes*sizeof(*s));
}
//#####################################################################
// Function Unpack_Copy
//#####################################################################
template<int n,class T> static typename std::enable_if<n!=1>::type
Unpack_Copy(const std::complex<T>* s,std::complex<T>* d,int* pack_sizes,int* ds,int* dd)
{
    for(int i=*pack_sizes-1;i>=0;i--)
        Unpack_Copy<n-1>(s+*ds*i,d+*dd*i,pack_sizes+1,ds+1,dd+1);
}
//#####################################################################
// Function Unpack_Conj
//#####################################################################
template<int n,class T> static typename std::enable_if<n==1>::type
Unpack_Conj(const std::complex<T>* s,std::complex<T>* d,int* sizes,int* st)
{
    for(int i=1;2*i<*sizes;i++)
        d[*sizes-i]=conj(s[i]);
}
//#####################################################################
// Function Unpack_Conj
//#####################################################################
template<int n,class T> static typename std::enable_if<n!=1>::type
Unpack_Conj(const std::complex<T>* s,std::complex<T>* d,int* sizes,int* st)
{
    Unpack_Conj<n-1>(s,d,sizes+1,st+1);
    for(int i=1;i<*sizes;i++)
        Unpack_Conj<n-1>(s+*st*i,d+*st*(*sizes-i),sizes+1,st+1);
}
//#####################################################################
// Function Pack
//#####################################################################
template<class T,int d> static void
Unpack(ARRAY<std::complex<T>,VECTOR<int,d> >& a)
{
    std::complex<T>* p=a.Get_Array_Pointer();
    int pack_size[d],pack_stride[d],full_stride[d];
    pack_size[d-1]=a.domain.max_corner(d-1)/2+1;
    pack_stride[d-1]=1;
    full_stride[d-1]=1;
    for(int i=d-2;i>=0;i--){
        pack_size[i]=a.domain.max_corner(i);
        pack_stride[i]=pack_stride[i+1]*pack_size[i+1];
        full_stride[i]=full_stride[i+1]*a.domain.max_corner(i+1);}
    Unpack_Copy<d>(p,p,pack_size,pack_stride,full_stride);
    Unpack_Conj<d>(p,p,&a.domain.max_corner[0],full_stride);
}
//#####################################################################
// Function Pack
//#####################################################################
template<class T> static void
Unpack(ARRAY<std::complex<T> >& a)
{
    std::complex<T>* p=a.Get_Array_Pointer();
    int pack_size[1]={a.m/2+1},stride[1]={1},size[1]={a.m};
    Unpack_Copy<1>(p,p,pack_size,stride,stride);
    Unpack_Conj<1>(p,p,size,stride);
}
template<class T,int d> void Assert_Zero_Start(const ARRAY<T,VECTOR<int,d> >& a)
{
    PHYSBAM_ASSERT(a.domain.min_corner==(VECTOR<int,d>()));
}
template<class T> void Assert_Zero_Start(const ARRAY<T>& a){}
//#####################################################################
// Function Plan
//#####################################################################
template<class TV,class INDEX,class I,class O> static void
Plan(FFT<TV>& fft,const ARRAY<I,INDEX>& in,
    const ARRAY<O,INDEX>& out,bool inverse,bool c2c)
{
    INDEX size=in.Size();
    PHYSBAM_ASSERT(size==out.Size());
    Assert_Zero_Start(in);
    Assert_Zero_Start(out);
    if(fft.counts[c2c][inverse]==size) return;
    fft.counts[c2c][inverse]=size;
    fft.plan[c2c][inverse]=plan_dft(sizeof(size)/sizeof(int),(int*)&size,
        conv(&in(INDEX())),conv(&out(INDEX())),
        inverse?FFTW_BACKWARD:FFTW_FORWARD,
        FFTW_ESTIMATE | FFTW_UNALIGNED);
}
template<class T> static bool Empty_Helper(ARRAY<T>& a){return !a.m;}
template<class T,int d> static bool Empty_Helper(ARRAY<T,VECTOR<int,d> >& a){return a.domain.max_corner.Contains(0);}
//#####################################################################
// Function Transform
//#####################################################################
template<class TV> void FFT<TV>::
Transform(ARRAY<T,INDEX>& in,ARRAY<C,INDEX>& out)
{
    if(Empty_Helper(in)) return;
    Plan(*this,in,out,false,false);
    execute_dft((typename FFT_POLICY<T>::PLAN)plan[0][0],
        conv(&in(INDEX())),conv(&out(INDEX())));
    Unpack(out);
}
//#####################################################################
// Function Transform
//#####################################################################
template<class TV> void FFT<TV>::
Transform(ARRAY<C,INDEX>& in,ARRAY<C,INDEX>& out)
{
    if(Empty_Helper(in)) return;
    Plan(*this,in,out,false,true);
    execute_dft((typename FFT_POLICY<T>::PLAN)plan[1][0],
        conv(&in(INDEX())),conv(&out(INDEX())));
}
template<class T> void Normalize_Helper(ARRAY<T>& a){a/=(T)a.m;}
template<class T,int d> void Normalize_Helper(ARRAY<T,VECTOR<int,d> >& a){a.array/=a.domain.max_corner.Product();}
//#####################################################################
// Function Inverse_Transform
//#####################################################################
template<class TV> void FFT<TV>::
Inverse_Transform(ARRAY<C,INDEX>& in,ARRAY<T,INDEX>& out)
{
    if(Empty_Helper(in)) return;
    Plan(*this,in,out,true,false);
    Pack(in);
    execute_dft((typename FFT_POLICY<T>::PLAN)plan[0][1],
        conv(&in(INDEX())),conv(&out(INDEX())));
    Normalize_Helper(out);
}
//#####################################################################
// Function Inverse_Transform
//#####################################################################
template<class TV> void FFT<TV>::
Inverse_Transform(ARRAY<C,INDEX>& in,ARRAY<C,INDEX>& out)
{
    if(Empty_Helper(in)) return;
    Plan(*this,in,out,true,true);
    execute_dft((typename FFT_POLICY<T>::PLAN)plan[1][1],
        conv(&in(INDEX())),conv(&out(INDEX())));
    Normalize_Helper(out);
}
#else
template<class TV> FFT<TV>::
~FFT()
{PHYSBAM_FATAL_ERROR();}
template<class TV> FFT<TV>::
FFT(FFT&& p)
{PHYSBAM_FATAL_ERROR();}
template<class TV> FFT<TV>& FFT<TV>::
operator=(FFT&& p)
{PHYSBAM_FATAL_ERROR();}
template<class TV> void FFT<TV>::
Transform(ARRAY<T,INDEX>& in,ARRAY<C,INDEX>& out)
{PHYSBAM_FATAL_ERROR();}
template<class TV> void FFT<TV>::
Transform(ARRAY<C,INDEX>& in,ARRAY<C,INDEX>& out)
{PHYSBAM_FATAL_ERROR();}
template<class TV> void FFT<TV>::
Inverse_Transform(ARRAY<C,INDEX>& in,ARRAY<T,INDEX>& out)
{PHYSBAM_FATAL_ERROR();}
template<class TV> void FFT<TV>::
Inverse_Transform(ARRAY<C,INDEX>& in,ARRAY<C,INDEX>& out)
{PHYSBAM_FATAL_ERROR();}
#endif
template struct FFT<double>;
template struct FFT<VECTOR<double,1> >;
template struct FFT<VECTOR<double,2> >;
template struct FFT<VECTOR<double,3> >;
template struct FFT<float>;
template struct FFT<VECTOR<float,1> >;
template struct FFT<VECTOR<float,2> >;
template struct FFT<VECTOR<float,3> >;
}
