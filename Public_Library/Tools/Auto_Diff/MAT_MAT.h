//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MAT_MAT
//##################################################################### 
#ifndef __MAT_MAT__
#define __MAT_MAT__

#include <Tools/Auto_Diff/COL_MAT.h>
#include <Tools/Auto_Diff/PRIMITIVE_MATRICES.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{
template<class TV,class H> struct MAT_MAT;

template<int nn,int m,int s> struct MM_FLAGS
{
    enum WA{n=nn,mask=m,sp_mask=s};
    STATIC_ASSERT(mask<(1<<(n*(n+1)/2)) && sp_mask<(1<<(n*(n+1)/2)));

    MM_FLAGS<n,m,s> operator+() const;
    MM_FLAGS<n,m|s,s> operator-() const;

    template<int m2,int s2> typename H_ADD<MM_FLAGS,MM_FLAGS<n,m2,s2> >::MM operator+(MM_FLAGS<n,m2,s2> x) const;
    template<int m2,int s2> typename H_ADD<MM_FLAGS,MM_FLAGS<n,s2|m2,s2> >::MM operator-(MM_FLAGS<n,m2,s2> x) const;
    MM_FLAGS<n,m|s,s> operator*(double a) const;
    MM_FLAGS<n,m|s,s> operator/(double a) const;
    MM_FLAGS<n-1,(mask>>n),(sp_mask>>n)> static Base();
    CM_FLAGS<n-1,((mask&((1<<n)-1))>>1),((sp_mask&(1<<n)-1)>>1)> static Col();
    template<int m2,int s2,int m3,int s3> MM_FLAGS<n+1,m2|(m3<<1)|(m<<(n+1)),s2|(s3<<1)|(s<<(n+1))> Append(M_FLAGS<m2,s2>,CM_FLAGS<n,m3,s3>);
    M_FLAGS<mask&1,sp_mask&1> static First();
};
typedef MM_FLAGS<0,0,0> MM_EMPTY;
template<int n,int m,int m2,int s,int s2> MM_FLAGS<n,m|m2,(~m|s)&(~m2|s2)&(m|m2|s|s2)> Choose(MM_FLAGS<n,m,s>,MM_FLAGS<n,m2,s2>);

template<int n,int m1,int m2> struct H_SYM_OUTER;
template<int n> struct H_SYM_OUTER<n,0,0> {enum WA {mask=0,sp_mask=0};};
template<int n,int m1,int m2> struct H_SYM_OUTER
{enum WA {mask=((m1&1)*m2)|((m2&1)*m1)|(H_SYM_OUTER<n-1,(m1>>1),(m2>>1)>::mask<<n),sp_mask=0};};

template<int n,int m1,int m2> MM_FLAGS<n,H_SYM_OUTER<n,m1,m2>::mask,H_SYM_OUTER<n,m1,m2>::sp_mask> Sym_Outer(CV_FLAGS<n,m1>,CV_FLAGS<n,m2>);

template<int n,int m> MM_FLAGS<n,H_SYM_OUTER<n,m,m>::mask,H_SYM_OUTER<n,m,m>::sp_mask> Outer(CV_FLAGS<n,m>);

template<int m,int s,int m2,int s2> M_FLAGS<(m|s)&(m2|s2),s&s2> Stt(M_FLAGS<m,s>,M_FLAGS<m2,s2>);

template<class H,class H2> struct STT_HELPER;
template<> struct STT_HELPER<CM_EMPTY,CM_EMPTY> {typedef MM_EMPTY TYPE;};

template<int n,int m,int s,int m2,int s2>
struct STT_HELPER<CM_FLAGS<n,m,s>,CM_FLAGS<n,m2,s2> >
{
    typedef typename STT_HELPER<decltype(CM_FLAGS<n,m,s>::Base()),decltype(CM_FLAGS<n,m2,s2>::Base())>::TYPE B;
    typedef decltype(Stt(CM_FLAGS<n,m,s>::First(),CM_FLAGS<n,m2,s2>::First())) F;
    typedef decltype(CM_FLAGS<n,m,s>::Base()*CM_FLAGS<n,m2,s2>::First()+CM_FLAGS<n,m2,s2>::Base()*CM_FLAGS<n,m,s>::First()) C;
    typedef decltype(B().Append(F(),C())) TYPE;
};

template<int n,int m,int s,int m2,int s2> typename STT_HELPER<CM_FLAGS<n,m,s>,CM_FLAGS<n,m2,s2> >::TYPE Stt(CM_FLAGS<n,m,s>,CM_FLAGS<n,m2,s2>);

template<int m,int s,int m2,int s2> M_FLAGS<m,s> Tst(M_FLAGS<m,s>);

template<class H> struct TST_HELPER;
template<> struct TST_HELPER<CM_EMPTY> {typedef MM_EMPTY TYPE;};

template<int n,int m,int s>
struct TST_HELPER<CM_FLAGS<n,m,s> >
{
    typedef typename TST_HELPER<decltype(CM_FLAGS<n,m,s>::Base())>::TYPE B;
    typedef decltype(Tst(CM_FLAGS<n,m,s>::First())) F;
    typedef decltype(CM_FLAGS<n,m,s>::Base()*CM_FLAGS<n,m,s>::First()) C;
    typedef decltype(B().Append(F(),C())) TYPE;
};

template<int n,int m,int s> typename TST_HELPER<CM_FLAGS<n,m,s> >::TYPE Tst(CM_FLAGS<n,m,s>);

#define MK_MM(TV,FLAGS) MAT_MAT<TV,decltype(FLAGS)>

template<class TV>
struct MAT_MAT<TV,MM_EMPTY>
{
    typedef MM_EMPTY FLAGS;
    typedef typename TV::SCALAR T;

    MATRIX<T,TV::m> operator()(int i,int j) const {return MATRIX<T,TV::m>();}

    SYMMETRIC_MATRIX<T,TV::m> operator()(int i) const {return SYMMETRIC_MATRIX<T,TV::m>();}

    void Neg(const MAT_MAT& u) {}

    void Scale(const MAT_MAT& u,T a) {}

    void Add(const MAT_MAT& u,const MAT_MAT& v) {}

    void Sub(const MAT_MAT& u,const MAT_MAT& v) {}

    void Fill_From(const MAT_MAT& u) {}
};

template<class TV,class H>
struct MAT_MAT:public MK_MM(TV,H::Base())
{
    typedef H FLAGS;
    typedef typename TV::SCALAR T;

    typedef MK_CM(TV,H::Col()) T_COL;
    typedef MK_SYM_MAT(TV,H::First()) T_DATA;
    T_DATA x;
    T_COL c;

    typedef MK_MM(TV,H::Base()) BASE;
    const BASE& Base() const {return static_cast<const BASE&>(*this);}
    BASE& Base() {return static_cast<BASE&>(*this);}

    SYMMETRIC_MATRIX<T,TV::m> operator()(int i) const {if(i==0) return Cast_Helper(x);return Base()(i-1);}

    // i>=j
    MATRIX<T,TV::m> operator()(int i,int j) const
    {
        if(j==0){
            if(i==0) return Cast_Helper(x);
            return c(i-1);}
        return Base()(i-1,j-1);
    }

    template<class H2>
    void Neg(const MAT_MAT<TV,H2>& u)
    {x=-u.x;c.Neg(u.c);Base().Neg(u.Base());}

    template<class H2>
    void Scale(const MAT_MAT<TV,H2>& u,T a)
    {x=u.x*a;c.Scale(u.c,a);Base().Scale(u.Base(),a);}

    template<class H2,class H3>
    void Add(const MAT_MAT<TV,H2>& u,const MAT_MAT<TV,H3>& v)
    {x=u.x+v.x;c.Add(u.c,v.c);Base().Add(u.Base(),v.Base());}

    template<class H2,class H3>
    void Sub(const MAT_MAT<TV,H2>& u,const MAT_MAT<TV,H3>& v)
    {x=u.x-v.x;c.Sub(u.c,v.c);Base().Sub(u.Base(),v.Base());}

    template<class H2>
    void Fill_From(const MAT_MAT<TV,H2>& u)
    {Fill_From_Helper(x,u.x);c.Fill_From(u.c);Base().Fill_From(u.Base());}
};

template<class TV> void
Outer_Product_Helper(MAT_MAT<TV,MM_EMPTY>& r,const COL_VEC<TV,CV_EMPTY>& u) {}

template<class TV,class H,class G> void
Outer_Product_Helper(MAT_MAT<TV,H>& r,const COL_VEC<TV,G>& u)
{
    Outer_Product_Helper(r.x,u.x);
    Outer_Product_Helper(r.c,u.Base(),u.x);
    Outer_Product_Helper(r.Base(),u.Base());
}

template<class TV> void
Symmetric_Outer_Product_Helper(MAT_MAT<TV,MM_EMPTY>& r,const COL_VEC<TV,CV_EMPTY>& u,const COL_VEC<TV,CV_EMPTY>& v) {}

template<class TV,class H,class G,class G2> void
Symmetric_Outer_Product_Helper(MAT_MAT<TV,H>& r,const COL_VEC<TV,G>& u,const COL_VEC<TV,G2>& v)
{
    Symmetric_Outer_Product_Helper(r.x,u.x,v.x);

    COL_MAT<TV,CM_FLAGS<G::n-1,((G::mask&((1<<G::n)-1))>>1)*(G2::mask&1),0> > c0;
    COL_MAT<TV,CM_FLAGS<G::n-1,((G2::mask&((1<<G::n)-1))>>1)*(G::mask&1),0> > c1;
    Outer_Product_Helper(c0,u.Base(),v.x);
    Outer_Product_Helper(c1,v.Base(),u.x);
    r.c.Add(c0,c1);

    Symmetric_Outer_Product_Helper(r.Base(),u.Base(),v.Base());
}

template<class TV> void Symmetric_Times_Transpose(MAT_MAT<TV,MM_EMPTY>& z,const COL_MAT<TV,CM_EMPTY>& u,const COL_MAT<TV,CM_EMPTY>& v) {}

template<class TV,class H,class H2,class H3> void
Symmetric_Times_Transpose(MAT_MAT<TV,H>& z,const COL_MAT<TV,H2>& u,const COL_MAT<TV,H3>& v)
{
    z.x=Symmetric_Times_Transpose(u.x,v.x);
    MK_CM(TV,H2::Base()*H3::First()) c0;
    MK_CM(TV,H3::Base()*H2::First()) c1;
    c0.Times(u.Base(),v.x.Transposed());
    c1.Times(v.Base(),u.x.Transposed());
    z.c.Add(c0,c1);
    Symmetric_Times_Transpose(z.Base(),u.Base(),v.Base());
}

template<class TV> void Times_Self_Transpose(MAT_MAT<TV,MM_EMPTY>& z,const COL_MAT<TV,CM_EMPTY>& u) {}

template<class TV,class H,class H2> void Times_Self_Transpose(MAT_MAT<TV,H>& z,const COL_MAT<TV,H2>& u)
{
    z.x=Times_Self_Transpose(u.x);
    z.c.Times(u.Base(),u.x.Transposed());
    Times_Self_Transpose(z.Base(),u.Base());
}

}
}


#endif
