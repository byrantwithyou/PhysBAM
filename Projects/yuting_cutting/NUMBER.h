//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYREGHT.txt.
//#####################################################################
// Class NUMBER
//#####################################################################
#ifndef __NUMBER__
#define __NUMBER__

namespace PhysBAM{
struct NUMBER;
inline NUMBER Magnitude_Squared(const NUMBER& a);
}

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Log/LOG.h>
#include <Tools/Utilities/TYPE_UTILITIES.h>
#include <cmath>
#include <string>
using std::abs;
namespace PhysBAM{

#define P(x) P_(x,#x)
#define NOE(x) NOE_(x,#x)
inline void P_(...){}
void P_(const NUMBER& n,const char* s);
inline void NOE_(...){}
void NOE_(NUMBER& n,const char* s);

struct NUMBER_COMP_HELPER
{
    ARRAY<ARRAY<bool> > number_tests_todo;
    ARRAY<bool> number_tests;
    int number_tests_used;
    ARRAY<std::string> comparisions;

    void Reset()
    {
        number_tests_todo.Remove_All();
        number_tests_todo.Resize(1);
        comparisions.Remove_All();
    }

    void Start_Test()
    {
        comparisions.Remove_All();
        number_tests=number_tests_todo.Pop();
        number_tests_used=0;
    }
};
extern NUMBER_COMP_HELPER number_comp_helper;

struct BOOL_NUMBER
{
    std::string val;

    operator bool() const
    {
        if(number_comp_helper.number_tests_used==number_comp_helper.number_tests.m){
            number_comp_helper.number_tests.Append(true);
            number_comp_helper.number_tests_todo.Append(number_comp_helper.number_tests);
            number_comp_helper.number_tests.Last()=false;}
        if(number_comp_helper.number_tests(number_comp_helper.number_tests_used++)){
            number_comp_helper.comparisions.Append(val);
            return true;}
        number_comp_helper.comparisions.Append("NOT("+val+")");
        return false;
    }
};

struct NUMBER
{
    std::string val;

    bool valid;
    bool pos;
    bool rel;
    double eps,eps2; // relative to value (if rel) or max mag (otherwise)
    double c,p; // max mag is c * L^p;
    bool inv;

    NUMBER(){val="0";pos=false;valid=false;rel=false;inv=false;eps=eps2=c=p=0;}
    NUMBER(int i){val=std::to_string(i);pos=(i>=0);valid=true;inv=false;rel=true;eps=eps2=p=0;c=fabs(i);}
    NUMBER(long long i){val=std::to_string(i);pos=(i>=0);valid=true;inv=false;rel=true;eps=eps2=p=0;c=fabs(i);}
    NUMBER(float i){val=std::to_string(i);pos=(i>=0);valid=true;inv=false;rel=true;eps=eps2=p=0;c=fabs(i);}
    NUMBER(double i){val=std::to_string(i);pos=(i>=0);valid=true;inv=false;rel=true;eps=eps2=p=0;c=fabs(i);}

    NUMBER& operator=(double i){return *this=NUMBER(i);}

    NUMBER operator- () const
    {
        NUMBER n=*this;
        n.val = "NEG(" + val + ")";
        n.pos=false;
        if(inv) n.valid=false;
        return n;
    }
    NUMBER operator* (const NUMBER& a) const
    {
        NUMBER n;
        n.val = "MUL(" + val + "," + a.val + ")";
        if(valid && a.valid)
        {
            n.valid=true;
            n.pos=(val==a.val || (pos && a.pos));
            if(rel && a.rel) n.rel=true;
            n.eps=1+eps+a.eps;
            n.eps2=(1+eps)*(1+a.eps)+eps2+a.eps2;
            n.c=c*a.c;
            n.p=p+a.p;
        }
        if(inv || a.inv) n.c=1;
        return n;
    }
    NUMBER operator/ (const NUMBER& a) const
    {
        NUMBER n;
        n.val = "DIV(" + val + "," + a.val + ")";
        if(valid && a.valid && rel && a.rel)
        {
            n.valid=true;
            n.rel=true;
            n.pos=(pos && a.pos);
            n.eps=1+eps+a.eps;
            n.eps2=(1+eps)*(1+a.eps)+eps2+a.eps2+a.eps*a.eps;
            n.c=c/a.c;
            n.p=p-a.p;
        }
        if(inv || a.inv) n.valid=false;
        n.inv=true;
        return n;
    }
    NUMBER operator+ (const NUMBER& a) const
    {
        NUMBER n;
        n.val = "ADD(" + val + "," + a.val + ")";
        if(valid && a.valid)
        {
            n.valid=true;
            n.pos=(pos && a.pos);
            if(p!=a.p) n.valid=false;
            n.c=c+a.c;
            n.p=p;
            n.rel=n.pos && rel && a.rel;
            if(n.rel)
            {
                n.eps=1+max(eps,a.eps);
                n.eps2=1+max(eps+eps2,a.eps+a.eps2);
            }
            else
            {
                n.eps=1+(c*eps+a.c*a.eps)/n.c;
                n.eps2=(c*(eps+eps2)+a.c*(a.eps+a.eps2)+1)/n.c;
            }
        }
        if(inv || a.inv) n.valid=false;
        return n;
    }
    NUMBER operator- (const NUMBER& a) const
    {
        NUMBER n;
        n.val = "SUB(" + val + "," + a.val + ")";
        if(valid && a.valid)
        {
            n.valid=true;
            n.pos=false;
            if(p!=a.p) n.valid=false;
            n.c=c+a.c;
            n.p=p;
            if(!eps && !eps2 && !a.eps && !a.eps2)
            {
                n.rel=true;
                n.eps=1;
                n.eps2=0;
            }
            else
            {
                n.rel=false;
                n.eps=1+(c*eps+a.c*a.eps)/n.c;
                n.eps2=(c*(eps+eps2)+a.c*(a.eps+a.eps2)+1)/n.c;
            }
        }
        if(inv || a.inv) n.valid=false;
        return n;
    }

    NUMBER& operator*= (const NUMBER& a) {return *this = *this * a;}
    NUMBER& operator/= (const NUMBER& a) {return *this = *this / a;}
    NUMBER& operator+= (const NUMBER& a) {return *this = *this + a;}
    NUMBER& operator-= (const NUMBER& a) {return *this = *this - a;}

    BOOL_NUMBER operator< (const NUMBER& a) const {BOOL_NUMBER n;n.val = "LT(" + val + "," + a.val + ")";return n;}
    BOOL_NUMBER operator> (const NUMBER& a) const {BOOL_NUMBER n;n.val = "GT(" + val + "," + a.val + ")";return n;}
    BOOL_NUMBER operator<= (const NUMBER& a) const {BOOL_NUMBER n;n.val = "LE(" + val + "," + a.val + ")";return n;}
    BOOL_NUMBER operator>= (const NUMBER& a) const {BOOL_NUMBER n;n.val = "GE(" + val + "," + a.val + ")";return n;}
    BOOL_NUMBER operator== (const NUMBER& a) const {BOOL_NUMBER n;n.val = "EQ(" + val + "," + a.val + ")";return n;}
    BOOL_NUMBER operator!= (const NUMBER& a) const {BOOL_NUMBER n;n.val = "NE(" + val + "," + a.val + ")";return n;}

    template<class RW> void Read(std::istream& input){}
    template<class RW> void Write(std::ostream& output) const{}

    explicit operator bool() const
    {
        if(number_comp_helper.number_tests_used==number_comp_helper.number_tests.m){
            number_comp_helper.number_tests.Append(true);
            number_comp_helper.number_tests_todo.Append(number_comp_helper.number_tests);
            number_comp_helper.number_tests.Last()=false;}
        if(number_comp_helper.number_tests(number_comp_helper.number_tests_used++)){
            number_comp_helper.comparisions.Append(val);
            return true;}
        number_comp_helper.comparisions.Append("NOT("+val+")");
        return false;
    }

    explicit operator int() const
    {
        number_comp_helper.comparisions.Append("Forced to int: "+val);
        return 0;
    }
};

inline void P_(const NUMBER& n,const char* s)
{
    if(n.valid) LOG::printf("%s : %s%s err: %.16g %.16g  mag: %.16g L^%.16g (%.16g %.16g)  :  %P\n",s,n.pos?"(+) ":"",n.rel?"rel ":"",n.eps,n.eps2,n.c,n.p,n.eps*n.c,n.eps2*n.c,n.val);
    else LOG::printf("%s : invalid  :  %P\n",s,n.val);
}
inline void NOE_(NUMBER& n,const char* s)
{
    n.valid=true;
    n.rel=true;
    n.eps=0;
    n.eps2=0;
    n.val=s;
}

template<> struct IS_SCALAR<NUMBER> {static const int value=true;};
template<> struct IS_FLOAT_OR_DOUBLE<NUMBER> {static const int value=true;};


template<class T> typename ENABLE_IF<is_same<NUMBER,T>::value,NUMBER>::TYPE  abs(const T& v)
{
    NUMBER n;
    n.val="ABS("+v.val+")";
    if(v.valid)
    {
        n.valid=true;
        n.pos=true;
        n.c=v.c;
        n.p=v.p;
        n.rel=v.rel;
        n.eps=v.eps;
        n.eps2=v.eps2;
    }
    if(v.inv) n.valid=false;
    return n;
}
template<class T> typename ENABLE_IF<is_same<NUMBER,T>::value,NUMBER>::TYPE  floor(const T& v){NUMBER n;n.val="FLOOR("+v.val+")";return n;}
template<class T> typename ENABLE_IF<is_same<NUMBER,T>::value,NUMBER>::TYPE  ceil(const T& v){NUMBER n;n.val="CEIL("+v.val+")";return n;}
template<class T> typename ENABLE_IF<is_same<NUMBER,T>::value,NUMBER>::TYPE  rint(const T& v){NUMBER n;n.val="RINT("+v.val+")";return n;}
template<class T> typename ENABLE_IF<is_same<NUMBER,T>::value,NUMBER>::TYPE  exp(const T& v){NUMBER n;n.val="EXP("+v.val+")";return n;}
template<class T> typename ENABLE_IF<is_same<NUMBER,T>::value,NUMBER>::TYPE  log(const T& v){NUMBER n;n.val="LOG("+v.val+")";return n;}
template<class T> typename ENABLE_IF<is_same<NUMBER,T>::value,NUMBER>::TYPE  sin(const T& v){NUMBER n;n.val="SIN("+v.val+")";return n;}
template<class T> typename ENABLE_IF<is_same<NUMBER,T>::value,NUMBER>::TYPE  cos(const T& v){NUMBER n;n.val="COS("+v.val+")";return n;}
template<class T> typename ENABLE_IF<is_same<NUMBER,T>::value,NUMBER>::TYPE  tan(const T& v){NUMBER n;n.val="TAN("+v.val+")";return n;}
template<class T> typename ENABLE_IF<is_same<NUMBER,T>::value,NUMBER>::TYPE  asin(const T& v){NUMBER n;n.val="ASIN("+v.val+")";return n;}
template<class T> typename ENABLE_IF<is_same<NUMBER,T>::value,NUMBER>::TYPE  acos(const T& v){NUMBER n;n.val="ACOS("+v.val+")";return n;}
template<class T> typename ENABLE_IF<is_same<NUMBER,T>::value,NUMBER>::TYPE  atan(const T& v){NUMBER n;n.val="ATAN("+v.val+")";return n;}
template<class T> typename ENABLE_IF<is_same<NUMBER,T>::value,NUMBER>::TYPE  sqrt(const T& v)
{
    NUMBER n;
    n.val="SQRT("+v.val+")";
    if(v.valid && v.pos && v.rel)
    {
        n.valid=true;
        n.pos=true;
        n.c=std::sqrt(v.c);
        n.p=v.p/2;
        n.rel=v.rel;
        n.eps=1+v.eps/2;
        n.eps2=1+v.eps2/2+fabs(v.eps/2-v.eps*v.eps/8);
    }
    if(v.inv) n.valid=false;
    return n;
}
template<class T> typename ENABLE_IF<is_same<NUMBER,T>::value,NUMBER>::TYPE  Inverse(const T& v){NUMBER n;n.val="INV("+v.val+")";return n;}
template<class T,class U> typename ENABLE_IF<is_same<NUMBER,T>::value||is_same<NUMBER,U>::value,NUMBER>::TYPE  min(const T& a,const U& b){NUMBER n;n.val="MIN("+a.val+","+b.val+")";return n;}
template<class T,class U> typename ENABLE_IF<is_same<NUMBER,T>::value||is_same<NUMBER,U>::value,NUMBER>::TYPE  max(const T& a,const U& b){NUMBER n;n.val="MAX("+a.val+","+b.val+")";return n;}
template<class T,class U> typename ENABLE_IF<is_same<NUMBER,T>::value||is_same<NUMBER,U>::value,NUMBER>::TYPE pow(const T& a,const U& b){NUMBER n;n.val="POW("+a.val+","+b.val+")";return n;}
template<class T,class U> typename ENABLE_IF<is_same<NUMBER,T>::value||is_same<NUMBER,U>::value,NUMBER>::TYPE atan2(const T& a,const U& b){NUMBER n;n.val="ATAN2("+a.val+","+b.val+")";return n;}
inline NUMBER clamp(const NUMBER& v,const NUMBER& vmin,const NUMBER& vmax){NUMBER n;n.val="MAX("+vmin.val+"MIN("+vmax.val+","+v.val+"))";return n;}

inline NUMBER operator* (double a,const NUMBER& b) {return NUMBER(a)*b;}
inline NUMBER operator/ (double a,const NUMBER& b) {return NUMBER(a)/b;}
inline NUMBER operator+ (double a,const NUMBER& b) {return NUMBER(a)+b;}
inline NUMBER operator- (double a,const NUMBER& b) {return NUMBER(a)-b;}

inline BOOL_NUMBER operator< (double a,const NUMBER& b) {return NUMBER(a)<b;}
inline BOOL_NUMBER operator> (double a,const NUMBER& b) {return NUMBER(a)>b;}
inline BOOL_NUMBER operator<= (double a,const NUMBER& b) {return NUMBER(a)<=b;}
inline BOOL_NUMBER operator>= (double a,const NUMBER& b) {return NUMBER(a)>=b;}
inline BOOL_NUMBER operator== (double a,const NUMBER& b) {return NUMBER(a)==b;}
inline BOOL_NUMBER operator!= (double a,const NUMBER& b) {return NUMBER(a)!=b;}

template<class T> typename ENABLE_IF<is_same<NUMBER,T>::value,NUMBER>::TYPE Magnitude_Squared(const T& a){return a*a;}
inline NUMBER Magnitude_Squared(const NUMBER& a){return a*a;}

inline std::ostream& operator<<(std::ostream& output,const NUMBER& n)
{output<<"("<<n.val<<")";return output;}

inline std::istream& operator>>(std::istream& input,NUMBER& box)
{return input;}

using std::atan;
using std::acos;
using std::asin;
using std::tan;
using std::cos;
using std::sin;
using std::abs;
using std::floor;
using std::ceil;
using std::exp;
using std::log;
using std::sin;
using std::cos;
using std::tan;
using std::asin;
using std::acos;
using std::atan;
using std::sqrt;
using std::pow;
using std::atan2;

}
#endif
