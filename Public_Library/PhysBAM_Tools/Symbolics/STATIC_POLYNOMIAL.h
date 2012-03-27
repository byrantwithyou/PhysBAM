//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STATIC_POLYNOMIAL
//#####################################################################
#ifndef __STATIC_POLYNOMIAL__
#define __STATIC_POLYNOMIAL__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/UNIFORM_ARRAY_ITERATOR.h>
#include <PhysBAM_Tools/Vectors/STATIC_TENSOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T,int rank,int d>
struct STATIC_POLYNOMIAL
{
    typedef VECTOR<int,rank> TV_INT;
    typedef VECTOR<T,rank> TV;
    TV_INT size;
    STATIC_TENSOR<T,rank,d+1> terms;

    void Set_Term(const TV_INT& power,T x)
    {size=size.Componentwise_Max(power);terms(power)=x;}

    template<int d2>
    STATIC_POLYNOMIAL<T,rank,(d>d2?d:d2)> operator+ (const STATIC_POLYNOMIAL<T,rank,d2>& p)
    {STATIC_POLYNOMIAL<T,rank,(d>d2?d:d2)> r;r+=*this;r+=p;return r;}

    template<int d2>
    STATIC_POLYNOMIAL<T,rank,(d>d2?d:d2)> operator- (const STATIC_POLYNOMIAL<T,rank,d2>& p)
    {STATIC_POLYNOMIAL<T,rank,(d>d2?d:d2)> r;r+=*this;r-=p;return r;}

    template<int d2>
    STATIC_POLYNOMIAL<T,rank,(d+d2)> operator* (const STATIC_POLYNOMIAL<T,rank,d2>& p)
    {STATIC_POLYNOMIAL<T,rank,(d+d2)> r;r.Multiply(*this,p,false);return r;}

    void Compress_Size()
    {
        RANGE<TV_INT> range(TV_INT(),size+1);
        size.Fill(0);
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next())
            if(terms(it.index))
                size=size.Componentwise_Max(it.index);
    }

    template<int d2,int d3> // Result must fit; no aliasing
    void Multiply(const STATIC_POLYNOMIAL<T,rank,d2>& a, const STATIC_POLYNOMIAL<T,rank,d3>& b, bool need_clear=true)
    {
        if(&a==this){
            STATIC_POLYNOMIAL<T,rank,d3> copy(a);
            if(&b==this) return Multiply(copy,copy);
            return Multiply(copy,b);}
        if(&b==this){
            STATIC_POLYNOMIAL<T,rank,d3> copy(b);
            return Multiply(a,copy);}

        size=a.size+b.size;
        PHYSBAM_ASSERT(size.All_Less_Equal(TV_INT()+d));
        RANGE<TV_INT> range(TV_INT(),a.size+1),range2(TV_INT(),b.size+1);
        if(need_clear){
            RANGE<TV_INT> range(TV_INT(),size+1);
            for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next()) terms(it.index)=0;}
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next()){
            T x=a.terms(it.index);
            if(!x) continue;
            for(UNIFORM_ARRAY_ITERATOR<rank> it2(range2);it2.Valid();it2.Next())
                terms(it.index+it2.index)+=x*b.terms(it2.index);}
    }

    template<int d2> // Result must fit
    STATIC_POLYNOMIAL& operator+=(const STATIC_POLYNOMIAL<T,rank,d2>& a)
    {
        PHYSBAM_ASSERT(a.size.All_Less_Equal(TV_INT()+d));
        RANGE<TV_INT> range(TV_INT(),a.size+1);
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next())
            terms(it.index)+=a.terms(it.index);
        size=size.Componentwise_Max(a.size);
        return *this;
    }

    template<int d2> // Result must fit
    STATIC_POLYNOMIAL& operator-=(const STATIC_POLYNOMIAL<T,rank,d2>& a)
    {
        PHYSBAM_ASSERT(a.size.All_Less_Equal(TV_INT()+d));
        RANGE<TV_INT> range(TV_INT(),a.size+1);
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next())
            terms(it.index)-=a.terms(it.index);
        size=size.Componentwise_Max(a.size);
        return *this;
    }

    STATIC_POLYNOMIAL Differentiate(int v) const
    {
        STATIC_POLYNOMIAL r;
        RANGE<TV_INT> range(TV_INT::Axis_Vector(v),size);
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next())
            r.terms(it.index-TV_INT::Axis_Vector(v))=terms(it.index)*it.index(v);
        r.size=size;
        r.size(v)--;
        return r;
    }

    STATIC_POLYNOMIAL<T,rank,d+1> Integrate(int v) const
    {
        STATIC_POLYNOMIAL<T,rank,d+1> r;
        RANGE<TV_INT> range(TV_INT(),size+1);
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next())
            r.terms(it.index+TV_INT::Axis_Vector(v))=terms(it.index)/(it.index(v)+1);
        r.size=size;
        r.size(v)++;
        return r;
    }

    void Negate_Variable(int v)
    {
        RANGE<TV_INT> range(TV_INT(),size+1);
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next())
            if(it.index(v)&1)
                terms(it.index)=-terms(it.index);
    }

    void Scale_Variables(const TV& scale)
    {for(int i=0;i<rank;i++) Scale_Variables(i,scale(i));}

    void Scale_Variables(int v,T scale)
    {
        T powers[d+1];
        powers[0]=1;
        for(int j=1;j<=size(v);j++)
            powers[j]=powers[j-1]*scale;
        RANGE<TV_INT> range(TV_INT(),size+1);
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next())
            terms(it.index)*=powers[it.index(v)];
    }

    void Shift(const TV& scale)
    {for(int i=0;i<rank;i++) Shift_Variables(i,scale(i));}

    void Shift_Variable(int v,T shift)
    {
        T table[d+1][d+1];
        table[0][0]=1;
        for(int i=1;i<=size(v);i++){
            table[i][0]=1;
            table[i][i]=table[i-1][i-1]*shift;
            for(int j=1;j<i;j++)
                table[i][j]=table[i-1][j-1]*shift+table[i-1][j];}
        RANGE<TV_INT> range(TV_INT(),size+1);
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next()){
            T x=terms(it.index);
            TV_INT index=it.index;
            for(int i=1;i<=it.index(v);i++){
                index(v)--;
                terms(index)+=x*table[it.index(v)][i];}}
    }

    void Exchange_Variables(int u,int v)
    {
        RANGE<TV_INT> range(TV_INT(),size+1);
        range.max_corner(u)=1;
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next())
            for(int j=0;j<it.index(v);j++){
                TV_INT index=it.index;
                index(u)=j;
                TV_INT index2=index;
                exchange(index2(u),index2(v));
                exchange(terms(index),terms(index2));}
        exchange(size(u),size(v));
    }

    T Definite_Integral(const RANGE<TV>& domain) const
    {
        T powers[2][rank][d+2];
        for(int i=0;i<rank;i++){
            powers[0][i][0]=1;
            powers[1][i][0]=1;
            for(int j=1;j<=size(i)+1;j++){
                powers[0][i][j]=powers[0][i][j-1]*domain.min_corner(i);
                powers[1][i][j]=powers[1][i][j-1]*domain.max_corner(i);}}
        T s=0;
        RANGE<TV_INT> range(TV_INT(),size+1);
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next()){
            TV_INT index=it.index+1;
            T t=terms(it.index)/index.Product();
            for(int j=0;j<rank;j++)
                t*=powers[1][j][index(j)]-powers[0][j][index(j)];
            s+=t;}
        return s;
    }

    T Value(const TV& point) const
    {
        T powers[rank][d+1];
        for(int i=0;i<rank;i++){
            powers[i][0]=1;
            for(int j=1;j<=size(i);j++)
                powers[i][j]=powers[i][j-1]*point(i);}
        T s=0;
        RANGE<TV_INT> range(TV_INT(),size+1);
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next()){
            T t=terms(it.index);
            for(int j=0;j<rank;j++)
                t*=powers[j][it.index(j)];
            s+=t;}
        return s;
    }

    // T Integrate_Over_Primitive(const VECTOR<TV,3>& vertices) const; // triangle
    // T Integrate_Over_Primitive(const VECTOR<TV,2>& vertices) const; // interval
    
    T Integrate_Over_Primitive(const VECTOR<TV,2>& vertices) const
    {
        STATIC_POLYNOMIAL copy(*this);
        
        copy.Shift(vertices(0));
        TV a=vertices(1)-vertices(0);
        
        T table[rank][d+1];
        for(int v=0;v<rank;v++){
            table[v][0]=1;
            for(int i=1;i<=size(v);i++) table[v][i]=table[v][i-1]*a(v);}
        
        typedef VECTOR<T,1> TV1;
        typedef VECTOR<int,1> TV_INT1;
        
        STATIC_POLYNOMIAL<T,1,d*rank> barycentric;
        RANGE<TV_INT> range(TV_INT(),size+1);
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next()){
            T scale=1;
            for(int v=0;v<rank;v++) scale*=table[v][it.index(v)];
            barycentric.terms(TV_INT1(it.index.Sum()))=terms(it.index)*scale;}
        barycentric.size=size.Sum();

        T integral=0;
        for(int i=0;i<barycentric.size(0);i++)
             integral+=barycentric.terms(TV_INT1(i))/(i+1);

        return integral*a.Magnitude();
    }

    T Integrate_Over_Primitive(const VECTOR<TV,3>& vertices) const
    {
        STATIC_POLYNOMIAL copy(*this);

        copy.Shift(vertices(0));
        TV a=vertices(1)-vertices(0);
        TV b=vertices(2)-vertices(0);
    
        T table[rank][d+1][d+1];
        for(int v=0;v<rank;v++){
            table[v][0][0]=1;
            for(int i=1;i<=size(v);i++){
                table[v][i][0]=table[v][i-1][0]*a(v);
                table[v][i][i]=table[v][i-1][i-1]*b(v);
                for(int j=1;j<i;j++)
                    table[v][i][j]=table[v][i-1][j-1]*b(v)+table[v][i-1][j]*a(v);}}
    
        typedef VECTOR<T,2> TV2;
        typedef VECTOR<int,2> TV_INT2;

        STATIC_POLYNOMIAL<T,2,d*rank> barycentric;
        RANGE<TV_INT> range(TV_INT(),size+1);
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next()){
            STATIC_POLYNOMIAL<T,2,0> monomial;
            monomial.Set_Term(TV_INT2(),copy.terms(it.index));
            for(int v=0;v<rank;v++){
                STATIC_POLYNOMIAL<T,2,d> variable;
                int n=it.index(v);
                for(int j=0;j<=n;j++)
                    variable.terms(TV_INT2(j,n-j))=table[v][n][j];
                variable.size=TV_INT2(n,n);
                monomial.Multiply(monomial,variable);}
            barycentric+=monomial;}
    
        int max_factorial=barycentric.Max_Power_Sum()+2;
        T factorial[20];
        factorial[0]=(T)1;
        for(int i=1;i<=max_factorial;i++) factorial[i]=factorial[i-1]*i;        

        T integral=0;
        for(int i=0;i<barycentric.terms.m;i++){
            TV_INT2 power=barycentric.terms(i).power;
            integral+=barycentric.terms(i).coeff*factorial[power(0)]*factorial[power(1)]/factorial[power(0)+power(1)+2];}

        return integral*TV::Cross_Product(a,b).Magnitude();
    }
};

template<class T,int rank,int d>
std::ostream& operator<< (std::ostream& o, const STATIC_POLYNOMIAL<T,rank,d>& p)
{
    bool first=true;
    o<<"(";
    RANGE<VECTOR<int,rank> > range(VECTOR<int,rank>(),p.size+1);
    for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next()){
        T c=p.terms(it.index);
        if(c==0) continue;
        if(c<0){o<<"-";c=-c;}
        else o<<(first?"":"+");
        if(c!=1 || it.index==VECTOR<int,rank>()) o<<c;
        first=false;

        for(int j=0;j<rank;j++)
            if(it.index(j)>0){
                o<<"abcdefghijklmnopqrstuvwxyz"[j];
                if(it.index(j)>1)
                    o<<'^'<<it.index(j);}}

    if(first) o<<"0";
    return o<<")";
}
}
#endif
