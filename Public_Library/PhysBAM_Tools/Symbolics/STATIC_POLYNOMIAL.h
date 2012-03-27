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

    STATIC_POLYNOMIAL();
    ~STATIC_POLYNOMIAL();

    template<int d2>
    STATIC_POLYNOMIAL<T,rank,(d>d2?d:d2)> operator+ (const STATIC_POLYNOMIAL<T,rank,d2>& p)
    {
        STATIC_POLYNOMIAL<T,rank,(d>d2?d:d2)> r;
        RANGE<TV_INT> range(TV_INT(),size+1),range2(TV_INT(),p.size+1);
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next())
            r.terms(it.index)+=terms(it.index);
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range2);it.Valid();it.Next())
            r.terms(it.index)+=p.terms(it.index);
        r.size=size.Componentwise_Max(p.size);
        return r;
    }

    template<int d2>
    STATIC_POLYNOMIAL<T,rank,(d>d2?d:d2)> operator- (const STATIC_POLYNOMIAL<T,rank,d2>& p)
    {
        STATIC_POLYNOMIAL<T,rank,(d>d2?d:d2)> r;
        RANGE<TV_INT> range(TV_INT(),size+1),range2(TV_INT(),p.size+1);
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next())
            r.terms(it.index)+=terms(it.index);
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range2);it.Valid();it.Next())
            r.terms(it.index)-=p.terms(it.index);
        r.size=size.Componentwise_Max(p.size);
        return r;
    }

    template<int d2>
    STATIC_POLYNOMIAL<T,rank,(d+d2)> operator* (const STATIC_POLYNOMIAL<T,rank,d2>& p)
    {
        STATIC_POLYNOMIAL<T,rank,(d+d2)> r;
        RANGE<TV_INT> range(TV_INT(),size+1),range2(TV_INT(),p.size+1);
        for(UNIFORM_ARRAY_ITERATOR<rank> it(range);it.Valid();it.Next())
            for(UNIFORM_ARRAY_ITERATOR<rank> it2(range2);it2.Valid();it2.Next())
                r.terms(it.index+it2.index)+=terms(it.index)*p.terms(it2.index);
        r.size=size+p.size;
        return r;
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
};

//template<class TV>
//std::ostream& operator<< (std::ostream& o, const STATIC_POLYNOMIAL<TV>& p);
}
#endif
