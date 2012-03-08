//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Symbolics/MULTIVARIATE_POLYNOMIAL.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MULTIVARIATE_POLYNOMIAL<TV>::
MULTIVARIATE_POLYNOMIAL()
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MULTIVARIATE_POLYNOMIAL<TV>::
~MULTIVARIATE_POLYNOMIAL()
{}
//#####################################################################
// Function Simplify
//#####################################################################
template<class TV> void MULTIVARIATE_POLYNOMIAL<TV>::
Simplify()
{
    terms.Coalesce();
    int k=0;
    for(int i=0;i<terms.m;i++) if(terms(i).coeff) terms(k++)=terms(i);
    terms.Resize(k);
}
//#####################################################################
// Function Max_Power
//#####################################################################
template<class TV> VECTOR<int,TV::m> MULTIVARIATE_POLYNOMIAL<TV>::
Max_Power() const
{
    TV_INT mp;
    for(int i=0;i<terms.m;i++) mp=mp.Componentwise_Max(mp,terms(i).power);
    return mp;
}
//#####################################################################
// Operator +=
//#####################################################################
template<class TV> MULTIVARIATE_POLYNOMIAL<TV>& MULTIVARIATE_POLYNOMIAL<TV>::
operator+= (const MULTIVARIATE_POLYNOMIAL& m)
{
    terms.Append_Elements(m.terms);
    Simplify();
    return *this;
}
//#####################################################################
// Operator -=
//#####################################################################
template<class TV> MULTIVARIATE_POLYNOMIAL<TV>& MULTIVARIATE_POLYNOMIAL<TV>::
operator-= (const MULTIVARIATE_POLYNOMIAL& m)
{
    int n=terms.m;
    terms.Append_Elements(m.terms);
    for(int i=n;i<terms.m;i++) terms(i).coeff=-terms(i).coeff;
    Simplify();
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class TV> MULTIVARIATE_POLYNOMIAL<TV>& MULTIVARIATE_POLYNOMIAL<TV>::
operator*= (const MULTIVARIATE_POLYNOMIAL& m)
{
    ARRAY<MULTIVARIATE_MONOMIAL<TV> > array;
    for(int i=0;i<terms.m;i++)
        for(int j=0;j<m.terms.m;j++)
            array.Append(MULTIVARIATE_MONOMIAL<TV>(terms(i).power+m.terms(j).power,terms(i).coeff*m.terms(j).coeff));
    array.Exchange(terms);
    Simplify();
    return *this;
}
//#####################################################################
// Function Differentiate
//#####################################################################
template<class TV> void MULTIVARIATE_POLYNOMIAL<TV>::
Differentiate(int v)
{
    int k=0;
    for(int i=0;i<terms.m;i++)
        if(terms(i).power(v)){
            terms(k)=terms(i);
            terms(k).coeff*=terms(k).power(v);
            terms(k++).power(v)--;}
    terms.Resize(k);
}
//#####################################################################
// Function Integrate
//#####################################################################
template<class TV> void MULTIVARIATE_POLYNOMIAL<TV>::
Integrate(int v)
{
    for(int i=0;i<terms.m;i++){
        terms(i).power(v)++;
        terms(i).coeff/=terms(i).power(v);}
}
//#####################################################################
// Function Negate_Variable
//#####################################################################
template<class TV> void MULTIVARIATE_POLYNOMIAL<TV>::
Negate_Variable(int v)
{
    for(int i=0;i<terms.m;i++)
        if(terms(i).power(v)&1)
            terms(i).coeff=-terms(i).coeff;
}
//#####################################################################
// Function Scale_Variables
//#####################################################################
template<class TV> void MULTIVARIATE_POLYNOMIAL<TV>::
Scale_Variables(const TV& scale)
{
    TV_INT max_power=Max_Power();
    assert(max_power.Max()<20);
    T powers[TV::m][20];
    for(int i=0;i<TV::m;i++){
        powers[i][0]=1;
        for(int j=1;j<=max_power(i);j++)
            powers[i][j]=powers[i][j-1]*scale(i);}
    for(int i=0;i<terms.m;i++)
        for(int j=0;j<TV::m;j++)
            terms(i).coeff*=powers[j][terms(i).power(j)];
}
//#####################################################################
// Function Shift_Variable
//#####################################################################
template<class TV> void MULTIVARIATE_POLYNOMIAL<TV>::
Shift_Variable(int v,T shift)
{
    int mx=Max_Power()(v);
    assert(mx<20);
    T table[20][20];
    table[0][0]=1;
    for(int i=1;i<=mx;i++){
        table[i][0]=1;
        table[i][i]=table[i-1][i-1]*shift;
        for(int j=1;j<mx;j++)
            table[i][j]=table[i-1][j-1]*shift+table[i-1][j];}
    for(int i=terms.m-1;i>=0;i--){
    for(int j=1;j<=terms(i).power(v);j++){
        terms.Append(MULTIVARIATE_MONOMIAL<TV>(terms(i).power,terms(i).coeff*table[terms(i).power(v)][j]));
        terms.Last().power(v)-=j;}}
    Simplify();
}
//#####################################################################
// Function Shift
//#####################################################################
template<class TV> void MULTIVARIATE_POLYNOMIAL<TV>::
Shift(const TV& shift)
{
    for(int i=0;i<TV::m;i++)
        Shift_Variable(i,shift(i));
}
//#####################################################################
// Function Exchange_Variables
//#####################################################################
template<class TV> void MULTIVARIATE_POLYNOMIAL<TV>::
Exchange_Variables(int u,int v)
{
    if(u==v) return;
    for(int i=0;i<terms.m;i++)
        exchange(terms(i).power(u),terms(i).power(v));
}
//#####################################################################
// Function Definite_Integral
//#####################################################################
template<class TV> typename TV::ELEMENT MULTIVARIATE_POLYNOMIAL<TV>::
Definite_Integral(RANGE<TV>& domain) const
{
    TV_INT max_power=Max_Power()+1;
    assert(max_power.Max()<20);
    T powers[2][TV::m][20];
    for(int i=0;i<TV::m;i++){
        powers[0][i][0]=1;
        powers[1][i][0]=1;
        for(int j=1;j<=max_power(i);j++){
            powers[0][i][j]=powers[0][i][j-1]*domain.min_corner(i);
            powers[1][i][j]=powers[1][i][j-1]*domain.max_corner(i);}}
    T s=0;
    for(int i=0;i<terms.m;i++){
        TV_INT p1=terms(i).power+1;
        T t=terms(i).coeff/p1.Product();
        for(int j=0;j<TV::m;j++)
            t*=powers[1][j][p1(j)]-powers[0][j][p1(j)];
        s+=t;}
    return s;
}
//#####################################################################
// Function Integrate_Over_Triangle
//#####################################################################
template<class TV> typename TV::ELEMENT MULTIVARIATE_POLYNOMIAL<TV>::
Integrate_Over_Triangle(const VECTOR<TV,3>& vertices) const
{
    
}
//#####################################################################
// Function operator<<
//#####################################################################
template<class TV> std::ostream& PhysBAM::
operator<< (std::ostream& o, const MULTIVARIATE_POLYNOMIAL<TV>& p)
{
    if(!p.terms.m) return o<<"(0)";
    o<<"(";
    for(int i=0;i<p.terms.m;i++){
        typename TV::ELEMENT c=p.terms(i).coeff;
        if(c<0){o<<(i>0?"-":"-");c=-c;}
        else o<<(i>0?"+":"");
        if(c!=1 || p.terms(i).power==VECTOR<int,TV::m>()) o<<c;

        for(int j=0;j<TV::m;j++)
            if(p.terms(i).power(j)>0){
                o<<"abcdefghijklmnopqrstuvwxyz"[j];
                if(p.terms(i).power(j)>1)
                    o<<'^'<<p.terms(i).power(j);}}
    return o<<")";
}
template class MULTIVARIATE_POLYNOMIAL<VECTOR<float,1> >;
template class MULTIVARIATE_POLYNOMIAL<VECTOR<float,2> >;
template class MULTIVARIATE_POLYNOMIAL<VECTOR<float,3> >;
template std::ostream& PhysBAM::operator<< <VECTOR<float,1> >(std::ostream&,MULTIVARIATE_POLYNOMIAL<VECTOR<float,1> > const&);
template std::ostream& PhysBAM::operator<< <VECTOR<float,2> >(std::ostream&,MULTIVARIATE_POLYNOMIAL<VECTOR<float,2> > const&);
template std::ostream& PhysBAM::operator<< <VECTOR<float,3> >(std::ostream&,MULTIVARIATE_POLYNOMIAL<VECTOR<float,3> > const&);
#ifndef COMPILATE_WITHOUT_DOUBLE_SUPPORT
template class MULTIVARIATE_POLYNOMIAL<VECTOR<double,1> >;
template class MULTIVARIATE_POLYNOMIAL<VECTOR<double,2> >;
template class MULTIVARIATE_POLYNOMIAL<VECTOR<double,3> >;
template std::ostream& PhysBAM::operator<< <VECTOR<double,1> >(std::ostream&,MULTIVARIATE_POLYNOMIAL<VECTOR<double,1> > const&);
template std::ostream& PhysBAM::operator<< <VECTOR<double,2> >(std::ostream&,MULTIVARIATE_POLYNOMIAL<VECTOR<double,2> > const&);
template std::ostream& PhysBAM::operator<< <VECTOR<double,3> >(std::ostream&,MULTIVARIATE_POLYNOMIAL<VECTOR<double,3> > const&);
#endif
