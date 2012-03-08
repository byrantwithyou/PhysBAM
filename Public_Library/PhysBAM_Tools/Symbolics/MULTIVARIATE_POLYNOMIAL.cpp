//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
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
// Function Max_Power_Sum
//#####################################################################
template<class TV> int MULTIVARIATE_POLYNOMIAL<TV>::
Max_Power_Sum() const
{
    int mps=0;
    for(int i=0;i<terms.m;i++) mps=max(mps,terms(i).Power_Sum());
    return mps;
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
// Function Shift_Variable
//#####################################################################
template<class TV> void MULTIVARIATE_POLYNOMIAL<TV>::
Shear(int v,int w,T shift) // v -> v + shift * w
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
            terms.Last().power(v)-=j;
            terms.Last().power(w)+=j;}}
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
Definite_Integral(const RANGE<TV>& domain) const
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
// Function Integrate_Over_Primitive
//#####################################################################
template<class TV> typename TV::ELEMENT MULTIVARIATE_POLYNOMIAL<TV>::
Integrate_Over_Primitive(const VECTOR<TV,3>& vertices) const
{
    MULTIVARIATE_POLYNOMIAL<TV> copy=(*this);

    copy.Shift(vertices(0));
    TV a=vertices(1)-vertices(0);
    TV b=vertices(1)-vertices(0);
    
    TV_INT mp=copy.Max_Power();
    assert(mp.Max()<20);
    T table[TV::m][20][20];
    for(int v=0;v<TV::m;v++){
        table[v][0][0]=1;
        for(int i=1;i<=mp(v);i++){
            table[v][i][0]=table[v][i-1][0]*a(v);
            table[v][i][i]=table[v][i-1][i-1]*b(v);
            for(int j=1;j<mp(v);j++)
                table[v][i][j]=table[v][i-1][j-1]*b(v)+table[v][i-1][j]*a(v);}}
    
    typedef VECTOR<T,2> TV2;
    typedef VECTOR<int,2> TV_INT2;

    MULTIVARIATE_POLYNOMIAL<TV2> barycentric;
    for(int i=0;i<copy.terms.m;i++){
        MULTIVARIATE_POLYNOMIAL<TV2> monomial;
        monomial.terms.Append(MULTIVARIATE_MONOMIAL<TV2>(TV_INT2(0,0),copy.terms(i).coeff));
        for(int v=0;v<TV::m;v++){
            MULTIVARIATE_POLYNOMIAL<TV2> variable;
            int n=copy.terms(i).power(v);
            for(int j=0;j<=n;j++){
                variable.terms.Append(MULTIVARIATE_MONOMIAL<TV2>(TV_INT2(j,n-j),table[v][n][j]));}
            monomial*=variable;}
        barycentric+=monomial;}
    
    int max_factorial=barycentric.Max_Power_Sum()+2;
    assert(max_factorial<20);
    T factorial[20];
    factorial[0]=(T)1;
    for(int i=1;i<=max_factorial;i++) factorial[i]=factorial[i-1]*i;        

    T integral=0;
    for(int i=0;i<barycentric.terms.m;i++){
        TV_INT2 power=barycentric.terms(i).power;
        integral+=barycentric.terms(i).coeff*factorial[power(0)]*factorial[power(1)]/factorial[power(0)+power(1)+2];}

    return integral*TV::Cross_Product(a,b).Magnitude();
}
//#####################################################################
// Function Integrate_Over_Primitive
//#####################################################################
template<class TV> typename TV::ELEMENT MULTIVARIATE_POLYNOMIAL<TV>::
Integrate_Over_Primitive(const VECTOR<TV,2>& vertices) const
{
    MULTIVARIATE_POLYNOMIAL<TV> copy=(*this);

    copy.Shift(vertices(0));
    TV a=vertices(1)-vertices(0);

    TV_INT mp=copy.Max_Power();
    assert(mp.Max()<20);
    T table[TV::m][20];
    for(int v=0;v<TV::m;v++){
        table[v][0]=1;
        for(int i=1;i<=mp(v);i++) table[v][i]=table[v][i-1]*a(v);}
    
    typedef VECTOR<T,1> TV1;
    typedef VECTOR<int,1> TV_INT1;

    MULTIVARIATE_POLYNOMIAL<TV1> barycentric;
    for(int i=0;i<copy.terms.m;i++){
        TV_INT power=copy.terms(i).power;
        int scale=1;
        for(int v=0;v<TV::m;v++) scale*=table[v][power(v)];
        barycentric.terms.Append(MULTIVARIATE_MONOMIAL<TV1>(TV_INT1(copy.terms(i).Power_Sum()),copy.terms(i).coeff*scale));}
    barycentric.Simplify();
    
    T integral=0;
    for(int i=0;i<barycentric.terms.m;i++){
        TV_INT1 power=barycentric.terms(i).power;
        integral+=barycentric.terms(i).coeff/(power(0)+1);}

    return integral*a.Magnitude();
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
