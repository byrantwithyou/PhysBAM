//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STATIC_POLYNOMIAL
//#####################################################################
#ifndef __STATIC_POLYNOMIAL__
#define __STATIC_POLYNOMIAL__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Math_Tools/RANGE_ITERATOR.h>
#include <Tools/Vectors/STATIC_TENSOR.h>
#include <Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T,int rank,int d> struct STATIC_POLYNOMIAL;

template<class T,int rank,int d>
bool Same_Polynomial(const STATIC_POLYNOMIAL<T,rank,d>& a,const STATIC_POLYNOMIAL<T,rank,d>& b) {return &a==&b;}

template<class T,int rank,int d,int d2>
bool Same_Polynomial(const STATIC_POLYNOMIAL<T,rank,d>& a,const STATIC_POLYNOMIAL<T,rank,d2>& b) {return false;}

struct QUADRITURE_RULE_TRI
{
    bool use_center;
    double center_weight;
    int number3;
    double location3[3];
    double weight3[3];
    int number6;
    double location6a[1],location6b[1];
    double weight6[1];
};

const QUADRITURE_RULE_TRI quadriture_rule_tri[9]={
    {1,1,0},
    {1,1,0},
    {0,0,1,{0.5},{1./3}},
    {1,-9./16,1,{0.2},{25./48}},
    {0,0,2,{0.44594849091596488,0.0915762135097707435},{0.223381589678011466,0.109951743655321868}},
    {1,0.225,2,{0.47014206410511511,0.10128650732345632},{0.132394152788506181,0.125939180544827152}},
    {0,0,2,{0.219429982549781601,0.480137964112216},{0.1713331241529814,0.0807310895930310},1,{0.0193717243612407884},{0.141619015923968181},{0.04063455979366063}},
    {1,0.197270749512171954,2,{0.0578921421877734927,1.10297642785225161},{0.0444948443100991969,0.754643787840106e-6},1,{0.0712675584422491723},{0.298212365719518663},{0.111540408937694457}},
    {1,-0.28341838511138821,3,{0.0337718440544803428,0.270347889165403623,0.476665439382152374},{0.0170909126716008895,0.218829882330447478,0.0699069619326525813},
     1,{0.0514643354866615260},{0.202706173746087101},{0.0609891857178808940}}
};

struct QUADRITURE_RULE_SEG
{
    bool use_center;
    double center_weight;
    int number;
    double location[3];
    double weight[3];
};

const QUADRITURE_RULE_SEG quadriture_rule_seg[12]={
    {1,1,0},
    {1,1,0},
    {0,0,1,{0.211324865405187118},{0.5}},
    {0,0,1,{0.211324865405187118},{0.5}},
    {1,0.444444444444444444,1,{0.112701665379258311},{0.277777777777777778}},
    {1,0.444444444444444444,1,{0.112701665379258311},{0.277777777777777778}},
    {0,0,2,{0.0694318442029737124,0.330009478207571866},{0.173927422568726928,0.326072577431273072}},
    {0,0,2,{0.0694318442029737124,0.330009478207571866},{0.173927422568726928,0.326072577431273072}},
    {1,0.284444444444444444,3,{0.0469100770306680036,0.230765344947158454},{0.118463442528094544,0.239314335249683234}},
    {1,0.284444444444444444,3,{0.0469100770306680036,0.230765344947158454},{0.118463442528094544,0.239314335249683234}},
    {0,0,3,{0.0337652428984239861,0.169395306766867739,0.380690406958401558},{0.0856622461895851725,0.180380786524069303,0.233956967286345524}},
    {0,0,3,{0.0337652428984239861,0.169395306766867739,0.380690406958401558},{0.0856622461895851725,0.180380786524069303,0.233956967286345524}}
};

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
    STATIC_POLYNOMIAL<T,rank,(d>d2?d:d2)> operator+ (const STATIC_POLYNOMIAL<T,rank,d2>& p) const
    {STATIC_POLYNOMIAL<T,rank,(d>d2?d:d2)> r;r=*this;r+=p;return r;}

    template<int d2>
    STATIC_POLYNOMIAL<T,rank,(d>d2?d:d2)> operator- (const STATIC_POLYNOMIAL<T,rank,d2>& p) const
    {STATIC_POLYNOMIAL<T,rank,(d>d2?d:d2)> r;r=*this;r-=p;return r;}

    template<int d2>
    STATIC_POLYNOMIAL<T,rank,(d+d2)> operator* (const STATIC_POLYNOMIAL<T,rank,d2>& p) const
    {STATIC_POLYNOMIAL<T,rank,(d+d2)> r;r.Multiply(*this,p,false);return r;}

    void Compress_Size()
    {
        RANGE<TV_INT> range(TV_INT(),size+1);
        size.Fill(0);
        for(RANGE_ITERATOR<rank> it(range);it.Valid();it.Next())
            if(terms(it.index))
                size=size.Componentwise_Max(it.index);
    }

    template<int d2,int d3> // Result must fit; no aliasing
    void Multiply(const STATIC_POLYNOMIAL<T,rank,d2>& a, const STATIC_POLYNOMIAL<T,rank,d3>& b, bool need_clear=true)
    {
        if(Same_Polynomial(a,*this)){
            STATIC_POLYNOMIAL<T,rank,d2> copy(a);
            if(Same_Polynomial(b,*this)) return Multiply(copy,copy);
            return Multiply(copy,b);}
        if(Same_Polynomial(b,*this)){
            STATIC_POLYNOMIAL<T,rank,d3> copy(b);
            return Multiply(a,copy);}

        size=a.size+b.size;
        PHYSBAM_ASSERT(size.All_Less_Equal(TV_INT()+d));
        RANGE<TV_INT> range(TV_INT(),a.size+1),range2(TV_INT(),b.size+1);
        if(need_clear){
            RANGE<TV_INT> range(TV_INT(),size+1);
            for(RANGE_ITERATOR<rank> it(range);it.Valid();it.Next()) terms(it.index)=0;}
        for(RANGE_ITERATOR<rank> it(range);it.Valid();it.Next()){
            T x=a.terms(it.index);
            if(!x) continue;
            for(RANGE_ITERATOR<rank> it2(range2);it2.Valid();it2.Next())
                terms(it.index+it2.index)+=x*b.terms(it2.index);}
    }

    template<int d2> // Result must fit
    STATIC_POLYNOMIAL& operator+=(const STATIC_POLYNOMIAL<T,rank,d2>& a)
    {
        PHYSBAM_ASSERT(a.size.All_Less_Equal(TV_INT()+d));
        RANGE<TV_INT> range(TV_INT(),a.size+1);
        for(RANGE_ITERATOR<rank> it(range);it.Valid();it.Next())
            terms(it.index)+=a.terms(it.index);
        size=size.Componentwise_Max(a.size);
        return *this;
    }

    template<int d2> // Result must fit
    STATIC_POLYNOMIAL& operator=(const STATIC_POLYNOMIAL<T,rank,d2>& a)
    {
        PHYSBAM_ASSERT(a.size.All_Less_Equal(TV_INT()+d));
        RANGE<TV_INT> range(TV_INT(),a.size+1);
        for(RANGE_ITERATOR<rank> it(range);it.Valid();it.Next())
            terms(it.index)=a.terms(it.index);
        size=a.size;
        return *this;
    }

    template<int d2> // Result must fit
    STATIC_POLYNOMIAL& operator-=(const STATIC_POLYNOMIAL<T,rank,d2>& a)
    {
        PHYSBAM_ASSERT(a.size.All_Less_Equal(TV_INT()+d));
        RANGE<TV_INT> range(TV_INT(),a.size+1);
        for(RANGE_ITERATOR<rank> it(range);it.Valid();it.Next())
            terms(it.index)-=a.terms(it.index);
        size=size.Componentwise_Max(a.size);
        return *this;
    }

    STATIC_POLYNOMIAL Differentiate(int v) const
    {
        STATIC_POLYNOMIAL r;
        RANGE<TV_INT> range(TV_INT::Axis_Vector(v),size+1);
        for(RANGE_ITERATOR<rank> it(range);it.Valid();it.Next())
            r.terms(it.index-TV_INT::Axis_Vector(v))=terms(it.index)*it.index(v);
        r.size=size;
        r.size(v)--;
        if(r.size(v)<0) r.size(v)=0;
        return r;
    }

    STATIC_POLYNOMIAL<T,rank,d+1> Integrate(int v) const
    {
        STATIC_POLYNOMIAL<T,rank,d+1> r;
        RANGE<TV_INT> range(TV_INT(),size+1);
        for(RANGE_ITERATOR<rank> it(range);it.Valid();it.Next())
            r.terms(it.index+TV_INT::Axis_Vector(v))=terms(it.index)/(it.index(v)+1);
        r.size=size;
        r.size(v)++;
        return r;
    }

    void Negate_Variable(int v)
    {
        RANGE<TV_INT> range(TV_INT(),size+1);
        for(RANGE_ITERATOR<rank> it(range);it.Valid();it.Next())
            if(it.index(v)&1)
                terms(it.index)=-terms(it.index);
    }

    void Scale(const TV& scale)
    {for(int i=0;i<rank;i++) Scale_Variable(i,scale(i));}

    void Scale_Variable(int v,T scale)
    {
        T powers[d+1];
        powers[0]=1;
        for(int j=1;j<=size(v);j++)
            powers[j]=powers[j-1]*scale;
        RANGE<TV_INT> range(TV_INT(),size+1);
        for(RANGE_ITERATOR<rank> it(range);it.Valid();it.Next())
            terms(it.index)*=powers[it.index(v)];
    }

    void Shift(const TV& scale)
    {for(int i=0;i<rank;i++) Shift_Variable(i,scale(i));}

    void Shift_Variable(int v,T shift)
    {
        T table[d+2][d+2]; // only need d+1, but use d+2 to work around 4.8 compiler warning
        table[0][0]=1;
        for(int i=1;i<=size(v) && i<d+1;i++){
            table[i][0]=1;
            table[i][i]=table[i-1][i-1]*shift;
            for(int j=1;j<i && j<d+1;j++)
                table[i][j]=table[i-1][j-1]*shift+table[i-1][j];}
        RANGE<TV_INT> range(TV_INT(),size+1);
        for(RANGE_ITERATOR<rank> it(range);it.Valid();it.Next()){
            T x=terms(it.index);
            TV_INT index=it.index;
            for(int i=1;i<=it.index(v);i++){
                index(v)--;
                terms(index)+=x*table[it.index(v)][i];}}
    }

    void Exchange_Variables(int u,int v)
    {
        if(u==v) return;
        if(size(u)>size(v)) exchange(u,v);
        RANGE<TV_INT> range(TV_INT(),size+1);
        range.max_corner(u)=1;
        for(RANGE_ITERATOR<rank> it(range);it.Valid();it.Next())
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
        for(RANGE_ITERATOR<rank> it(range);it.Valid();it.Next()){
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
        for(RANGE_ITERATOR<rank> it(range);it.Valid();it.Next()){
            T t=terms(it.index);
            for(int j=0;j<rank;j++)
                t*=powers[j][it.index(j)];
            s+=t;}
        return s;
    }
    
    T Integrate_Over_Primitive(const VECTOR<TV,2>& vertices) const
    {
        STATIC_POLYNOMIAL copy(*this);
        
        copy.Shift(vertices(0));
        TV a=vertices(1)-vertices(0);
        
        T table[rank][d+1];
        for(int v=0;v<rank;v++){
            table[v][0]=1;
            for(int i=1;i<=copy.size(v);i++) table[v][i]=table[v][i-1]*a(v);}
        
        typedef VECTOR<int,1> TV_INT1;
        
        STATIC_POLYNOMIAL<T,1,d*rank> barycentric;
        RANGE<TV_INT> range(TV_INT(),copy.size+1);
        for(RANGE_ITERATOR<rank> it(range);it.Valid();it.Next()){
            T scale=1;
            for(int v=0;v<rank;v++) scale*=table[v][it.index(v)];
            barycentric.terms(TV_INT1(it.index.Sum()))+=copy.terms(it.index)*scale;}
        barycentric.size=TV_INT1(size.Sum());

        T integral=0;
        for(int i=0;i<=barycentric.size(0);i++)
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
    
        typedef VECTOR<int,2> TV_INT2;

        STATIC_POLYNOMIAL<T,2,d*rank> barycentric;
        RANGE<TV_INT> range(TV_INT(),size+1);
        for(RANGE_ITERATOR<rank> it(range);it.Valid();it.Next()){
            STATIC_POLYNOMIAL<T,2,rank*d> monomial;
            monomial.Set_Term(TV_INT2(),copy.terms(it.index));
            for(int v=0;v<rank;v++){
                STATIC_POLYNOMIAL<T,2,d> variable;
                int n=it.index(v);
                for(int j=0;j<=n;j++)
                    variable.terms(TV_INT2(j,n-j))=table[v][n][j];
                variable.size=TV_INT2(n,n);
                monomial.Multiply(monomial,variable);}
            barycentric+=monomial;}
    
        int max_factorial=barycentric.size.Sum()+2;
        T factorial[2*rank*d+2];
        factorial[0]=(T)1;
        for(int i=1;i<=max_factorial;i++) factorial[i]=factorial[i-1]*i;        

        T integral=0;
        RANGE<TV_INT2> range2(TV_INT2(),barycentric.size+1);
        for(RANGE_ITERATOR<2> it(range2);it.Valid();it.Next()){
            integral+=barycentric.terms(it.index)*factorial[it.index(0)]*factorial[it.index(1)]/factorial[it.index(0)+it.index(1)+2];}

        return integral*TV::Cross_Product(a,b).Magnitude();
    }

    static T Eval(const TV& pt,const TV_INT& p)
    {
        T r=1;
        for(int i=0;i<TV::m;i++)
            for(int j=0;j<p(i);j++)
                r*=pt(i);
        return r;
    }

    static T Quadrature_Over_Primitive(const VECTOR<TV,2>& vertices,const TV_INT& p)
    {
        int o=p.Sum();
        PHYSBAM_ASSERT((unsigned)o<sizeof(quadriture_rule_seg)/sizeof(*quadriture_rule_seg));
        TV a=vertices(1)-vertices(0);
        if(!o) return a.Magnitude();
        T r=0;
        
        const QUADRITURE_RULE_SEG& rule=quadriture_rule_seg[o];
        if(rule.use_center) r+=rule.center_weight*Eval(vertices.Sum()/2,p);
        for(int i=0;i<rule.number;i++){
            T s=0;
            s+=Eval(vertices.x+(T)rule.location[i]*a,p);
            s+=Eval(vertices.y-(T)rule.location[i]*a,p); 
            r+=rule.weight[i]*s;}
        return r*a.Magnitude();
    }

    T Quadrature_Over_Primitive(const VECTOR<TV,2>& vertices) const
    {
        int o=size.Sum();
        PHYSBAM_ASSERT((unsigned)o<sizeof(quadriture_rule_seg)/sizeof(*quadriture_rule_seg));
        TV a=vertices(1)-vertices(0);
        if(!o) return terms(TV_INT())*a.Magnitude();
        T r=0;
        
        const QUADRITURE_RULE_SEG& rule=quadriture_rule_seg[o];
        if(rule.use_center) r+=rule.center_weight*Value(vertices.Sum()/2);
        for(int i=0;i<rule.number;i++){
            T s=0;
            s+=Value(vertices.x+(T)rule.location[i]*a);
            s+=Value(vertices.y-(T)rule.location[i]*a); 
            r+=rule.weight[i]*s;}
        return r*a.Magnitude();
    }

    static T Quadrature_Over_Primitive(const VECTOR<TV,3>& vertices,const TV_INT& p)
    {
        int o=p.Sum();
        PHYSBAM_ASSERT((unsigned)o<sizeof(quadriture_rule_tri)/sizeof(*quadriture_rule_tri));
        TV a=vertices(1)-vertices(0);
        TV b=vertices(2)-vertices(0);
        if(!o) return TV::Cross_Product(a,b).Magnitude()/2;
        T r=0;
        
        const QUADRITURE_RULE_TRI& rule=quadriture_rule_tri[o];
        if(rule.use_center) r+=rule.center_weight*Eval(vertices.Sum()/3,p);
        for(int i=0;i<rule.number3;i++){
            T s=0;
            s+=Eval(vertices.x+(T)rule.location3[i]*(vertices.y-vertices.x)+(T)rule.location3[i]*(vertices.z-vertices.x),p);
            s+=Eval(vertices.y+(T)rule.location3[i]*(vertices.z-vertices.y)+(T)rule.location3[i]*(vertices.x-vertices.y),p);
            s+=Eval(vertices.z+(T)rule.location3[i]*(vertices.x-vertices.z)+(T)rule.location3[i]*(vertices.y-vertices.z),p);
            r+=rule.weight3[i]*s;}
        for(int i=0;i<rule.number6;i++){
            T s=0;
            s+=Eval(vertices.x+(T)rule.location6a[i]*(vertices.y-vertices.x)+(T)rule.location6b[i]*(vertices.z-vertices.x),p);
            s+=Eval(vertices.y+(T)rule.location6a[i]*(vertices.z-vertices.y)+(T)rule.location6b[i]*(vertices.x-vertices.y),p);
            s+=Eval(vertices.z+(T)rule.location6a[i]*(vertices.x-vertices.z)+(T)rule.location6b[i]*(vertices.y-vertices.z),p);
            s+=Eval(vertices.x+(T)rule.location6b[i]*(vertices.y-vertices.x)+(T)rule.location6a[i]*(vertices.z-vertices.x),p);
            s+=Eval(vertices.y+(T)rule.location6b[i]*(vertices.z-vertices.y)+(T)rule.location6a[i]*(vertices.x-vertices.y),p);
            s+=Eval(vertices.z+(T)rule.location6b[i]*(vertices.x-vertices.z)+(T)rule.location6a[i]*(vertices.y-vertices.z),p);
            r+=rule.weight6[i]*s;}
        return r*TV::Cross_Product(a,b).Magnitude()/2;
    }

    T Quadrature_Over_Primitive(const VECTOR<TV,3>& vertices) const
    {
        int o=size.Sum();
        PHYSBAM_ASSERT((unsigned)o<sizeof(quadriture_rule_seg)/sizeof(*quadriture_rule_seg));
        TV a=vertices(1)-vertices(0);
        TV b=vertices(2)-vertices(0);
        if(!o) return terms(TV_INT())*TV::Cross_Product(a,b).Magnitude()/2;
        T r=0;
        const QUADRITURE_RULE_TRI& rule=quadriture_rule_tri[o];
        if(rule.use_center) r+=rule.center_weight*Value(vertices.Sum()/3);
        for(int i=0;i<rule.number3;i++){
            T s=0;
            s+=Value(vertices.x+(T)rule.location3[i]*(vertices.y-vertices.x)+(T)rule.location3[i]*(vertices.z-vertices.x));
            s+=Value(vertices.y+(T)rule.location3[i]*(vertices.z-vertices.y)+(T)rule.location3[i]*(vertices.x-vertices.y));
            s+=Value(vertices.z+(T)rule.location3[i]*(vertices.x-vertices.z)+(T)rule.location3[i]*(vertices.y-vertices.z));
            r+=rule.weight3[i]*s;}
        for(int i=0;i<rule.number6;i++){
            T s=0;
            s+=Value(vertices.x+(T)rule.location6a[i]*(vertices.y-vertices.x)+(T)rule.location6b[i]*(vertices.z-vertices.x));
            s+=Value(vertices.y+(T)rule.location6a[i]*(vertices.z-vertices.y)+(T)rule.location6b[i]*(vertices.x-vertices.y));
            s+=Value(vertices.z+(T)rule.location6a[i]*(vertices.x-vertices.z)+(T)rule.location6b[i]*(vertices.y-vertices.z));
            s+=Value(vertices.x+(T)rule.location6b[i]*(vertices.y-vertices.x)+(T)rule.location6a[i]*(vertices.z-vertices.x));
            s+=Value(vertices.y+(T)rule.location6b[i]*(vertices.z-vertices.y)+(T)rule.location6a[i]*(vertices.x-vertices.y));
            s+=Value(vertices.z+(T)rule.location6b[i]*(vertices.x-vertices.z)+(T)rule.location6a[i]*(vertices.y-vertices.z));
            r+=rule.weight6[i]*s;}
        return r*TV::Cross_Product(a,b).Magnitude()/2;
    }
};

template<class T,int rank,int d>
std::ostream& operator<< (std::ostream& o, const STATIC_POLYNOMIAL<T,rank,d>& p)
{
    bool first=true;
    o<<"(";
    RANGE<VECTOR<int,rank> > range(VECTOR<int,rank>(),p.size+1);
    for(RANGE_ITERATOR<rank> it(range);it.Valid();it.Next()){
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
