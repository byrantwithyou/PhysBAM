//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTIAL_UNIT
//#####################################################################
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Log/LOG.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include "PARTIAL_UNIT.h"
#include <boost/format.hpp>
namespace PhysBAM{
namespace UNITS{
//#####################################################################
// Stream output
//#####################################################################
struct MONOMIAL
{
    std::string x;
    rational<int> p;

    MONOMIAL()
    {}

    MONOMIAL(const std::string& x,const rational<int>& p)
        :x(x),p(p)
    {}
};

inline std::ostream& operator<<(std::ostream& output,const MONOMIAL& m)
{
    output<<m.x; 
    if(m.p!=1){
        output<<'^'<<m.p.numerator();
        if(m.p.denominator()!=1) output<<'/'<<m.p.denominator();}
    return output;
}

std::ostream& operator<<(std::ostream& output,const SIMPLE_UNIT& unit)
{ARRAY<MONOMIAL> monomials;
if(unit.m!=0) monomials.Append(MONOMIAL("m",unit.m));
if(unit.kg!=0) monomials.Append(MONOMIAL("kg",unit.kg));
if(unit.s!=0) monomials.Append(MONOMIAL("s",unit.s));
if(!monomials.m) output<<"()";
else{output<<monomials(1);for(int i=2;i<=monomials.m;i++) output<<'.'<<monomials(i);}
return output;}

std::ostream& operator<<(std::ostream& output,const PARTIAL_UNIT& unit)
{ARRAY<MONOMIAL> monomials;
if(unit.a.m!=0) monomials.Append(MONOMIAL("m",unit.a.m));
if(unit.a.kg!=0) monomials.Append(MONOMIAL("kg",unit.a.kg));
if(unit.a.s!=0) monomials.Append(MONOMIAL("s",unit.a.s));
if(unit.power!=0) monomials.Append(MONOMIAL(STRING_UTILITIES::string_sprintf("x%d",unit.x),unit.power));
if(!monomials.m) output<<"()";
else{output<<monomials(1);for(int i=2;i<=monomials.m;i++) output<<'.'<<monomials(i);}
return output;}

//#####################################################################
// Class UNIT_ENVIRONMENT
//#####################################################################
namespace{
class UNIT_ENVIRONMENT:public NONCOPYABLE
{
private:
    struct ENTRY:PARTIAL_UNIT 
    {
        ENTRY()
            :PARTIAL_UNIT(PARTIAL_UNIT::One())
        {}

        ENTRY(const SIMPLE_UNIT& u)
            :PARTIAL_UNIT(u)
        {}

        ENTRY(const PARTIAL_UNIT& u)
            :PARTIAL_UNIT(u)
        {}
    };
    mutable HASHTABLE<VARIABLE,ENTRY> env;
    mutable VARIABLE next_variable;
public:

    UNIT_ENVIRONMENT()
        :next_variable(1)
    {}

    static UNIT_ENVIRONMENT& Singleton()
    {static UNIT_ENVIRONMENT environment;return environment;}

    PARTIAL_UNIT New() const
    {if(next_variable==0) PHYSBAM_FATAL_ERROR("Ran out of unit variables");
    return PARTIAL_UNIT(SIMPLE_UNIT::One(),next_variable++,1);}

    void Find(PARTIAL_UNIT& u) const
    {//LOG::cout<<"FIND "<<u<<std::endl;
    if(u.power==0 || !env.Contains(u.x)) return;
    PARTIAL_UNIT& ux=env.Get(u.x);Find(ux);
    u=u.a*pow(ux,u.power);}

    void Unify(PARTIAL_UNIT u1,PARTIAL_UNIT u2) const
    {if(u1.a.Is_Degenerate() || u2.a.Is_Degenerate()) return;
    Find(u1);Find(u2);
    if(u1.power==0 && u2.power==0){
        if(u1.a!=u2.a) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Can't unify %s with %s",u1,u2));}
    else if(u1.x==u2.x){
        if(u1.power==u2.power){
            if(u1.a!=u2.a) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Can't unify %s with %s",u1,u2));}
        else env.Insert(u1.x,pow(u2.a/u1.a,1/(u1.power-u2.power)));}
    else if(u2.power==0 || (u1.power!=0 && abs(u1.power)<abs(u2.power)))
        env.Insert(u1.x,pow(u2/u1.a,1/u1.power));
    else env.Insert(u2.x,pow(u1/u2.a,1/u2.power));}

    void Unify_One(PARTIAL_UNIT u) const
    {if(u.a.Is_Degenerate()) return;
    Find(u);
    if(u.power==0){
        if(u.a!=SIMPLE_UNIT::One()) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Can't unify %s with 1",u));}
    else env.Insert(u.x,pow(u.a.Inverse(),1/u.power));}
};
}
//#####################################################################
// Function New
//#####################################################################
PARTIAL_UNIT PARTIAL_UNIT::New()
{
    return UNIT_ENVIRONMENT::Singleton().New();
}
//#####################################################################
// Function Unify
//#####################################################################
void Unify(const PARTIAL_UNIT& u1,const PARTIAL_UNIT& u2)
{
    UNIT_ENVIRONMENT::Singleton().Unify(u1,u2);
}
//#####################################################################
// Function Unify_One
//#####################################################################
void Unify_One(const PARTIAL_UNIT& u)
{
    UNIT_ENVIRONMENT::Singleton().Unify_One(u);
}
//#####################################################################
// Function Find
//#####################################################################
void PARTIAL_UNIT::Find() const
{
    UNIT_ENVIRONMENT::Singleton().Find(const_cast<PARTIAL_UNIT&>(*this));
}
//#####################################################################
}
}
