//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MULTIVARIATE_POLYNOMIAL
//#####################################################################
#ifndef __MULTIVARIATE_POLYNOMIAL__
#define __MULTIVARIATE_POLYNOMIAL__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV>
struct MULTIVARIATE_MONOMIAL
{
    typedef VECTOR<int,TV::m> TV_INT;
    typedef typename TV::ELEMENT T;

    TV_INT power;
    T coeff;

    MULTIVARIATE_MONOMIAL(): coeff(0) {}
    MULTIVARIATE_MONOMIAL(const TV_INT& p, T c): power(p), coeff(c) {}

    bool operator< (const MULTIVARIATE_MONOMIAL& m) const
    {
        return LEXICOGRAPHIC_COMPARE()(m.power,power); // Leading monomial first
    }

    void Merge(MULTIVARIATE_MONOMIAL& m)
    {
        assert(power==m.power);
        coeff+=m.coeff;
    }
};

template<class TV>
struct MULTIVARIATE_POLYNOMIAL
{
    typedef VECTOR<int,TV::m> TV_INT;
    typedef typename TV::ELEMENT T;

    ARRAY<MULTIVARIATE_MONOMIAL<TV> > terms;

    MULTIVARIATE_POLYNOMIAL();
    ~MULTIVARIATE_POLYNOMIAL();

    void Simplify();
    TV_INT Max_Power() const;
    MULTIVARIATE_POLYNOMIAL& operator+= (const MULTIVARIATE_POLYNOMIAL& m);
    MULTIVARIATE_POLYNOMIAL& operator-= (const MULTIVARIATE_POLYNOMIAL& m);
    MULTIVARIATE_POLYNOMIAL& operator*= (const MULTIVARIATE_POLYNOMIAL& m);
    void Differentiate(int v);
    void Integrate(int v);
    void Negate_Variable(int v);
    void Scale_Variables(const TV& scale);
    void Shift_Variable(int v,T shift);
    void Shift(const TV& shift);
    void Exchange_Variables(int u,int v);
    T Definite_Integral(RANGE<TV>& domain) const;
    T Integrate_Over_Triangle(const VECTOR<TV,3>& vertices) const;
};

template<class TV>
std::ostream& operator<< (std::ostream& o, const MULTIVARIATE_POLYNOMIAL<TV>& p);
}
#endif
