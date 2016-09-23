//#####################################################################
// Copyright 2002-2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DYNAMIC_IMPLICIT_SURFACE
//#####################################################################
#ifndef __DYNAMIC_IMPLICIT_SURFACE__
#define __DYNAMIC_IMPLICIT_SURFACE__

#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <string>
namespace PhysBAM{

template<class T>
struct DYNAMIC_OP
{
    typedef VECTOR<T,3> TV;
    ARRAY<DYNAMIC_OP<T>*> A;

    T S;
    TV D;
    MATRIX<T,3> H;

    virtual ~DYNAMIC_OP() {}
    virtual DYNAMIC_OP<T>* New() const=0;
    virtual void Compute(const TV& X)=0;
};

template<class T_input>
class DYNAMIC_IMPLICIT_SURFACE:public IMPLICIT_OBJECT<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;
public:
    typedef IMPLICIT_OBJECT<TV> BASE;
    using BASE::box;

    HASHTABLE<std::string,DYNAMIC_OP<T>*> token_lookup;
    HASHTABLE<std::string,DYNAMIC_OP<T>*> dict;
    mutable ARRAY<DYNAMIC_OP<T>*> eval_list;

    T minimum_cell_size;

    DYNAMIC_IMPLICIT_SURFACE(const RANGE<TV>& range,const T minimum_cell_size=0)
        :minimum_cell_size(minimum_cell_size)
    {
        box=range;
        Initialize_Tokens();
    }

    void Update_Box() override
    {}

    virtual ~DYNAMIC_IMPLICIT_SURFACE();

    void Initialize_Tokens();
    std::string Next_Token(const char*& str) const;
    DYNAMIC_OP<T>* Parse(const char*& str);
    void Eval(const TV& X) const;
    void Set_Expression(const std::string& expr);

    T Minimum_Cell_Size() const override
    {return minimum_cell_size;}

    T Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const override
    {return minimum_cell_size;}

    T operator()(const TV& location) const override
    {return Extended_Phi(location);}

    T Extended_Phi(const TV& location) const override
    {Eval(location);//LOG::cout<<location<<" -> "<<eval_list.Last()->S<<std::endl;
return eval_list.Last()->S;}

    TV Normal(const TV& location,const int aggregate=-1) const override
    {Eval(location);return eval_list.Last()->D.Normalized();}

    TV Extended_Normal(const TV& location,const int aggregate=-1) const override
    {return Normal(location);}

    VECTOR<T,2> Principal_Curvatures(const TV& X) const override;

    bool Lazy_Inside(const TV& location,const T contour_value=0) const override
    {return Extended_Phi(location)<=contour_value;}

    bool Lazy_Inside_And_Value(const TV& location,T& phi_value,const T contour_value=0) const override
    {phi_value=Extended_Phi(location);return phi_value<=contour_value;}

    bool Lazy_Outside(const TV& location,const T contour_value=0) const override
    {return Extended_Phi(location)>=contour_value;}

    virtual void Read(TYPED_ISTREAM& input) override {PHYSBAM_FATAL_ERROR();}
    virtual void Write(TYPED_OSTREAM& output) const override {PHYSBAM_FATAL_ERROR();}

//#####################################################################
};
}
#endif
