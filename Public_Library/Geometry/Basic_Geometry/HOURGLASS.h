//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HOURGLASS
//##################################################################### 
#ifndef __HOURGLASS__
#define __HOURGLASS__

#include <Tools/Math_Tools/constants.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/AUTODIFF_LEVELSET.h>
namespace PhysBAM{

template<class TV> class HOURGLASS;

template<class T>
class HOURGLASS<VECTOR<T,2> >:public AUTODIFF_LEVELSET<VECTOR<T,2>,HOURGLASS<VECTOR<T,2> > >
{
    typedef VECTOR<T,2> TV;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef TV VECTOR_T;

    MATRIX<T,TV::m> R;
    TV axis,center;
    T o,q,r0,r1;
    TV N,C0,C1,v0;

    HOURGLASS(const TV& axis,const TV& center,T bulb_radius,T neck_radius,T height,T neck_width);

    template<class IN>
    auto Raw_Phi_Helper(const IN& Z) const
    {
        auto p=Z.x.Dot(v0);
        decltype(min(min(r1-(Z-C1).Magnitude(),(Z-C0).Magnitude()-r0),r1-(Z-C1).Dot(N))) ret;
        if(p<o) ret.Fill_From(r1-(Z-C1).Magnitude());
        else if(p>q) ret.Fill_From((Z-C0).Magnitude()-r0);
        else ret.Fill_From(r1-(Z-C1).Dot(N));
        return ret;
    }

    template<class IN>
    auto Raw_Phi(const IN& X) const
    {
        auto Z=R*(X.x-center);
        DIAGONAL_MATRIX<T,2> dm(sign_nonzero(Z.x),sign_nonzero(Z.y));
        return Raw_Phi_Helper((dm*R)*(X-center));
    }

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,R,center,o,q,N,C0,C1,v0,r0,r1);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,R,center,o,q,N,C0,C1,v0,r0,r1);}

//#####################################################################
    RANGE<TV> Bounding_Box() const;
    static std::string Name() {return "HOURGLASS<VECTOR<T,2> >";}
//#####################################################################
};   

template<class T>
class HOURGLASS<VECTOR<T,3> >:public AUTODIFF_LEVELSET<VECTOR<T,3>,HOURGLASS<VECTOR<T,3> > >
{
    typedef VECTOR<T,3> TV;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef TV VECTOR_T;

    TV axis,center;
    MATRIX<T,2,3> M;
    HOURGLASS<VECTOR<T,2> > hg2;
    RANGE<TV> bounding_box;

    HOURGLASS(const TV& axis,const TV& center,T bulb_radius,T neck_radius,T height,T neck_width);

    template<class IN>
    auto Raw_Phi(const IN& X) const
    {
        auto Z=X-center;
        return hg2.Raw_Phi_Helper(VECTOR<T,2>(0,1)*abs(Z.Dot(axis))+VECTOR<T,2>(1,0)*(M*Z).Magnitude());
    }

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,axis,center,M,hg2);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,axis,center,M,hg2);}

//#####################################################################
    RANGE<TV> Bounding_Box() const {return bounding_box;}
    static std::string Name() {return "HOURGLASS<VECTOR<T,3> >";}
//#####################################################################
};   
template<class T> class TRIANGULATED_SURFACE;
template<class T> class SEGMENTED_CURVE_2D;
namespace TESSELLATION{
template<class T> SEGMENTED_CURVE_2D<T>* Tessellate_Boundary(const HOURGLASS<VECTOR<T,2> >& hourglass,int axis_div=32);
template<class T> TRIANGULATED_SURFACE<T>* Tessellate_Boundary(const HOURGLASS<VECTOR<T,3> >& hourglass,int axis_div=32,int circ_div=16);
}
}
#endif
