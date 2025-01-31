//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Geometry/Analytic_Tests/ANALYTIC_IMPLICIT_SURFACE_LEVELSET.h>
namespace PhysBAM{
//#####################################################################
// Function Find_Closest_Point
//#####################################################################
template<class TV> VECTOR<typename TV::SCALAR,TV::m+1> ANALYTIC_IMPLICIT_SURFACE_LEVELSET<TV>::
Find_Closest_Point(const TV& X,T t) const
{
    TV w=Closest_Point_Estimate(X,t);
    VECTOR<T,TV::m+1> z(w.Append(0));
    for(int i=0;i<100;i++){
        TV Z(z.Remove_Index(TV::m)),dg=df(Z,t);
        T L=z(TV::m),g=f(Z,t);
        VECTOR<T,TV::m+1> G=((Z-X)*2+L*dg).Append(g),Hcol=dg.Append(0);
        MATRIX<T,TV::m+1> H;
        H.Set_Submatrix(0,0,L*ddf(Z,t)+2);
        H.Set_Row(TV::m,Hcol);
        H.Set_Column(TV::m,Hcol);
        z-=H.In_Place_PLU_Solve(G);
        if(G.Magnitude_Squared()<1e-25) break;}
    return z;
}
//#####################################################################
// Function Phi
//#####################################################################
template<class TV> typename TV::SCALAR ANALYTIC_IMPLICIT_SURFACE_LEVELSET<TV>::
phi2(const TV& X,T t) const
{
    return (X-Find_Closest_Point(X,t).Remove_Index(TV::m)).Magnitude()*sign(f(X,t));
}
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV ANALYTIC_IMPLICIT_SURFACE_LEVELSET<TV>::
N2(const TV& X,T t) const
{
    return (X-Find_Closest_Point(X,t).Remove_Index(TV::m)).Normalized()*sign(f(X,t));
}
template class ANALYTIC_IMPLICIT_SURFACE_LEVELSET<VECTOR<float,1> >;
template class ANALYTIC_IMPLICIT_SURFACE_LEVELSET<VECTOR<float,2> >;
template class ANALYTIC_IMPLICIT_SURFACE_LEVELSET<VECTOR<float,3> >;
template class ANALYTIC_IMPLICIT_SURFACE_LEVELSET<VECTOR<double,1> >;
template class ANALYTIC_IMPLICIT_SURFACE_LEVELSET<VECTOR<double,2> >;
template class ANALYTIC_IMPLICIT_SURFACE_LEVELSET<VECTOR<double,3> >;
}
