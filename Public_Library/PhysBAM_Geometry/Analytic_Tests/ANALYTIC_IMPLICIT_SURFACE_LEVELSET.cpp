//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_IMPLICIT_SURFACE_LEVELSET.h>
namespace PhysBAM{
//#####################################################################
// Function Find_Closest_Point
//#####################################################################
template<class TV> VECTOR<typename TV::SCALAR,TV::m+1> ANALYTIC_IMPLICIT_SURFACE_LEVELSET<TV>::
Find_Closest_Point(const TV& X) const
{
    TV w=Closest_Point_Estimate(X);
    VECTOR<T,TV::m+1> z(w.Append(0));
    for(int i=0;i<100;i++){
        TV Z(z.Remove_Index(TV::m)),dg=df(Z);
        T L=z(TV::m),g=f(Z);
        VECTOR<T,TV::m+1> G=((Z-X)*2+L*dg).Append(g),Hcol=dg.Append(0);
        MATRIX<T,TV::m+1> H;
        H.Set_Submatrix(0,0,L*ddf(Z)+2);
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
Phi(const TV& X) const
{
    return (X-Find_Closest_Point(X).Remove_Index(TV::m)).Magnitude()*sign(f(X));
}
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV ANALYTIC_IMPLICIT_SURFACE_LEVELSET<TV>::
Normal(const TV& X) const
{
    return (X-Find_Closest_Point(X).Remove_Index(TV::m)).Normalized();
}
template class ANALYTIC_IMPLICIT_SURFACE_LEVELSET<VECTOR<float,1> >;
template class ANALYTIC_IMPLICIT_SURFACE_LEVELSET<VECTOR<float,2> >;
template class ANALYTIC_IMPLICIT_SURFACE_LEVELSET<VECTOR<float,3> >;
template class ANALYTIC_IMPLICIT_SURFACE_LEVELSET<VECTOR<double,1> >;
template class ANALYTIC_IMPLICIT_SURFACE_LEVELSET<VECTOR<double,2> >;
template class ANALYTIC_IMPLICIT_SURFACE_LEVELSET<VECTOR<double,3> >;
}
