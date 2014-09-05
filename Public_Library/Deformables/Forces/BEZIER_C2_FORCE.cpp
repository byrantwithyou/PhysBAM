//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Deformables/Forces/BEZIER_C2_FORCE.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BEZIER_C2_FORCE<TV>::
BEZIER_C2_FORCE(DEFORMABLE_PARTICLES<TV>& particles,const BEZIER_SPLINE<TV,3>& spline_input,
    T stiffness_input)
    :DEFORMABLES_FORCES<TV>(particles),spline(spline_input),stiffness(stiffness_input)
{
    X0=spline.particles.X;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BEZIER_C2_FORCE<TV>::
~BEZIER_C2_FORCE()
{
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void BEZIER_C2_FORCE<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(int i=0;i<data.m;i++)
        F.Subset(data(i).pts)-=data(i).ge;

    for(int i=0;i<2;i++)
        F.Subset(end_pts[i])-=end_ge[i];
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void BEZIER_C2_FORCE<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void BEZIER_C2_FORCE<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T scale,const T time) const
{
    for(int i=0;i<data.m;i++){
        const MATRIX<MATRIX<T,TV::m>,4>& he=data(i).he;
        VECTOR<TV,4> v(V.Subset(data(i).pts)*scale),f;
        for(int j=0;j<4;j++)
            for(int k=0;k<4;k++)
                f(j)+=he(j,k)*v(k);
        F.Subset(data(i).pts)-=f;}

    for(int i=0;i<2;i++){
        const MATRIX<MATRIX<T,TV::m>,3>& he=end_he[i];
        VECTOR<TV,3> v(V.Subset(end_pts[i])*scale),f;
        for(int j=0;j<3;j++)
            for(int k=0;k<3;k++)
                f(j)+=he(j,k)*v(k);
        F.Subset(end_pts[i])-=f;}
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void BEZIER_C2_FORCE<TV>::
Initialize_CFL(ARRAY_VIEW<typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA> frequency)
{
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR BEZIER_C2_FORCE<TV>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void BEZIER_C2_FORCE<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void BEZIER_C2_FORCE<TV>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    MATRIX<T,2> rot(0,1,-1,0);
    pe=0;
    data.Resize(spline.control_points.m-1);
    for(int i=0;i<spline.control_points.m-1;i++){
        const VECTOR<int,4>& nodes0=spline.control_points(i);
        const VECTOR<int,4>& nodes1=spline.control_points(i+1);
        DATA& dat=data(i);
        dat.pts=VECTOR<int,4>(nodes0(1),nodes0(2),nodes1(1),nodes1(2));
        auto Z0=From_Var<4,0>(particles.X(nodes0(1)));
        auto Z1=From_Var<4,1>(particles.X(nodes0(2)));
        auto Z2=From_Var<4,2>(particles.X(nodes1(1)));
        auto Z3=From_Var<4,3>(particles.X(nodes1(2)));
        auto dev=Z0-Z1*2+Z2*2-Z3;
        auto new_pe=stiffness/2*dev.Magnitude_Squared();
        pe+=new_pe.x;
        Extract(dat.ge,new_pe.dx);
        Extract(dat.he,new_pe.ddx);}

    const VECTOR<int,4>& n0=spline.control_points(0);
    const VECTOR<int,4>& n1=spline.control_points.Last();
    end_pts[0]={n0(0),n0(1),n0(2)};
    end_pts[1]={n1(3),n1(2),n1(1)};
    for(int i=0;i<2;i++){
        auto Z0=From_Var<3,0>(particles.X(end_pts[i][0]));
        auto Z1=From_Var<3,1>(particles.X(end_pts[i][1]));
        auto Z2=From_Var<3,2>(particles.X(end_pts[i][2]));
        auto dev=Z0-Z1*2+Z2;
        auto new_pe=stiffness/2*dev.Magnitude_Squared();
        pe+=new_pe.x;
        Extract(end_ge[i],new_pe.dx);
        Extract(end_he[i],new_pe.ddx);}
    printf("PE: %g\n",pe);
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR BEZIER_C2_FORCE<TV>::
Potential_Energy(const T time) const
{
    return pe;
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void BEZIER_C2_FORCE<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
template class BEZIER_C2_FORCE<VECTOR<float,2> >;
template class BEZIER_C2_FORCE<VECTOR<double,2> >;
}
