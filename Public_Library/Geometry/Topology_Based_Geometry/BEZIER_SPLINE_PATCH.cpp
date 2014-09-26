//#####################################################################
// Copyright 2014
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BEZIER_SPLINE_PATCH
//##################################################################### 
#include <Tools/Arrays/IDENTITY_ARRAY.h>
#include <Tools/Matrices/BANDED_MATRIX.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Topology_Based_Geometry/BEZIER_SPLINE_PATCH.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> BEZIER_SPLINE_PATCH<TV,d>::
BEZIER_SPLINE_PATCH()
    :need_destroy_particles(true),particles(*new GEOMETRY_PARTICLES<TV>)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> BEZIER_SPLINE_PATCH<TV,d>::
BEZIER_SPLINE_PATCH(GEOMETRY_PARTICLES<TV>& particles)
    :need_destroy_particles(false),particles(particles)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> BEZIER_SPLINE_PATCH<TV,d>::
~BEZIER_SPLINE_PATCH()
{
    if(need_destroy_particles) delete &particles;
}
//#####################################################################
// Function Append_Particles_And_Create_Copy
//#####################################################################
template<class TV,int d> BEZIER_SPLINE_PATCH<TV,d>* BEZIER_SPLINE_PATCH<TV,d>::
Append_Particles_And_Create_Copy(GEOMETRY_PARTICLES<TV>& new_particles,ARRAY<int>* particle_indices) const
{
    BEZIER_SPLINE_PATCH* bs=new BEZIER_SPLINE_PATCH(new_particles);
    int offset=new_particles.Size();
    new_particles.Append(particles);
    if(particle_indices) particle_indices->Append_Elements(IDENTITY_ARRAY<>(particles.Size())+offset);
    bs->control_points=control_points+offset;
    return bs;
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV,int d> void BEZIER_SPLINE_PATCH<TV,d>::
Read(TYPED_ISTREAM& input)
{
    int size;
    Read_Binary(input,control_points,size);
    particles.Clean_Memory(); // strip everything away except for position
    particles.Resize(size);
    if(input.type.use_doubles) Read_Binary_Array<double>(input.stream,particles.X.Get_Array_Pointer(),size);
    else Read_Binary_Array<float>(input.stream,particles.X.Get_Array_Pointer(),size);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV,int d> void BEZIER_SPLINE_PATCH<TV,d>::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,control_points,particles.X);
}
//#####################################################################
// Function Smooth_Fit
//#####################################################################
template<class TV> void PhysBAM::
Smooth_Fit(BEZIER_SPLINE_PATCH<TV,3>& bs,ARRAY_VIEW<TV> X)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Smooth_Fit_Loop
//#####################################################################
template<class TV> void PhysBAM::
Smooth_Fit_Loop(BEZIER_SPLINE_PATCH<TV,3>& bs,ARRAY_VIEW<TV> X)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Create_Segmented_Curve
//#####################################################################
template<class TV,int d> typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,1>::OBJECT* PhysBAM::
Create_Segmented_Curve(const BEZIER_SPLINE_PATCH<TV,d>& spline,bool same_particles)
{
    /////TODO
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int d> BEZIER_SPLINE_PATCH<TV,d>* BEZIER_SPLINE_PATCH<TV,d>::
Create()
{
    return new BEZIER_SPLINE_PATCH;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int d> BEZIER_SPLINE_PATCH<TV,d>* BEZIER_SPLINE_PATCH<TV,d>::
Create(GEOMETRY_PARTICLES<TV>& particles)
{
    return new BEZIER_SPLINE_PATCH(particles);
}
namespace{
const int comb[9][9]=
{
    {1},
    {1,1},
    {1,2,1},
    {1,3,3,1},
    {1,4,6,4,1},
    {1,5,10,10,5,1},
    {1,6,15,20,15,6,1},
    {1,7,21,35,35,21,7,1},
    {1,8,28,56,70,56,28,8,1}
};
}
//#####################################################################
// Function Evaluate
//#####################################################################
template<class TV,int d> TV BEZIER_SPLINE_PATCH<TV,d>::
Evaluate(int id,T s,T t) const
{
    VECTOR<TV,(d+1)*(d+1)> X(particles.X.Subset(control_points(id)));
    T pow_t[d+1]={1};
    T pow_u[d+1]={1};
    T pow_s[d+1]={1};
    T pow_r[d+1]={1};
    T u=1-t; T r=1-s;
    for(int i=0;i<d;i++){
        pow_t[i+1]=pow_t[i]*t;
        pow_u[i+1]=pow_u[i]*u;
        pow_s[i+1]=pow_s[i]*s;
        pow_r[i+1]=pow_r[i]*r;}

    TV z;
    for(int i=0;i<d+1;i++)
        for(int j=0;j<d+1;j++)
            z+=comb[d][i]*comb[d][j]*pow_t[j]*pow_u[d-j]*pow_s[i]*pow_r[d-i]*X((d+1)*i+j);
    return z;
}
//#####################################################################
// Function Name
//#####################################################################
template<class TV, int d> std::string BEZIER_SPLINE_PATCH<TV,d>::
Name() const
{
    return Static_Name();
}
//#####################################################################
// Function Static_Name
//#####################################################################
template<class TV, int d> std::string BEZIER_SPLINE_PATCH<TV,d>::
Static_Name()
{
    return STRING_UTILITIES::string_sprintf("BEZIER_SPLINE_PATCH<VECTOR<T,%d> ,%d>",TV::dimension,d);
}
//#####################################################################
// Function Extension
//#####################################################################
template<class TV, int d> std::string BEZIER_SPLINE_PATCH<TV,d>::
Extension() const
{
    return Static_Extension();
}
//#####################################################################
// Function Static_Extension
//#####################################################################
template<class TV, int d> std::string BEZIER_SPLINE_PATCH<TV,d>::
Static_Extension()
{
    return "";
}
namespace PhysBAM{
template class BEZIER_SPLINE_PATCH<VECTOR<float,2>,3>;
template class BEZIER_SPLINE_PATCH<VECTOR<float,3>,3>;
template class BEZIER_SPLINE_PATCH<VECTOR<double,2>,3>;
template class BEZIER_SPLINE_PATCH<VECTOR<double,3>,3>;
template void Smooth_Fit<VECTOR<float,2> >(BEZIER_SPLINE_PATCH<VECTOR<float,2>,3>&,ARRAY_VIEW<VECTOR<float,2>,int>);
template void Smooth_Fit<VECTOR<double,2> >(BEZIER_SPLINE_PATCH<VECTOR<double,2>,3>&,ARRAY_VIEW<VECTOR<double,2>,int>);
}
