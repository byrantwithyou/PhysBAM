//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/KD_TREE.h>
#include <Core/Matrices/FRAME.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Geometry/Basic_Geometry/CONE.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Basic_Geometry/HOURGLASS.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Basic_Geometry/TORUS.h>
#include <Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_DILATE.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_INTERSECTION.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_INVERT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_UNION.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Deformables/Collisions_And_Interactions/PINNING_FORCE.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_IMPLICIT_SPHERE.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <fstream>
#include "STANDARD_TESTS_3D.h"
namespace PhysBAM{
//#####################################################################
// Function Initialize_Implicit_Surface
//
// This was copied from DEFORMABLES_STANDARD_TESTS.cpp
// TODO: put this function somewhere more convenient (maybe as a constructor of LEVELSET_IMPLICIT_OBJECT?)
//#####################################################################
template<class T> LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >*
Initialize_Implicit_Surface(TRIANGULATED_SURFACE<T>& surface,int max_res)
{
    typedef VECTOR<int,3> TV_INT;
    LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >& undeformed_levelset=*LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >::Create();
    surface.Update_Bounding_Box();
    RANGE<VECTOR<T,3> > box=*surface.bounding_box;
    GRID<VECTOR<T,3> >& ls_grid=undeformed_levelset.levelset.grid;
    ARRAY<T,TV_INT>& phi=undeformed_levelset.levelset.phi;
    ls_grid=GRID<VECTOR<T,3> >::Create_Grid_Given_Cell_Size(box,box.Edge_Lengths().Max()/max_res,false,5);
    phi.Resize(ls_grid.Domain_Indices(3));
    LEVELSET_MAKER_UNIFORM<VECTOR<T,3> >::Compute_Level_Set(surface,ls_grid,3,phi);
    undeformed_levelset.Update_Box();
    return &undeformed_levelset;
}
//#####################################################################
// Function Levelset_From_File
//
// TODO: put this function where it belongs
//#####################################################################
template<class T> LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* 
Levelset_From_File(const std::string& filename,int max_resolution=200)
{
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    FILE_UTILITIES::Read_From_File(STREAM_TYPE(0.f),filename,*surface);
    LOG::printf("Read mesh: %d triangle, %d particles\n",surface->mesh.elements.m,surface->particles.number);
    surface->mesh.Initialize_Adjacent_Elements();
    surface->mesh.Initialize_Neighbor_Nodes();
    surface->mesh.Initialize_Incident_Elements();
    surface->Update_Bounding_Box();
    surface->Initialize_Hierarchy();
    surface->Update_Triangle_List();
    return Initialize_Implicit_Surface(*surface,max_resolution);
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,3> >::
STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :STANDARD_TESTS_BASE<TV>(stream_type_input,parse_args)
{
    parse_args.Parse();
    if(!this->override_output_directory) output_directory=LOG::sprintf("Test_3d_%i",test_number);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,3> >::
~STANDARD_TESTS()
{
    if(destroy) destroy();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Write_Output_Files(const int frame)
{
    if(write_output_files) write_output_files(frame);
    BASE::Write_Output_Files(frame);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Read_Output_Files(const int frame)
{
    if(read_output_files) read_output_files(frame);
    BASE::Read_Output_Files(frame);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Initialize()
{
    switch(test_number)
    {
        case 1:{ // rotating sphere
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            SPHERE<TV> sphere(TV(.5,.5,.5)*m,.3*m);
            VECTOR<T,3> angular_velocity(TV(0.4,0,0)/s);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return angular_velocity.Cross(X-sphere.center);},
                [=](const TV&){return MATRIX<T,3>::Cross_Product_Matrix(angular_velocity);}
                ,density,particles_per_cell);
            T total_mass=particles.mass.Sum();
            TV total_momentum=particles.V.Weighted_Sum(particles.mass);
            TV dV=total_momentum/total_mass;
            particles.V-=dV;
        } break;
        case 2:{ // Oscillating sphere
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            SPHERE<TV> sphere(TV(.5,.5,.5)*m,.3*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,[=](const TV&){return MATRIX<T,3>();},
                density,particles_per_cell);
        } break;
        case 3:{ // Freefall sphere
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            SPHERE<TV> sphere(TV(.5,.5,.5)*m,.3*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,[=](const TV&){return MATRIX<T,3>();},
                density,particles_per_cell);
        } break;
        default: PHYSBAM_FATAL_ERROR("test number not implemented");
    }
    if(forced_collision_type!=-1)
        for(int i=0;i<collision_objects.m;i++)
            collision_objects(i)->type=(COLLISION_TYPE)forced_collision_type;
}
//#####################################################################
// Function Begin_Frame
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Begin_Frame(const int frame)
{
    if(begin_frame) begin_frame(frame);
}
//#####################################################################
// Function End_Frame
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
End_Frame(const int frame)
{
    if(end_frame) end_frame(frame);
}
//#####################################################################

// Function Begin_Time_Step
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Begin_Time_Step(const T time)
{
    if(begin_time_step) begin_time_step(time);
}
//#####################################################################
// Function End_Time_Step
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
End_Time_Step(const T time)
{
    if(end_time_step) end_time_step(time);
}

template class STANDARD_TESTS<VECTOR<float,3> >;
template class STANDARD_TESTS<VECTOR<double,3> >;
}
