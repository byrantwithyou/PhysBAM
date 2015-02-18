//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/PARTICLE_GRID_FORCES.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_EXAMPLE<TV>::
MPM_EXAMPLE(const STREAM_TYPE stream_type)
    :stream_type(stream_type),particles(*new MPM_PARTICLES<TV>),debug_particles(*new DEBUG_PARTICLES<TV>),
    rhs(*new MPM_KRYLOV_VECTOR<TV>(valid_grid_indices)),weights(0),
    gather_scatter(*new GATHER_SCATTER<TV>(simulated_particles,weights)),initial_time(0),last_frame(100),
    write_substeps_level(-1),substeps_delay_frame(-1),write_output_files(true),output_directory("output"),
    restart(0),dt(0),time(0),frame_dt((T)1/24),min_dt(0),max_dt(frame_dt),order(2),ghost(3),
    use_reduced_rasterization(false),use_affine(false),flip(0),cfl(1),newton_tolerance(-100),newton_iterations(-100),
    solver_tolerance(-100),solver_iterations(-100)
{
    PHYSBAM_ASSERT(grid.Is_MAC_Grid());
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_EXAMPLE<TV>::
~MPM_EXAMPLE()
{
    delete &particles;
    delete &debug_particles;
    delete &rhs;
    delete &weights;
    delete &gather_scatter;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    if(!write_output_files) return;
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);

    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",grid);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/particles",output_directory.c_str(),frame),particles);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/velocities",output_directory.c_str(),frame),velocity);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/density",output_directory.c_str(),frame),mass);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/restart_data",output_directory.c_str(),frame),time);

    debug_particles.Write_Debug_Particles(stream_type,output_directory,frame);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/particles",output_directory.c_str(),frame),particles);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/restart_data",output_directory.c_str(),frame),time);
}
//#####################################################################
// Function Precompute_Forces
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Precompute_Forces(const T time)
{
    for(int i=0;i<forces.m;i++)
        forces(i)->Precompute(time);
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_EXAMPLE<TV>::
Potential_Energy(const T time) const
{
    typename TV::SCALAR pe=0;
    for(int i=0;i<forces.m;i++)
        pe+=forces(i)->Potential_Energy(time);
    return pe;
}
//#####################################################################
// Function Add_Forces
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const
{
    for(int i=0;i<forces.m;i++)
        forces(i)->Add_Forces(F,time);
}
//#####################################################################
// Function Add_Hessian_Times
//#####################################################################
template<class TV> void MPM_EXAMPLE<TV>::
Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& V,const T time) const
{
    for(int i=0;i<forces.m;i++)
        forces(i)->Add_Hessian_Times(F,V,time);
}
//#####################################################################
namespace PhysBAM{
template class MPM_EXAMPLE<VECTOR<float,2> >;
template class MPM_EXAMPLE<VECTOR<float,3> >;
template class MPM_EXAMPLE<VECTOR<double,2> >;
template class MPM_EXAMPLE<VECTOR<double,3> >;
}
