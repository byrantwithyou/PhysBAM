//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Deformables/Forces/LAGGED_FORCE.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MICROPOLAR_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_FORCE_HELPER.h>
#include <Hybrid_Methods/Forces/MPM_PLASTICITY_MODEL.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS_SPLINE.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_MICROPOLAR_EXAMPLE<TV>::
MPM_MICROPOLAR_EXAMPLE(const STREAM_TYPE stream_type)
    :stream_type(stream_type),particles(*new MPM_PARTICLES<TV>),
    debug_particles(*new DEBUG_PARTICLES<TV>),
    gather_scatter(*new GATHER_SCATTER<TV>(grid,simulated_particles)),
    force_helper(*new MPM_FORCE_HELPER<TV>(particles,quad_F_coeff))
{
    debug_particles.debug_particles.template Add_Array<T>("display_size");
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_MICROPOLAR_EXAMPLE<TV>::
~MPM_MICROPOLAR_EXAMPLE()
{
    delete &particles;
    delete &debug_particles;
    delete weights;
    delete &gather_scatter;
    delete &force_helper;
    forces.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void MPM_MICROPOLAR_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    std::string f=LOG::sprintf("%d",frame);
    if(this->use_test_output){
        std::string file=LOG::sprintf("%s/%s-%03d.txt",output_directory.c_str(),test_output_prefix.c_str(),frame);
        OCTAVE_OUTPUT<T> oo(file.c_str());
        oo.Write("X",particles.X.Flattened());
        oo.Write("V",particles.V.Flattened());
        oo.Write("u",velocity.array.Flattened());}

#pragma omp parallel
#pragma omp single
    {
#pragma omp task
        Write_To_File(stream_type,output_directory+"/common/grid",grid);
#pragma omp task
        Write_To_File(stream_type,LOG::sprintf("%s/%d/restart_data",output_directory.c_str(),frame),time);

        if(!only_write_particles){
#pragma omp task
            Write_To_File(stream_type,LOG::sprintf("%s/%d/centered_velocities",output_directory.c_str(),frame),velocity);
#pragma omp task
            Write_To_File(stream_type,LOG::sprintf("%s/%d/mpm_particles",output_directory.c_str(),frame),particles);
#pragma omp task
            {
                GRID<TV> ghost_grid(grid.numbers_of_cells+2*ghost,grid.Ghost_Domain(ghost),true);
                debug_particles.Write_Debug_Particles(stream_type,output_directory,frame);
            }
        }
    }
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void MPM_MICROPOLAR_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=LOG::sprintf("%d",frame);
    Read_From_File(stream_type,LOG::sprintf("%s/%d/mpm_particles",output_directory.c_str(),frame),particles);
    Read_From_File(stream_type,LOG::sprintf("%s/%d/restart_data",output_directory.c_str(),frame),time);
}
//#####################################################################
// Function Capture_Stress
//#####################################################################
template<class TV> void MPM_MICROPOLAR_EXAMPLE<TV>::
Capture_Stress()
{
    force_helper.Fn=particles.F;
    if(particles.store_S) force_helper.Sn=particles.S;
}
//#####################################################################
// Function Precompute_Forces
//#####################################################################
template<class TV> void MPM_MICROPOLAR_EXAMPLE<TV>::
Precompute_Forces(const T time,const T dt,const bool update_hessian)
{
    for(int i=0;i<forces.m;i++)
        forces(i)->Precompute(time,dt,true,update_hessian);
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MICROPOLAR_EXAMPLE<TV>::
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
template<class TV> void MPM_MICROPOLAR_EXAMPLE<TV>::
Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const
{
    for(int i=0;i<forces.m;i++)
        forces(i)->Add_Forces(F,time);
}
//#####################################################################
// Function Add_Hessian_Times
//#####################################################################
template<class TV> void MPM_MICROPOLAR_EXAMPLE<TV>::
Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& V,const T time) const
{
    for(int i=0;i<forces.m;i++)
        forces(i)->Add_Hessian_Times(F,V,time);
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class TV> int MPM_MICROPOLAR_EXAMPLE<TV>::
Add_Force(PARTICLE_GRID_FORCES<TV>& force)
{
    return forces.Append(&force);
}
//#####################################################################
// Function Set_Weights
//#####################################################################
template<class TV> void MPM_MICROPOLAR_EXAMPLE<TV>::
Set_Weights(int order)
{
    delete weights;
    if(order==1) gather_scatter.weights=weights=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,1>(grid,threads);
    else if(order==2) gather_scatter.weights=weights=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,2>(grid,threads);
    else if(order==3) gather_scatter.weights=weights=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,3>(grid,threads);
    else PHYSBAM_FATAL_ERROR("Unrecognized interpolation order");
}
//#####################################################################
namespace PhysBAM{
template class MPM_MICROPOLAR_EXAMPLE<VECTOR<float,2> >;
template class MPM_MICROPOLAR_EXAMPLE<VECTOR<float,3> >;
template class MPM_MICROPOLAR_EXAMPLE<VECTOR<double,2> >;
template class MPM_MICROPOLAR_EXAMPLE<VECTOR<double,3> >;
}
