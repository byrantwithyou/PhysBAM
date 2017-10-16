//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Boundaries/BOUNDARY_MAC_GRID_PERIODIC.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MAC_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS_SPLINE.h>
#include <Hybrid_Methods/Projection/MPM_PROJECTION_SYSTEM.h>
#include <Hybrid_Methods/Projection/MPM_PROJECTION_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_MAC_EXAMPLE<TV>::
MPM_MAC_EXAMPLE(const STREAM_TYPE stream_type)
    :stream_type(stream_type),particles(*new MPM_PARTICLES<TV>),
    projection_system(*new MPM_PROJECTION_SYSTEM<TV>),
    sol(*new MPM_PROJECTION_VECTOR<TV>),rhs(*new MPM_PROJECTION_VECTOR<TV>),
    ghost(3),periodic_boundary(*new BOUNDARY_MAC_GRID_PERIODIC<TV,T>),extrap_type('r'),extrap_boundary(false),
    clamp_particles(false),
    use_affine(true),lag_Dp(false),flip(0),initial_time(0),
    last_frame(100),write_substeps_level(-1),substeps_delay_frame(-1),
    output_directory("output"),data_directory("../../Public_Data"),use_test_output(false),
    restart(0),dt(0),time(0),frame_dt((T)1/24),min_dt(0),max_dt(frame_dt),
    only_write_particles(false),cfl(1),
    solver_tolerance(std::numeric_limits<T>::epsilon()*10),solver_iterations(1000),
    threads(1),use_particle_volumes(false),use_constant_density(true),
    move_mass_inside(false),move_mass_inside_nearest(false),
    use_phi(false),rk_particle_order(0),
    use_massless_particles(false),use_multiphase_projection(false),use_bump(false),use_reseeding(false),
    use_periodic_test_shift(false),use_viscosity(false),
    debug_particles(*new DEBUG_PARTICLES<TV>),
    test_system(false),print_matrix(false)
{
    bc_type.Fill(BC_SLIP);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_MAC_EXAMPLE<TV>::
~MPM_MAC_EXAMPLE()
{
    delete &particles;
    delete &debug_particles;
    for(int i=0;i<TV::m;i++) delete weights(i);
    collision_objects.Delete_Pointers_And_Clean_Memory();
    fluid_walls.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_MAC_EXAMPLE<TV>::PHASE::
PHASE()
    :gather_scatter(0),levelset(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_MAC_EXAMPLE<TV>::PHASE::
~PHASE()
{
    delete gather_scatter;
    delete levelset;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void MPM_MAC_EXAMPLE<TV>::PHASE::
Initialize(const GRID<TV>& grid,const VECTOR<PARTICLE_GRID_WEIGHTS<TV>*,TV::m>& weights,
    int ghost,int threads)
{
    gather_scatter=new GATHER_SCATTER<TV>(grid,simulated_particles);
    gather_scatter->threads=threads;
    gather_scatter->face_weights=weights;
    gather_scatter->weights=weights(0);
    
    mass.Resize(grid.Domain_Indices(ghost));
    volume.Resize(grid.Domain_Indices(ghost));
    velocity.Resize(grid.Domain_Indices(ghost));
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void MPM_MAC_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    std::string f=LOG::sprintf("%d",frame);
    if(this->use_test_output){
        std::string file=LOG::sprintf("%s/%s-%03d.txt",output_directory.c_str(),test_output_prefix.c_str(),frame);
        OCTAVE_OUTPUT<T> oo(file.c_str());
        oo.Write("X",particles.X.Flattened());
        oo.Write("V",particles.V.Flattened());
        for(PHASE_ID i(0);i<phases.m;i++)
            oo.Write("u",phases(i).velocity.array);}

#pragma omp parallel
#pragma omp single
    {
#pragma omp task
        Write_To_File(stream_type,output_directory+"/common/grid",grid);
#pragma omp task
        Write_To_File(stream_type,LOG::sprintf("%s/%d/restart_data",output_directory.c_str(),frame),time);
#pragma omp task
        Write_To_File(stream_type,LOG::sprintf("%s/%d/mpm_particles",output_directory.c_str(),frame),particles);
#pragma omp task
        if(phases(PHASE_ID(0)).levelset)
            Write_To_File(stream_type,LOG::sprintf("%s/%d/levelset",output_directory.c_str(),frame),*phases(PHASE_ID(0)).levelset);

        if(!only_write_particles){
#pragma omp task
            Write_To_File(stream_type,LOG::sprintf("%s/%d/mac_velocities",output_directory.c_str(),frame),phases(PHASE_ID()).velocity);
#pragma omp task
            {
                GRID<TV> ghost_grid(grid.numbers_of_cells+2*ghost,grid.Ghost_Domain(ghost),true);
                for(int i=0;i<collision_objects.m;i++)
                    if(IMPLICIT_OBJECT<TV>* io=collision_objects(i)->Get_Implicit_Object(time))
                        Dump_Levelset(ghost_grid,*io,VECTOR<T,3>(0.7,0.3,0.3));
                debug_particles.Write_Debug_Particles(stream_type,output_directory,frame);
            }
        }
    }
    for(int i=0;i<write_output_files.m;i++) write_output_files(i)(frame);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void MPM_MAC_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=LOG::sprintf("%d",frame);
    Read_From_File(stream_type,LOG::sprintf("%s/%d/mpm_particles",output_directory.c_str(),frame),particles);
    Read_From_File(stream_type,LOG::sprintf("%s/%d/restart_data",output_directory.c_str(),frame),time);
    for(int i=0;i<read_output_files.m;i++) read_output_files(i)(frame);
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_EXAMPLE<TV>::
Potential_Energy(const T time) const
{
    // TODO:
    return 0;
}
//#####################################################################
// Function Apply_Forces
//#####################################################################
template<class TV> void MPM_MAC_EXAMPLE<TV>::
Apply_Forces(const T time)
{
    for(PHASE_ID p(0);p<phases.m;p++){
        PHASE& ph=phases(p);
        for(int i=0;i<ph.valid_flat_indices.m;i++){
            int k=ph.valid_flat_indices(i);
            FACE_INDEX<TV::m> f=ph.valid_indices(i);
            TV af=Compute_Analytic_Force(p,grid.Face(f),time);
            ph.velocity.array(k)+=dt*af(f.axis);
            ph.velocity.array(k)+=dt*gravity(f.axis);}}
}
//#####################################################################
// Function Compute_Analytic_Force
//#####################################################################
template<class TV> TV MPM_MAC_EXAMPLE<TV>::
Compute_Analytic_Force(PHASE_ID p,const TV& X,T time) const
{
    return TV();
}
//#####################################################################
// Function Set_Weights
//#####################################################################
template<class TV> void MPM_MAC_EXAMPLE<TV>::
Set_Weights(int order)
{
    for(int i=0;i<TV::m;i++){
        GRID<TV> face_grid=grid.Get_Face_MAC_Grid(i);
        if(order==1)
            weights(i)=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,1>(face_grid);
        else if(order==2)
            weights(i)=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,2>(face_grid);
        else if(order==3)
            weights(i)=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,3>(face_grid);
        else PHYSBAM_FATAL_ERROR("Unrecognized interpolation order");}
}
//#####################################################################
// Function Total_Particle_Linear_Momentum
//#####################################################################
template<class TV> TV MPM_MAC_EXAMPLE<TV>::
Total_Particle_Linear_Momentum() const
{
    TV result;
    for(PHASE_ID i(0);i<phases.m;i++)
#pragma omp parallel for
        for(int t=0;t<threads;t++){
            int a=t*phases(i).simulated_particles.m/threads;
            int b=(t+1)*phases(i).simulated_particles.m/threads;
            TV result_local;
            for(int k=a;k<b;k++){
                int p=phases(i).simulated_particles(k);
                result_local+=particles.mass(p)*particles.V(p);}
#pragma omp critical
            result+=result_local;}
    return result;
}
//#####################################################################
// Function Total_Grid_Linear_Momentum
//#####################################################################
template<class TV> TV MPM_MAC_EXAMPLE<TV>::
Total_Grid_Linear_Momentum() const
{
    TV result;
    for(PHASE_ID i(0);i<phases.m;i++){
        const PHASE& ph=phases(i);
#pragma omp parallel for
        for(int t=0;t<threads;t++){
            int a=t*ph.valid_flat_indices.m/threads;
            int b=(t+1)*ph.valid_flat_indices.m/threads;
            TV result_local;
            for(int k=a;k<b;k++){
                int j=ph.valid_flat_indices(k);
                result_local(ph.valid_indices(k).axis)+=ph.mass.array(j)*ph.velocity.array(j);}
#pragma omp critical
            result+=result_local;}}
    return result;
}
//#####################################################################
// Function Total_Particle_Angular_Momentum
//#####################################################################
template<class TV> typename TV::SPIN MPM_MAC_EXAMPLE<TV>::
Total_Particle_Angular_Momentum() const
{
    typename TV::SPIN result;
    for(PHASE_ID i(0);i<phases.m;i++){
        const PHASE& ph=phases(i);
#pragma omp parallel for
        for(int t=0;t<threads;t++){
            int a=t*ph.simulated_particles.m/threads;
            int b=(t+1)*ph.simulated_particles.m/threads;
            typename TV::SPIN result_local;
            for(int k=a;k<b;k++){
                int p=ph.simulated_particles(k);
                result_local+=particles.mass(p)*particles.X(p).Cross(particles.V(p));
                if(particles.store_B) result_local-=particles.mass(p)*particles.B(p).Contract_Permutation_Tensor();}
#pragma omp critical
            result+=result_local;}}
    return result;
}
//#####################################################################
// Function Total_Grid_Angular_Momentum
//#####################################################################
template<class TV> typename TV::SPIN MPM_MAC_EXAMPLE<TV>::
Total_Grid_Angular_Momentum(T dt) const
{
    typename TV::SPIN result;
    for(PHASE_ID i(0);i<phases.m;i++){
        const PHASE& ph=phases(i);
#pragma omp parallel for
        for(int t=0;t<threads;t++){
            int a=t*ph.valid_flat_indices.m/threads;
            int b=(t+1)*ph.valid_flat_indices.m/threads;
            typename TV::SPIN result_local;
            for(int k=a;k<b;k++){
                int i=ph.valid_flat_indices(k);
                TV X=location.array(i);
                result_local+=ph.mass.array(i)*TV::Cross_Product(X,ph.velocity.array(i)*TV::Axis_Vector(ph.valid_indices(k).axis));}
#pragma omp critical
            result+=result_local;}}
    return result;
}
//#####################################################################
// Function Total_Grid_Kinetic_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_EXAMPLE<TV>::
Total_Grid_Kinetic_Energy() const
{
    T result=0;
    for(PHASE_ID i(0);i<phases.m;i++){
        const PHASE& ph=phases(i);
#pragma omp parallel for reduction(+:result)
        for(int i=0;i<ph.valid_flat_indices.m;i++){
            int j=ph.valid_flat_indices(i);
            result+=(T).5*ph.mass.array(j)*sqr(ph.velocity.array(j));}}
    return result;
}
//#####################################################################
// Function Total_Particle_Kinetic_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_EXAMPLE<TV>::
Total_Particle_Kinetic_Energy() const
{
    T result=0;
    for(PHASE_ID i(0);i<phases.m;i++)
#pragma omp parallel for reduction(+:result)
        for(int k=0;k<phases(i).simulated_particles.m;k++){
            int p=phases(i).simulated_particles(k);
            T result_local=particles.mass(p)/2*particles.V(p).Magnitude_Squared();
            if(particles.store_B)
                for(int a=0;a<TV::m;a++){
                    SYMMETRIC_MATRIX<T,TV::m> D=lag_Dp?Dp_inv(a)(p).Inverse():Dp_inv(a)(p);
                    TV b=particles.B(p).Row(a);
                    result_local+=particles.mass(p)/2*b.Dot(D*b);}
            result+=result_local;}
    return result;
}
//#####################################################################
// Function Average_Particle_Mass
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_EXAMPLE<TV>::
Average_Particle_Mass() const
{
    T result=0;
    for(PHASE_ID i(0);i<phases.m;i++)
#pragma omp parallel for reduction(+:result)
        for(int k=0;k<phases(i).simulated_particles.m;k++)
            result+=particles.mass(phases(i).simulated_particles(k));
    return result/(T)particles.number;
}
//#####################################################################
// Function Add_Collision_Object
//#####################################################################
template<class TV> void MPM_MAC_EXAMPLE<TV>::
Add_Collision_Object(IMPLICIT_OBJECT<TV>* io,COLLISION_TYPE type,T friction,std::function<FRAME<TV>(T)> func_frame,std::function<TWIST<TV>(T)> func_twist)
{
    collision_objects.Append(new MPM_COLLISION_IMPLICIT_OBJECT<TV>(io,type,friction,func_frame,func_twist));
}
//#####################################################################
// Function Add_Fluid_Wall
//#####################################################################
template<class TV> void MPM_MAC_EXAMPLE<TV>::
Add_Fluid_Wall(IMPLICIT_OBJECT<TV>* io)
{
    fluid_walls.Append(io);
}
//#####################################################################
// Function Add_Callbacks
//#####################################################################
template<class TV> void MPM_MAC_EXAMPLE<TV>::
Add_Callbacks(bool is_begin,const char* func_name,std::function<void()> func)
{
    time_step_callbacks.Get_Or_Insert(func_name).y(is_begin).Append(func);
}
//#####################################################################
// Function Print_Grid_Stats
//#####################################################################
template<class TV> void MPM_MAC_EXAMPLE<TV>::
Print_Grid_Stats(const char* str)
{
    typename TV::SPIN am=Total_Grid_Angular_Momentum(dt);
    TV lm=Total_Grid_Linear_Momentum();
    T ke=Total_Grid_Kinetic_Energy();
    T pe=Potential_Energy(time);
    T te=ke+pe;
    LOG::cout<<str<<" linear  "<<"time " <<time<<" value "<<lm<<"  diff "<<(lm-last_linear_momentum)<<std::endl;
    LOG::cout<<str<<" angular "<<"time " <<time<<" value "<<am<<"  diff "<<(am-last_angular_momentum)<<std::endl;
    LOG::cout<<str<<" ke "<<"time " <<time<<" value "<<ke<<"  diff "<<(ke-last_grid_ke)<<std::endl;
    LOG::cout<<str<<" pe "<<"time " <<time<<" value "<<pe<<std::endl;
    LOG::cout<<str<<" total energy "<<"time " <<time<<" value "<<te<<" diff "<<(te-last_grid_te)<<std::endl;
    last_linear_momentum=lm;
    last_angular_momentum=am;
    last_grid_ke=ke;
    last_grid_te=te;
}
//#####################################################################
// Function Print_Grid_Stats
//#####################################################################
template<class TV> void MPM_MAC_EXAMPLE<TV>::
Print_Particle_Stats(const char* str)
{
    typename TV::SPIN am=Total_Particle_Angular_Momentum();
    TV lm=Total_Particle_Linear_Momentum();
    T ke=Total_Particle_Kinetic_Energy();
    T pe=Potential_Energy(time);
    T te=ke+pe;
    LOG::cout<<str<<" linear  "<<"time " <<time<<" value "<<lm<<"  diff "<<(lm-last_linear_momentum)<<std::endl;
    LOG::cout<<str<<" angular "<<"time " <<time<<" value "<<am<<"  diff "<<(am-last_angular_momentum)<<std::endl;
    LOG::cout<<str<<" ke "<<"time " <<time<<ke<<"  diff "<<(ke-last_part_ke)<<std::endl;
    LOG::cout<<str<<" pe "<<"time " <<time<<" value "<<pe<<std::endl;
    LOG::cout<<str<<" total energy "<<"time " <<time<<" value "<<te<<" diff "<<(te-last_part_te)<<std::endl;
    last_linear_momentum=lm;
    last_angular_momentum=am;
    last_part_ke=ke;
    last_part_te=te;
}
//#####################################################################
// Dump_Grid_ShiftTest
//#####################################################################
template<class TV> void MPM_MAC_EXAMPLE<TV>::
Dump_Grid_ShiftTest(const std::string& var_name,const ARRAY<T,FACE_INDEX<TV::m> >& arr)
{
    ARRAY<T,FACE_INDEX<TV::m> > arr_shift(arr.domain_indices);
    TV_INT id_shift;
    for(FACE_ITERATOR<TV> fit(grid);fit.Valid();fit.Next()){
        id_shift=wrap(fit.face.index+periodic_test_shift,TV_INT(),grid.counts);
        FACE_INDEX<TV::m> fid_shift(fit.face.axis,id_shift);
        arr_shift(fit.Full_Index())=arr(fid_shift);}
    LOG::printf("%s\n%P\n",var_name,arr_shift.array);
}
//#####################################################################
namespace PhysBAM{
template class MPM_MAC_EXAMPLE<VECTOR<float,2> >;
template class MPM_MAC_EXAMPLE<VECTOR<float,3> >;
template class MPM_MAC_EXAMPLE<VECTOR<double,2> >;
template class MPM_MAC_EXAMPLE<VECTOR<double,3> >;
}
