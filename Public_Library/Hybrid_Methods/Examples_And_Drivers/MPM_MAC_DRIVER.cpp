//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/FINE_TIMER.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Core/Matrices/SPARSE_MATRIX_THREADED_CONSTRUCTION.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Tools/Parallel_Computation/APPEND_HOLDER.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/CELL_ITERATOR_THREADED.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR_THREADED.h>
#include <Grid_PDE/Poisson/PROJECTION_UNIFORM.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MAC_DRIVER.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MAC_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_PLASTIC_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Forces/MPM_PLASTICITY_MODEL.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <Hybrid_Methods/Projection/MPM_PROJECTION_SYSTEM.h>
#include <Hybrid_Methods/Projection/MPM_PROJECTION_VECTOR.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_VECTOR.h>
#include <Hybrid_Methods/System/MPM_OBJECTIVE.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
using namespace PhysBAM;
namespace{
    template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
    {
        ((MPM_MAC_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
    }
};
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_MAC_DRIVER<TV>::
MPM_MAC_DRIVER(MPM_MAC_EXAMPLE<TV>& example)
    :example(example)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_MAC_DRIVER<TV>::
~MPM_MAC_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Execute_Main_Program
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
    FINE_TIMER::Dump_Timing_Info();
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Initialize()
{
    TIMER_SCOPE_FUNC;
    LOG::cout<<std::setprecision(16)<<std::endl;
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.substeps_delay_frame<0?example.write_substeps_level:-1);

    // setup time
    output_number=current_frame=example.restart;

    example.Initialize();

    PHYSBAM_ASSERT(example.grid.Is_MAC_Grid());
    if(example.restart)
        example.Read_Output_Files(example.restart);

    example.particles.Store_B(example.use_affine);
    example.particles.Store_C(false);
    PHYSBAM_ASSERT(!example.particles.store_B || !example.particles.store_C);

    for(PHASE_ID i(0);i<example.phases.m;i++)
        example.phases(i).Initialize(example.grid,example.ghost,example.threads);

    RANGE<TV_INT> range(example.grid.Cell_Indices(example.ghost));
    example.location.Resize(example.grid,example.ghost);
    for(FACE_ITERATOR<TV> it(example.grid,example.ghost);it.Valid();it.Next())
        example.location(it.Full_Index())=it.Location();

    Update_Simulated_Particles();

    if(!example.restart) Write_Output_Files(0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after init",0,1);
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Advance_One_Time_Step()
{
    TIMER_SCOPE_FUNC;
    if(example.begin_time_step) example.begin_time_step(example.time);

    Update_Simulated_Particles();
    Print_Particle_Stats("particle state",example.dt);
    Update_Particle_Weights();
    Prepare_Scatter();
    Particle_To_Grid();
    Print_Grid_Stats("after particle to grid",example.dt);
    Print_Energy_Stats("after particle to grid");
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after particle to grid",0,1);
    Apply_Forces();
    Print_Grid_Stats("after forces",example.dt);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after forces",0,1);
    Pressure_Projection();
    Print_Grid_Stats("after projection",example.dt);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after projection",0,1);
    Grid_To_Particle();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after grid to particle",0,1);

    if(example.end_time_step) example.end_time_step(example.time);
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    TIMER_SCOPE_FUNC;
    for(;current_frame<frame;current_frame++){
        LOG::SCOPE scope("FRAME","frame %d",current_frame+1);
        if(example.begin_frame) example.begin_frame(current_frame);
        if(example.substeps_delay_frame==current_frame)
            DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);
        T time_at_frame=example.time+example.frame_dt;
        bool done=false;
        for(int substep=0;!done;substep++){
            LOG::SCOPE scope("SUBSTEP","substep %d",substep+1);
            example.dt=Compute_Dt();
            example.dt=clamp(example.dt,example.min_dt,example.max_dt);
            T next_time=example.time+example.dt;
            if(next_time>time_at_frame){
                next_time=time_at_frame;
                done=true;}
            else if(next_time+example.dt>time_at_frame) next_time=(example.time+time_at_frame)/2;
            example.dt=next_time-example.time;
            LOG::cout<<"substep dt: "<<example.dt<<std::endl;

            Advance_One_Time_Step();
            example.time=next_time;}
        if(example.end_frame) example.end_frame(current_frame);
        Write_Output_Files(++output_number);}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    TIMER_SCOPE_FUNC;
    if(level<=example.write_substeps_level){
        example.frame_title=title;
        LOG::printf("Writing substep [%s]: output_number=%i, time=%g, frame=%i, substep=%i\n",
            example.frame_title,output_number+1,example.time,current_frame,substep);
        Write_Output_Files(++output_number);
        example.frame_title="";}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Write_Output_Files(const int frame)
{
    TIMER_SCOPE_FUNC;
    LOG::SCOPE scope("Write_Output_Files");
    FILE_UTILITIES::Create_Directory(example.output_directory);
    FILE_UTILITIES::Create_Directory(example.output_directory+LOG::sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+LOG::sprintf("/%d/frame_title",frame),example.frame_title);
    if(frame==0)
        FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/first_frame",frame,"\n");
    example.Write_Output_Files(frame);
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
// Function Update_Particle_Weights
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Update_Particle_Weights()
{
    TIMER_SCOPE_FUNC;
    for(int i=0;i<TV::m;++i)
        example.weights(i)->Update(example.particles.X);
}
//#####################################################################
// Function Particle_To_Grid
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Particle_To_Grid(PHASE& ph) const
{
    TIMER_SCOPE_FUNC;
    const MPM_PARTICLES<TV>& particles=example.particles;

#pragma omp parallel for
    for(int i=0;i<ph.mass.array.m;i++){
        ph.mass.array(i)=0;
        ph.volume.array(i)=0;
        ph.velocity.array(i)=0;}

    T scale=1;
    if(particles.store_B && !example.weights(0)->use_gradient_transfer)
        scale=example.weights(0)->Constant_Scalar_Inverse_Dp();
    bool use_gradient=false;
    ARRAY_VIEW<MATRIX<T,TV::m> > dV(particles.B);
    ph.gather_scatter->template Scatter<int>(true,
        [this,scale,&particles,use_gradient,dV,&ph](int p,const PARTICLE_GRID_FACE_ITERATOR<TV>& it,int data)
        {
            T w=it.Weight();
            FACE_INDEX<TV::m> index=it.Index();
            ph.mass(index)+=w*particles.mass(p);
            ph.volume(index)+=w*particles.volume(p);
            T V=particles.V(p)(index.axis);
            if(example.use_affine){
                V+=dV(p).Row(index.axis).Dot(scale*(example.grid.Face(index)-particles.X(p)));}
            ph.velocity(index)+=particles.mass(p)*w*V;
        });

    ph.valid_indices.Remove_All();
    ph.valid_flat_indices.Remove_All();

    APPEND_HOLDER<int> flat_h(ph.valid_flat_indices);
    APPEND_HOLDER<FACE_INDEX<TV::m> > indices_h(ph.valid_indices);
#pragma omp parallel
    {
#pragma omp single
        {
            flat_h.Init();
            indices_h.Init();
        }
#pragma omp barrier
        ARRAY<int>& flat_t=flat_h.Array();
        ARRAY<FACE_INDEX<TV::m> >& indices_t=indices_h.Array();
        for(FACE_ITERATOR_THREADED<TV> it(example.grid,0);it.Valid();it.Next()){
            int i=ph.mass.Standard_Index(it.Full_Index());
            if(ph.mass.array(i)){
                flat_t.Append(i);
                indices_t.Append(it.Full_Index());
                ph.velocity.array(i)/=ph.mass.array(i);}
            else ph.velocity.array(i)=0;}
    }
    flat_h.Combine();
    indices_h.Combine();
    if(example.flip) ph.velocity_save=ph.velocity;
}
//#####################################################################
// Function Particle_To_Grid
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Particle_To_Grid()
{
    for(PHASE_ID i(0);i<example.phases.m;i++)
        Particle_To_Grid(example.phases(i));
}
//#####################################################################
// Function Grid_To_Particle
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Grid_To_Particle(const PHASE& ph)
{
    TIMER_SCOPE_FUNC;
    struct HELPER
    {
        TV V;
        TV flip_V;
        MATRIX<T,TV::m> B;
    };

    MPM_PARTICLES<TV>& particles=example.particles;
    T dt=example.dt;
    bool use_flip=example.flip;
    ph.gather_scatter->template Gather<HELPER>(true,
        [](int p,HELPER& h){h=HELPER();},
        [this,dt,&particles,use_flip,&ph](int p,const PARTICLE_GRID_FACE_ITERATOR<TV>& it,HELPER& h)
        {
            T w=it.Weight();
            FACE_INDEX<TV::m> index=it.Index();
            T V=ph.velocity(index);
            h.V(index.axis)+=w*V;
            if(use_flip) h.flip_V(index.axis)+=w*(V-ph.velocity_save(index));
            if(example.use_affine){
                TV Z=example.grid.Face(index);
                h.B.Add_Row(index.axis,w*V*(Z-particles.X(p)));}
        },
        [this,dt,&particles,use_flip](int p,HELPER& h)
        {
            if(particles.store_B) particles.B(p)=h.B;
            particles.X(p)+=dt*h.V;
            if(use_flip) particles.V(p)=(1-example.flip)*h.V+example.flip*(particles.V(p)+h.flip_V);
            else particles.V(p)=h.V;

            if(!example.grid.domain.Lazy_Inside(particles.X(p))) particles.valid(p)=false;
        });
}
//#####################################################################
// Function Grid_To_Particle
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Grid_To_Particle()
{
    for(PHASE_ID i(0);i<example.phases.m;i++)
        Grid_To_Particle(example.phases(i));
}
//#####################################################################
// Function Compute_Volume_For_Face
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_DRIVER<TV>::
Compute_Volume_For_Face(const FACE_INDEX<TV::m>& face) const
{
    TIMER_SCOPE_FUNC;
    T total_volume=0;
    int num=0;
    int width=example.weights(face.axis)->Order()+1;
    RANGE<TV_INT> range(TV_INT(),TV_INT()+width*2);
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
        TV X=example.grid.dX*(TV(it.index-width)*(T).5+(T).25);
        TV Y=example.grid.Face(face)+X;
        bool in=false;
        for(int i=0;i<example.collision_objects.m;i++){
            MPM_COLLISION_OBJECT<TV>* o=example.collision_objects(i);
            if(o->Phi(Y,example.time)<0){
                Add_Debug_Particle(example.grid.Face(face),VECTOR<T,3>(.5,.5,.5));
                in=true;
                break;}}
        if(in) continue;
        T weight=example.weights(face.axis)->Weight(X);
        num++;
        total_volume+=weight;}
    return total_volume;
}
//#####################################################################
// Function Compute_Poisson_Matrix
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Compute_Poisson_Matrix()
{
    TIMER_SCOPE_FUNC;
    VECTOR<typename PARTICLE_GRID_WEIGHTS<TV>::SCRATCH,TV::m> scratch;
    for(int i=0;i<TV::m;i++)
        example.weights(i)->Compute(example.grid.Face(FACE_INDEX<TV::m>(i,TV_INT())),scratch(i),false);
//    ARRAY<T,FACE_INDEX<TV::m> > face_fraction(example.grid,0,false); // Fraction of region of influence that is inside
    ARRAY<bool,FACE_INDEX<TV::m> > psi_N(example.grid,3,false);

#pragma omp parallel
    for(FACE_ITERATOR_THREADED<TV> it(example.grid,3);it.Valid();it.Next()){
        TV X=it.Location();
        bool N=false;
        for(int i=0;i<example.collision_objects.m;i++){
            MPM_COLLISION_OBJECT<TV>* o=example.collision_objects(i);
            if(o->Phi(X,example.time)<0){
                N=true;
                example.phases(PHASE_ID()).velocity(it.Full_Index())=o->Velocity(X,example.time)(it.axis);
                break;}}
        psi_N(it.Full_Index())=N;}

// #pragma omp parallel
//     for(FACE_ITERATOR_THREADED<TV> it(example.grid);it.Valid();it.Next()){
//         T ff=0;
//         for(PARTICLE_GRID_ITERATOR<TV> pit(scratch(it.axis));pit.Valid();pit.Next())
//             if(!psi_N(FACE_INDEX<TV::m>(it.axis,it.index+pit.Index())))
//                 ff+=pit.Weight();
//         face_fraction(it.Full_Index())=ff;}

    ARRAY<int,TV_INT> cell_index(example.grid.Domain_Indices(1),false);
    int next_cell=0;
    for(CELL_ITERATOR<TV> it(example.grid,0);it.Valid();it.Next()){
        bool dirichlet=false,all_N=true;
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> face(a,it.index);
            for(int s=0;s<2;s++){
                if(!psi_N(face)){
                    all_N=false;
                    if(!example.phases(PHASE_ID()).mass(face)) dirichlet=true;}
                face.index(a)++;}}
        cell_index(it.index)=(all_N || dirichlet)?-1:next_cell++;}
    for(CELL_ITERATOR<TV> it(example.grid,1,GRID<TV>::GHOST_REGION);it.Valid();it.Next())
        cell_index(it.index)=-1;

    example.projection_system.A.Reset(next_cell);
    ARRAY<int> tmp0,tmp1;
#pragma omp parallel
    {
        SPARSE_MATRIX_THREADED_CONSTRUCTION<T> helper(example.projection_system.A,tmp0,tmp1);
        for(CELL_ITERATOR_THREADED<TV> it(example.grid,0);it.Valid();it.Next()){
            int center_index=cell_index(it.index);
            if(center_index<0) continue;
    
            T diag=0;
            helper.Start_Row();
            for(int a=0;a<TV::m;a++){
                FACE_INDEX<TV::m> face(a,it.index);
                for(int s=0;s<2;s++){
                    TV_INT cell=it.index;
                    cell(a)+=2*s-1;
                    assert(face.Cell_Index(s)==cell);
                    if(!psi_N(face)){
                        T entry=sqr(example.grid.one_over_dX(a))/example.phases(PHASE_ID()).mass(face);
                        if(example.use_particle_volumes) entry*=example.phases(PHASE_ID()).volume(face);
                        diag+=entry;
                        int ci=cell_index(cell);
                        if(ci>=0) helper.Add_Entry(ci,-entry);}
                    face.index(a)++;}}
            helper.Add_Entry(center_index,diag);}
        helper.Finish();
    }
    if(example.projection_system.use_preconditioner)
        example.projection_system.A.Construct_Incomplete_Cholesky_Factorization();

    example.projection_system.mass.Remove_All();
    example.projection_system.faces.Remove_All();
    example.projection_system.gradient.Reset(next_cell);
    APPEND_HOLDER<T> mass_h(example.projection_system.mass);
    APPEND_HOLDER<FACE_INDEX<TV::m> > faces_h(example.projection_system.faces);
    ARRAY<int> tmp2,tmp3;
#pragma omp parallel
    {
#pragma omp single
        {
            mass_h.Init();
            faces_h.Init();
        }
#pragma omp barrier
        ARRAY<T>& mass_t=mass_h.Array();
        ARRAY<FACE_INDEX<TV::m> >& faces_t=faces_h.Array();
        SPARSE_MATRIX_THREADED_CONSTRUCTION<T> G_helper(example.projection_system.gradient,tmp2,tmp3);
        for(FACE_ITERATOR_THREADED<TV> it(example.grid);it.Valid();it.Next()){
            if(psi_N(it.Full_Index())) continue;
            T mass=example.phases(PHASE_ID()).mass(it.Full_Index());
            if(!mass) continue;
            int c0=cell_index(it.First_Cell_Index());
            int c1=cell_index(it.Second_Cell_Index());
            if(c0<0 && c1<0) continue;
            G_helper.Start_Row(); // cannot start a row if no elements in it
            if(c0>=0) G_helper.Add_Entry(c0,-example.grid.one_over_dX(it.axis));
            if(c1>=0) G_helper.Add_Entry(c1,example.grid.one_over_dX(it.axis));
            faces_t.Append(it.Full_Index());
            if(example.use_particle_volumes) mass/=example.phases(PHASE_ID()).volume(it.Full_Index());
            mass_t.Append(mass);
        }
        G_helper.Finish();
    }
    mass_h.Combine();
    faces_h.Combine();
}
//#####################################################################
// Function Pressure_Projection
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Pressure_Projection()
{
    TIMER_SCOPE_FUNC;
    Compute_Poisson_Matrix();
    example.sol.v.Resize(example.projection_system.A.m);
    example.rhs.v.Resize(example.projection_system.A.m);

    ARRAY<T> tmp(example.projection_system.gradient.m);
    for(int i=0;i<tmp.m;i++)
        tmp(i)=example.phases(PHASE_ID()).velocity(example.projection_system.faces(i));
    example.projection_system.gradient.Transpose_Times(tmp,example.rhs.v);
    
    if(example.test_system) example.projection_system.Test_System(example.sol);
    
    CONJUGATE_GRADIENT<T> cg;
    cg.finish_before_indefiniteness=true;
    cg.relative_tolerance=true;
    bool converged=cg.Solve(example.projection_system,example.sol,example.rhs,
        example.av,example.solver_tolerance,0,example.solver_iterations);
    if(!converged) LOG::printf("SOLVER DID NOT CONVERGE.\n");

    example.projection_system.gradient.Times(example.sol.v,tmp);
    for(int i=0;i<tmp.m;i++)
        example.phases(PHASE_ID()).velocity(example.projection_system.faces(i))-=tmp(i)/example.projection_system.mass(i);
}
//#####################################################################
// Function Apply_Forces
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Apply_Forces()
{
    TIMER_SCOPE_FUNC;
    for(PHASE_ID p(0);p<example.phases.m;p++){
        PHASE& ph=example.phases(p);
        for(int i=0;i<ph.valid_flat_indices.m;i++){
            int k=ph.valid_flat_indices(i);
            ph.velocity.array(k)+=example.dt*example.gravity(ph.valid_indices(i).axis);}}
}
//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_DRIVER<TV>::
Compute_Dt() const
{
    TIMER_SCOPE_FUNC;
    T critical_speed=example.cfl*example.grid.dX.Min()/example.max_dt;
    T v=Grid_V_Upper_Bound();
    return (v>critical_speed)?(example.cfl*example.grid.dX.Min()/v):example.max_dt;
}
//#####################################################################
// Function Max_Particle_Speed
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_DRIVER<TV>::
Max_Particle_Speed() const
{
    TIMER_SCOPE_FUNC;
    T v2=0;
#pragma omp parallel for reduction(max:v2)
    for(int k=0;k<example.simulated_particles.m;k++){
        int p=example.simulated_particles(k);
        v2=max(v2,example.particles.V(p).Magnitude_Squared());}
    return sqrt(v2);
}
//#####################################################################
// Function Grid_V_Upper_Bound
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_DRIVER<TV>::
Grid_V_Upper_Bound() const
{
    TIMER_SCOPE_FUNC;
    if(!example.use_affine || !example.weights(0)->constant_scalar_inertia_tensor) return Max_Particle_Speed();
    T result=0;
    T xi=(T)6*sqrt((T)TV::m)*example.grid.one_over_dX.Min();
#pragma omp parallel for reduction(max:result)
    for(int k=0;k<example.simulated_particles.m;k++){
        int p=example.simulated_particles(k);
        T v=example.particles.V(p).Magnitude();
        if(example.particles.store_B) v+=example.particles.B(p).Frobenius_Norm()*xi;
        result=max(result,v);}
    return result;
}
//#####################################################################
// Function Update_Simulated_Particles
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Update_Simulated_Particles()
{
    TIMER_SCOPE_FUNC;
    example.simulated_particles.Remove_All();
    for(int p=0;p<example.particles.number;p++)
        if(example.particles.valid(p))
            example.simulated_particles.Append(p);
}
//#####################################################################
// Function Print_Grid_Stats
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Print_Grid_Stats(const char* str,T dt)
{
    TIMER_SCOPE_FUNC;
    if(!example.print_stats) return;
    typename TV::SPIN am=example.Total_Grid_Angular_Momentum(dt);
    TV lm=example.Total_Grid_Linear_Momentum();
    T ke=example.Total_Grid_Kinetic_Energy();
    LOG::cout<<str<<" linear  "<<"time " <<example.time<<" value "<<lm<<"  diff "<<(lm-example.last_linear_momentum)<<std::endl;
    LOG::cout<<str<<" angular "<<"time " <<example.time<<" value "<<am<<"  diff "<<(am-example.last_angular_momentum)<<std::endl;
    LOG::cout<<str<<" ke "<<"time " <<example.time<<" value "<<ke<<"  diff "<<(ke-example.last_grid_ke)<<std::endl;
    example.last_linear_momentum=lm;
    example.last_angular_momentum=am;
    example.last_grid_ke=ke;
}
//#####################################################################
// Function Print_Grid_Stats
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Print_Particle_Stats(const char* str,T dt)
{
    TIMER_SCOPE_FUNC;
    if(!example.print_stats) return;
    typename TV::SPIN am=example.Total_Particle_Angular_Momentum();
    TV lm=example.Total_Particle_Linear_Momentum();
    // T ke=example.Total_Particle_Kinetic_Energy();
    LOG::cout<<str<<" linear  "<<"time " <<example.time<<" value "<<lm<<"  diff "<<(lm-example.last_linear_momentum)<<std::endl;
    LOG::cout<<str<<" angular "<<"time " <<example.time<<" value "<<am<<"  diff "<<(am-example.last_angular_momentum)<<std::endl;
    // LOG::cout<<str<<" ke "<<ke<<"  diff "<<(ke-example.last_grid_ke)<<std::endl;
    example.last_linear_momentum=lm;
    example.last_angular_momentum=am;
    // example.last_grid_ke=ke;
}
//#####################################################################
// Function Print_Energy_Stats
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Print_Energy_Stats(const char* str)
{
    TIMER_SCOPE_FUNC;
    if(!example.print_stats) return;
    T ke=example.Total_Grid_Kinetic_Energy();
    T ke2=example.Total_Particle_Kinetic_Energy();
    T pe=example.Potential_Energy(example.time);
    T te=ke+pe;
    LOG::cout<<str<<" kinetic  "<<"time " <<example.time<<" value "<<ke<<std::endl;
    LOG::cout<<str<<" potential "<<"time " <<example.time<<" value "<<pe<<std::endl;
    LOG::cout<<str<<" total energy "<<"time " <<example.time<<" value "<<te<<" diff "<<(te-example.last_te)<<std::endl;
    LOG::cout<<str<<" particle total energy "<<"time " <<example.time<<" value "<<(ke2+pe)<<std::endl;
    example.last_te=te;
}
//#####################################################################
// Function Prepare_Scatter
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Prepare_Scatter()
{
    for(PHASE_ID i(0);i<example.phases.m;i++)
        example.phases(i).gather_scatter->Prepare_Scatter(example.particles);
}
//#####################################################################
namespace PhysBAM{
template class MPM_MAC_DRIVER<VECTOR<float,2> >;
template class MPM_MAC_DRIVER<VECTOR<float,3> >;
template class MPM_MAC_DRIVER<VECTOR<double,2> >;
template class MPM_MAC_DRIVER<VECTOR<double,3> >;
}
