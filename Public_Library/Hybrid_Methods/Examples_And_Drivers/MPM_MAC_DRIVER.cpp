//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Poisson/PROJECTION_UNIFORM.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MAC_DRIVER.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MAC_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_PLASTIC_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <Hybrid_Methods/MPM_PLASTICITY_MODEL.h>
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
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Initialize()
{
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

    example.mass.Resize(example.grid.Domain_Indices(example.ghost));
    example.volume.Resize(example.grid.Domain_Indices(example.ghost));
    example.density.Resize(example.grid.Domain_Indices(example.ghost));
    example.velocity.Resize(example.grid.Domain_Indices(example.ghost));
    example.volume.Resize(example.grid.Domain_Indices(example.ghost));

    RANGE<TV_INT> range(example.grid.Cell_Indices(example.ghost));
    example.location.Resize(example.grid,example.ghost);
    for(FACE_ITERATOR<TV> it(example.grid,example.ghost);it.Valid();it.Next())
        example.location(it.Full_Index())=it.Location();

    // TODO
    // for(CELL_ITERATOR<TV> iterator(example.grid,example.ghost,GRID<TV>::GHOST_REGION);iterator.Valid();iterator.Next()){
    //     for(int k=0;k<fluid_walls.m;k++){
    //         if(fluid_walls(k)->Lazy_Inside(iterator.Location())){
    //             cell_soli....
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
    example.Begin_Time_Step(example.time);

    Update_Simulated_Particles();
    Print_Particle_Stats("particle state",example.dt);
    Update_Particle_Weights();
    example.gather_scatter.Prepare_Scatter(example.particles);
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

    example.End_Time_Step(example.time);
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    for(;current_frame<frame;current_frame++){
        LOG::SCOPE scope("FRAME","frame %d",current_frame+1);
        example.Begin_Frame(current_frame);
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
        example.End_Frame(current_frame);
        Write_Output_Files(++output_number);}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    if(level<=example.write_substeps_level){
        example.frame_title=title;
        LOG::printf("Writing substep [%s]: output_number=%i, time=%g, frame=%i, substep=%i\n",
            example.frame_title,output_number+1,time,current_frame,substep);
        Write_Output_Files(++output_number);
        example.frame_title="";}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Write_Output_Files(const int frame)
{
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
    for(int i=0;i<TV::m;++i)
        example.weights(i)->Update(example.particles.X);
}
//#####################################################################
// Function Particle_To_Grid
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Particle_To_Grid()
{
    MPM_PARTICLES<TV>& particles=example.particles;

#pragma omp parallel for
    for(int i=0;i<example.mass.array.m;i++){
        example.mass.array(i)=0;
        example.volume.array(i)=0;
        example.velocity.array(i)=0;}

    T scale=1;
    if(particles.store_B && !example.weights(0)->use_gradient_transfer)
        scale=example.weights(0)->Constant_Scalar_Inverse_Dp();
    bool use_gradient=false;
    ARRAY_VIEW<MATRIX<T,TV::m> > dV(particles.B);
    example.gather_scatter.template Scatter<int>(true,
        [this,scale,&particles,use_gradient,dV](int p,const PARTICLE_GRID_FACE_ITERATOR<TV>& it,int data)
        {
            T w=it.Weight();
            FACE_INDEX<TV::m> index=it.Index();
            example.mass(index)+=w*particles.mass(p);
            example.volume(index)+=w*particles.volume(p);
            T V=particles.V(p)(index.axis);
            if(example.use_affine){
                V+=dV(p).Row(index.axis).Dot(scale*(example.grid.Face(index)-particles.X(p)));}
            example.velocity(index)+=particles.mass(p)*w*V;
        });

    example.valid_indices.Remove_All();
    example.valid_flat_indices.Remove_All();

    for(FACE_ITERATOR<TV> it(example.grid,0);it.Valid();it.Next()){
        int i=example.mass.Standard_Index(it.Full_Index());
        if(example.mass.array(i)){
            example.valid_flat_indices.Append(i);
            example.valid_indices.Append(it.Full_Index());
            example.velocity.array(i)/=example.mass.array(i);}
        else example.velocity.array(i)=0;}
}
//#####################################################################
// Function Grid_To_Particle
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Grid_To_Particle()
{
    struct HELPER
    {
        TV V;
        MATRIX<T,TV::m> B;
    };

    MPM_PARTICLES<TV>& particles=example.particles;
    T dt=example.dt;
    example.gather_scatter.template Gather<HELPER>(true,
        [](int p,HELPER& h){h=HELPER();},
        [this,dt,&particles](int p,const PARTICLE_GRID_FACE_ITERATOR<TV>& it,HELPER& h)
        {
            T w=it.Weight();
            FACE_INDEX<TV::m> index=it.Index();
            T V=example.velocity(index);
            h.V(index.axis)+=w*V;
            if(example.use_affine){
                TV Z=example.grid.Face(index);
                h.B.Add_Row(index.axis,w*V*(Z-particles.X(p)));}
        },
        [this,dt,&particles](int p,HELPER& h)
        {
            if(particles.store_B) particles.B(p)=h.B;
            particles.X(p)+=dt*h.V;
            particles.V(p)=h.V;

            if(!example.grid.domain.Lazy_Inside(particles.X(p))) particles.valid(p)=false;
        });
}
//#####################################################################
// Function Compute_Volume_For_Face
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_DRIVER<TV>::
Compute_Volume_For_Face(const FACE_INDEX<TV::m>& face) const
{
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
    VECTOR<typename PARTICLE_GRID_WEIGHTS<TV>::SCRATCH,TV::m> scratch;
    for(int i=0;i<TV::m;i++)
        example.weights(i)->Compute(example.grid.Face(FACE_INDEX<TV::m>(i,TV_INT())),scratch(i),false);
    ARRAY<T,FACE_INDEX<TV::m> > face_fraction(example.grid,3); // Fraction of region of influence that is inside
    ARRAY<bool,FACE_INDEX<TV::m> > psi_N(example.grid);
    for(FACE_ITERATOR<TV> it(example.grid);it.Valid();it.Next()){
        TV X=it.Location();
        for(int i=0;i<example.collision_objects.m;i++){
            MPM_COLLISION_OBJECT<TV>* o=example.collision_objects(i);
            if(o->Phi(X,example.time)<0){
                psi_N(it.Full_Index())=true;
                example.velocity(it.Full_Index())=o->Velocity(X,example.time)(it.axis);
                break;}}
        if(!psi_N(it.Full_Index()))
            for(PARTICLE_GRID_ITERATOR<TV> pit(scratch(it.axis));pit.Valid();pit.Next())
                face_fraction(FACE_INDEX<TV::m>(it.axis,it.index+pit.Index()))+=pit.Weight();}

    ARRAY<int,TV_INT> cell_index(example.grid.Domain_Indices(1));
    int next_cell=0;
    for(CELL_ITERATOR<TV> it(example.grid,0);it.Valid();it.Next()){
        bool dirichlet=false,all_N=true;
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> face(a,it.index);
            for(int s=0;s<2;s++){
                if(!psi_N(face)){
                    all_N=false;
                    if(!example.mass(face)) dirichlet=true;}
                face.index(a)++;}}
        cell_index(it.index)=(all_N || dirichlet)?-1:next_cell++;}
    for(CELL_ITERATOR<TV> it(example.grid,1,GRID<TV>::GHOST_REGION);it.Valid();it.Next())
        cell_index(it.index)=-1;

    SPARSE_MATRIX_FLAT_MXN<T>& A=example.projection_system.A;
    A.Reset(next_cell);
    for(CELL_ITERATOR<TV> it(example.grid,0);it.Valid();it.Next()){
        int center_index=cell_index(it.index);
        if(center_index<0) continue;
    
        T diag=0;
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> face(a,it.index);
            for(int s=0;s<2;s++){
                TV_INT cell=it.index;
                cell(a)+=2*s-1;
                assert(face.Cell_Index(s)==cell);
                if(!psi_N(face)){
                    T entry=sqr(example.grid.one_over_dX(a))*Compute_Volume_For_Face(face)/example.mass(face);
                    diag+=entry;
                    int ci=cell_index(cell);
                    if(ci>=0) A.Append_Entry_To_Current_Row(ci,-entry);}
                face.index(a)++;}}
        A.Append_Entry_To_Current_Row(center_index,diag);
        A.Finish_Row();}
    A.Sort_Entries();
    A.Construct_Incomplete_Cholesky_Factorization();

    SPARSE_MATRIX_FLAT_MXN<T>& G=example.projection_system.gradient;
    ARRAY<T>& M=example.projection_system.mass;
    ARRAY<FACE_INDEX<TV::m> >& faces=example.projection_system.faces;
    M.Remove_All();
    faces.Remove_All();
    G.Reset(next_cell);
    for(FACE_ITERATOR<TV> it(example.grid);it.Valid();it.Next()){
        if(psi_N(it.Full_Index())) continue;
        T mass=example.mass(it.Full_Index());
        if(!mass) continue;
        int c0=cell_index(it.First_Cell_Index());
        int c1=cell_index(it.Second_Cell_Index());
        if(c0>=0) G.Append_Entry_To_Current_Row(c0,-example.grid.one_over_dX(it.axis));
        if(c1>=0) G.Append_Entry_To_Current_Row(c1,example.grid.one_over_dX(it.axis));
        if(c0<0 && c1<0) continue;
        faces.Append(it.Full_Index());
        M.Append(mass/Compute_Volume_For_Face(it.Full_Index()));
        G.Finish_Row();}
    G.Sort_Entries();
}
//#####################################################################
// Function Pressure_Projection
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Pressure_Projection()
{
    Compute_Poisson_Matrix();
    example.sol.v.Resize(example.projection_system.A.m);
    example.rhs.v.Resize(example.projection_system.A.m);

    ARRAY<T> tmp(example.projection_system.gradient.m);
    for(int i=0;i<tmp.m;i++)
        tmp(i)=example.velocity(example.projection_system.faces(i));
    example.projection_system.gradient.Transpose_Times(tmp,example.rhs.v);
    
    example.projection_system.Test_System(example.sol);
    
    CONJUGATE_GRADIENT<T> cg;
    cg.finish_before_indefiniteness=true;
    cg.relative_tolerance=true;
    cg.print_diagnostics=true;
    cg.print_residuals=true;
    bool converged=cg.Solve(example.projection_system,example.sol,example.rhs,
        example.av,example.solver_tolerance,0,example.solver_iterations);
    if(!converged) LOG::printf("SOLVER DID NOT CONVERGE.\n");

    example.projection_system.gradient.Times(example.sol.v,tmp);
    for(int i=0;i<tmp.m;i++)
        example.velocity(example.projection_system.faces(i))-=tmp(i)/example.projection_system.mass(i);
}
//#####################################################################
// Function Apply_Forces
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Apply_Forces()
{
    for(int i=0;i<example.valid_flat_indices.m;i++){
        int k=example.valid_flat_indices(i);
        example.velocity.array(k)+=example.dt*example.gravity(example.valid_indices(i).axis);}
}
//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_DRIVER<TV>::
Compute_Dt() const
{
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
    example.simulated_particles.Remove_All();
    for(int p=0;p<example.particles.number;p++)
        if(example.particles.valid(p))
            example.simulated_particles.Append(p);
    // example.simulated_particles=IDENTITY_ARRAY<>(example.particles.X.m);
    example.particle_is_simulated.Remove_All();
    example.particle_is_simulated.Resize(example.particles.X.m);
    example.particle_is_simulated.Subset(example.simulated_particles).Fill(true);
}
//#####################################################################
// Function Print_Grid_Stats
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Print_Grid_Stats(const char* str,T dt)
{
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
namespace PhysBAM{
template class MPM_MAC_DRIVER<VECTOR<float,2> >;
template class MPM_MAC_DRIVER<VECTOR<float,3> >;
template class MPM_MAC_DRIVER<VECTOR<double,2> >;
template class MPM_MAC_DRIVER<VECTOR<double,3> >;
}
