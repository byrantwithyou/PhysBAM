//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/FINE_TIMER.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Core/Matrices/SPARSE_MATRIX_THREADED_CONSTRUCTION.h>
#include <Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <Tools/Parallel_Computation/APPEND_HOLDER.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/CELL_ITERATOR_THREADED.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR_THREADED.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_MAC.h>
#include <Grid_PDE/Poisson/PROJECTION_UNIFORM.h>
#include <Geometry/Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <Geometry/Level_Sets/LEVELSET.h>
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

    for(PHASE_ID i(0);i<example.phases.m;i++){
        PHASE& ph=example.phases(i);
        ph.Initialize(example.grid,example.weights,example.ghost,example.threads);
        ph.phi.Resize(example.grid.Domain_Indices(example.ghost));
        if(example.use_phi) ph.levelset=new LEVELSET<TV>(example.grid,ph.phi,example.ghost);}

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
    Build_Level_Sets();
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
    Create_Directory(example.output_directory);
    Create_Directory(example.output_directory+LOG::sprintf("/%d",frame));
    Create_Directory(example.output_directory+"/common");
    Write_To_Text_File(example.output_directory+LOG::sprintf("/%d/frame_title",frame),example.frame_title);
    if(frame==0)
        Write_To_Text_File(example.output_directory+"/common/first_frame",frame,"\n");
    example.Write_Output_Files(frame);
    Write_To_Text_File(example.output_directory+"/common/last_frame",frame,"\n");
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
    // TODO: move mass, momentum, volume inside
    if(example.move_mass_inside)
        Move_Mass_Momentum_Inside(ph);
    else if(example.move_mass_inside_nearest)
        Move_Mass_Momentum_Inside_Nearest(ph);

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
// Function Build_Level_Sets
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Build_Level_Sets()
{
    if(!example.use_phi) return;
    for(PHASE_ID i(0);i<example.phases.m;i++)
        Build_Level_Sets(example.phases(i));
    if(Value(example.phases.m)==1) return;
    for(CELL_ITERATOR<TV> it(example.grid,example.ghost);it.Valid();it.Next()){
        T min1=FLT_MAX,min2=FLT_MAX;
        for(PHASE_ID i(0);i<example.phases.m;i++){
            T p=example.phases(i).phi(it.index);
            if(p<min1){min2=min1;min1=p;}
            else if(p<min2) min2=p;}
        T shift=(T).5*(min2+min1);
        for(PHASE_ID i(0);i<example.phases.m;i++)
            example.phases(i).phi(it.index)-=shift;}
}
//#####################################################################
// Function Build_Level_Sets
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Build_Level_Sets(PHASE& ph)
{
    PHYSBAM_ASSERT(ph.levelset);
    T dx=example.grid.dX.Max();
ph.phi.array.Fill(3*dx);
    RANGE<TV_INT> grid_domain=example.grid.Domain_Indices(example.ghost);
    for(int k=0;k<ph.simulated_particles.m;k++){
        int p=ph.simulated_particles(k);
        T r=0.36*dx;
        T influence_bound=r+dx*(T)1.1;
        TV X=example.particles.X(p);
        RANGE<TV> bound(X-influence_bound,X+influence_bound);
        RANGE<TV_INT> range=example.grid.Clamp_To_Cell(bound).Intersect(grid_domain);
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            T d=(X-example.grid.Center(it.index)).Magnitude();
            ph.phi(it.index)=min(ph.phi(it.index),d-r);}}

    ARRAY<TV_INT> seed_indices;
    for(FACE_ITERATOR<TV> it(example.grid,example.ghost,GRID<TV>::INTERIOR_REGION);it.Valid();it.Next()){
        TV_INT a=it.First_Cell_Index(),b=it.Second_Cell_Index();
        if(ph.phi(a)<0){if(ph.phi(b)>=0) seed_indices.Append(b);}
        else if(ph.phi(b)<0) seed_indices.Append(a);}

    FAST_MARCHING_METHOD_UNIFORM<TV> fmm(*ph.levelset,example.ghost);
    fmm.Fast_Marching_Method(ph.phi,3*dx,&seed_indices);
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

    LINEAR_INTERPOLATION_MAC<TV,T> li(example.grid);
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
        [this,dt,&particles,use_flip,&li,&ph](int p,HELPER& h)
        {
            if(particles.store_B) particles.B(p)=h.B;
            if(!example.rk_particle_order) particles.X(p)+=dt*h.V;
            else // RK step particles with linear interpolation
                for(RUNGEKUTTA<TV> rk(particles.X(p),example.rk_particle_order,dt,0);rk.Valid();rk.Next())
                    particles.X(p)+=dt*li.Clamped_To_Array(ph.velocity,particles.X(p));
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
// Function Move_Mass_Momentum_Inside
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Move_Mass_Momentum_Inside(PHASE& ph) const
{
    for(FACE_ITERATOR<TV> fit(example.grid,2);fit.Valid();fit.Next()){
        TV loc=fit.Location();
        T mass=ph.mass(fit.Full_Index());
        if(!mass) continue;
        T momentum=ph.velocity(fit.Full_Index());
        T volume=0;
        if(example.use_particle_volumes)
            volume=ph.volume(fit.Full_Index());
        for(int i=0;i<example.collision_objects.m;i++){
            MPM_COLLISION_OBJECT<TV>* o=example.collision_objects(i);
            if(o->Phi(loc,example.time)>=0) continue; // outside collision object
            Add_Debug_Particle(loc,VECTOR<T,3>(0,0,1)); // DEBUG: the point mass will be moved out from
            IMPLICIT_OBJECT<TV>* io=o->Get_Implicit_Object(example.time);
            TV loc_bd=io->Closest_Point_On_Boundary(loc);
            TV loc_reflect=loc+2*(loc_bd-loc); 
            Add_Debug_Particle(loc_reflect,VECTOR<T,3>(1,0,0)); // DEBUG: reflection across boundary
            FACE_INDEX<TV::m> first_face_index=example.grid.Face(loc_reflect,fit.axis); // get the smallest face index
            TV first_face_loc=example.grid.Face(first_face_index);
            TV offset=(loc_reflect-first_face_loc)*example.grid.one_over_dX;
            VECTOR<T,1<<TV::m> weights;
            for(int j=0;j<1<<TV::m;j++){
                TV w=offset,loc=first_face_loc;
                for(int k=0;k<TV::m;k++)
                    if(j&(1<<k)){
                        w(k)=1-w(k);
                        loc(k)+=example.grid.dX(k);}
                weights(j)=w.Product();
                if(o->Phi(loc,example.time)>0){
                    Add_Debug_Particle(loc,VECTOR<T,3>(1,1,0));/* DEBUG: mass moved to*/}
                else{
                    weights(j)=0;
                    Add_Debug_Particle(loc,VECTOR<T,3>(0,1,1));/* DEBUG: mass won't be moved to*/}}
            // scale weights if not all nodes are inside
            weights/=weights.Sum();
            for(int j=0;j<1<<TV::m;j++){
                if(weights(j)){
                    FACE_INDEX<TV::m> face(first_face_index);
                    for(int k=0;k<TV::m;k++)
                        if(j&(1<<k))
                            face.index(k)++;
                    // move mass momentum, inside
                    ph.mass(face)+=weights(j)*mass;
                    ph.velocity(face)+=weights(j)*momentum;
                    if(example.use_particle_volumes) 
                        ph.volume(face)+=weights(j)*volume;}}
            PHYSBAM_DEBUG_WRITE_SUBSTEP("faces mass moved to",0,1);
            // clear mass,momentum outside
            ph.mass(fit.Full_Index())=0;
            ph.velocity(fit.Full_Index())=0;
            if(example.use_particle_volumes)
                ph.volume(fit.Full_Index())=0;
            break;}}
}
//#####################################################################
// Function Move_Mass_Momentum_Inside
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Move_Mass_Momentum_Inside_Nearest(PHASE& ph) const
{
    INTERPOLATED_COLOR_MAP<T> color_map;
    color_map.Initialize_Colors(0,ph.mass.array.Max(),false,true,false);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("purge",0,1);
    for(FACE_ITERATOR<TV> it(example.grid,2);it.Valid();it.Next())
        if(ph.mass(it.Full_Index()))
            Add_Debug_Particle(it.Location(),color_map(ph.mass(it.Full_Index())));
    PHYSBAM_DEBUG_WRITE_SUBSTEP("masses before ",0,1);

    for(FACE_ITERATOR<TV> fit(example.grid,2);fit.Valid();fit.Next()){
        TV loc=fit.Location();
        T mass=ph.mass(fit.Full_Index());
        if(!mass) continue;
        T momentum=ph.velocity(fit.Full_Index());
        T volume=0;
        if(example.use_particle_volumes)
            volume=ph.volume(fit.Full_Index());
        for(int i=0;i<example.collision_objects.m;i++){
            MPM_COLLISION_OBJECT<TV>* o=example.collision_objects(i);
            if(o->Phi(loc,example.time)>=0) continue; // outside collision object
            IMPLICIT_OBJECT<TV>* io=o->Get_Implicit_Object(example.time);
            TV loc_bd=io->Closest_Point_On_Boundary(loc);
            TV loc_reflect=loc+2*(loc_bd-loc); 
            T closest_distance=FLT_MAX;
            FACE_INDEX<TV::m> first_face_index=example.grid.Face(loc_reflect,fit.axis); // get the smallest face index
            FACE_INDEX<TV::m> closest_face(first_face_index);
            for(int j=0;j<1<<TV::m;j++){
                FACE_INDEX<TV::m> f(first_face_index);
                for(int k=0;k<TV::m;k++)
                    if(j&(1<<k))
                        f.index(k)++;
                T dist2=(example.grid.Face(f)-loc).Magnitude_Squared();
                if(ph.mass(f) && dist2<closest_distance){
                    closest_distance=dist2;
                    closest_face=f;}}
            Add_Debug_Object(VECTOR<TV,2>(loc,example.grid.Face(closest_face)),VECTOR<T,3>(1,1,0));
            ph.mass(closest_face)+=ph.mass(fit.Full_Index());
            ph.mass(fit.Full_Index())=0;
            ph.velocity(closest_face)+=ph.velocity(fit.Full_Index());
            ph.velocity(fit.Full_Index())=0;
            if(example.use_particle_volumes){
                ph.volume(closest_face)+=ph.volume(fit.Full_Index());
                ph.volume(fit.Full_Index())=0;}
            break;}}

    PHYSBAM_DEBUG_WRITE_SUBSTEP("purge",0,1);
    for(FACE_ITERATOR<TV> it(example.grid,2);it.Valid();it.Next())
        if(ph.mass(it.Full_Index()))
            Add_Debug_Particle(it.Location(),color_map(ph.mass(it.Full_Index())));
    PHYSBAM_DEBUG_WRITE_SUBSTEP("masses after",0,1);
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
    ARRAY<bool,FACE_INDEX<TV::m> > psi_N(example.grid,3,false);

#pragma omp parallel
    for(FACE_ITERATOR_THREADED<TV> it(example.grid,3);it.Valid();it.Next()){
        TV X=it.Location();
        bool N=false;
        for(int i=0;i<example.collision_objects.m;i++){
            MPM_COLLISION_OBJECT<TV>* o=example.collision_objects(i);
            if(o->Phi(X,example.time)<0){
                N=true;
                T v=o->Velocity(X,example.time)(it.axis);
                for(PHASE_ID p(0);p<example.phases.m;p++)
                    if(example.phases(p).mass(it.Full_Index()))
                        example.phases(p).velocity(it.Full_Index())=v;
                break;}}
        psi_N(it.Full_Index())=N;}

    example.projection_system.dc_present=false;
    ARRAY<int,TV_INT> cell_index(example.grid.Domain_Indices(1),true,-1);
    ARRAY<PHASE_ID,TV_INT> cell_phase(example.grid.Domain_Indices(1),false);
    int next_cell=0;
    for(CELL_ITERATOR<TV> it(example.grid);it.Valid();it.Next()){
        bool dirichlet=false,all_N=true;
        PHASE_ID phase(-1);
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> face(a,it.index);
            for(int s=0;s<2;s++){
                if(!psi_N(face)){
                    all_N=false;
                    bool has_mass=false;
                    for(PHASE_ID p(0);p<example.phases.m;p++)
                        if(example.phases(p).mass(face)){
                            phase=p;
                            has_mass=true;}
                    if(!has_mass){
                        dirichlet=true;
                        example.projection_system.dc_present=true;}}
                face.index(a)++;}}
        if(!all_N && !dirichlet){
            cell_index(it.index)=next_cell++;
            cell_phase(it.index)=phase;}}

    example.projection_system.A.Reset(next_cell);
    ARRAY<int> tmp0,tmp1;
#pragma omp parallel
    {
        SPARSE_MATRIX_THREADED_CONSTRUCTION<T> helper(example.projection_system.A,tmp0,tmp1);
        for(CELL_ITERATOR_THREADED<TV> it(example.grid,0);it.Valid();it.Next()){
            int center_index=cell_index(it.index);
            if(center_index<0) continue;
            PHASE_ID p=cell_phase(it.index);

            T diag=0;
            helper.Start_Row();
            for(int a=0;a<TV::m;a++){
                FACE_INDEX<TV::m> face(a,it.index);
                for(int s=0;s<2;s++){
                    TV_INT cell=it.index;
                    cell(a)+=2*s-1;
                    assert(face.Cell_Index(s)==cell);
                    if(!psi_N(face)){
                        T entry=sqr(example.grid.one_over_dX(a));
                        if(example.use_massless_particles)
                            entry/=example.phases(p).density;
                        else{
                            entry/=example.phases(p).mass(face);
                            if(example.use_particle_volumes)
                                entry*=example.phases(p).volume(face);}
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
    example.projection_system.phases.Remove_All();
    example.projection_system.gradient.Reset(next_cell);
    APPEND_HOLDER<T> mass_h(example.projection_system.mass);
    APPEND_HOLDER<FACE_INDEX<TV::m> > faces_h(example.projection_system.faces);
    APPEND_HOLDER<PHASE_ID> phases_h(example.projection_system.phases);
    ARRAY<int> tmp2,tmp3;
#pragma omp parallel
    {
#pragma omp single
        {
            mass_h.Init();
            faces_h.Init();
            phases_h.Init();
        }
#pragma omp barrier
        ARRAY<T>& mass_t=mass_h.Array();
        ARRAY<FACE_INDEX<TV::m> >& faces_t=faces_h.Array();
        ARRAY<PHASE_ID>& phases_t=phases_h.Array();
        SPARSE_MATRIX_THREADED_CONSTRUCTION<T> G_helper(example.projection_system.gradient,tmp2,tmp3);
        for(FACE_ITERATOR_THREADED<TV> it(example.grid);it.Valid();it.Next()){
            if(psi_N(it.Full_Index())) continue;
            for(PHASE_ID p(0);p<example.phases.m;p++){
                T mass=0;
                if(example.use_massless_particles) mass=example.phases(p).density;
                else{
                    mass=example.phases(p).mass(it.Full_Index());
                    if(!mass) continue;
                    if(example.use_particle_volumes)
                        mass/=example.phases(p).volume(it.Full_Index());}
                int c0=cell_index(it.First_Cell_Index());
                int c1=cell_index(it.Second_Cell_Index());
                if(c0<0 && c1<0) continue;
                G_helper.Start_Row(); // cannot start a row if no elements in it
                if(c0>=0) G_helper.Add_Entry(c0,-example.grid.one_over_dX(it.axis));
                if(c1>=0) G_helper.Add_Entry(c1,example.grid.one_over_dX(it.axis));
                faces_t.Append(it.Full_Index());
                phases_t.Append(p);
                mass_t.Append(mass);}}
        G_helper.Finish();
    }
    mass_h.Combine();
    faces_h.Combine();
    phases_h.Combine();
}
//#####################################################################
// Function Pressure_Projection
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Pressure_Projection()
{
    static int solve_id=-1;
    TIMER_SCOPE_FUNC;
    
    Compute_Poisson_Matrix();
    example.sol.v.Resize(example.projection_system.A.m);
    example.rhs.v.Resize(example.projection_system.A.m);

    ARRAY<T> tmp(example.projection_system.gradient.m);
    for(int i=0;i<tmp.m;i++)
        tmp(i)=example.phases(example.projection_system.phases(i)).velocity(example.projection_system.faces(i));
    example.projection_system.gradient.Transpose_Times(tmp,example.rhs.v);

    if(!example.projection_system.dc_present)
        example.projection_system.Compute_Ones_Nullspace();
    
    solve_id++;
    if(example.test_system){
        example.projection_system.Test_System(example.sol);
        example.projection_system.Test();}
    if(example.print_matrix){
        LOG::cout<<"solve id: "<<solve_id<<std::endl;
        KRYLOV_SOLVER<T>::Ensure_Size(example.av,example.sol,2);
        const MPM_PROJECTION_SYSTEM<TV>& system=example.projection_system;
        OCTAVE_OUTPUT<T>(LOG::sprintf("M-%i.txt",solve_id).c_str()).Write("M",system,*example.av(0),*example.av(1));
        OCTAVE_OUTPUT<T>(LOG::sprintf("P-%i.txt",solve_id).c_str()).Write_Projection("P",system,*example.av(0));
        OCTAVE_OUTPUT<T>(LOG::sprintf("b-%i.txt",solve_id).c_str()).Write("b",example.rhs);}

    CONJUGATE_GRADIENT<T> cg;
    cg.finish_before_indefiniteness=true;
    cg.relative_tolerance=true;
    bool converged=cg.Solve(example.projection_system,example.sol,example.rhs,
        example.av,example.solver_tolerance,0,example.solver_iterations);
    if(!converged) LOG::printf("SOLVER DID NOT CONVERGE.\n");

    if(example.print_matrix)
        OCTAVE_OUTPUT<T>(LOG::sprintf("x-%i.txt",solve_id).c_str()).Write("x",example.sol);

    example.projection_system.gradient.Times(example.sol.v,tmp);
    for(int i=0;i<tmp.m;i++)
        example.phases(example.projection_system.phases(i)).velocity(example.projection_system.faces(i))-=tmp(i)/example.projection_system.mass(i);
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
Max_Particle_Speed(const PHASE& ph) const
{
    T v2=0;
#pragma omp parallel for reduction(max:v2)
    for(int k=0;k<ph.simulated_particles.m;k++){
        int p=ph.simulated_particles(k);
        v2=max(v2,example.particles.V(p).Magnitude_Squared());}
    return sqrt(v2);
}
//#####################################################################
// Function Max_Particle_Speed
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_DRIVER<TV>::
Max_Particle_Speed() const
{
    TIMER_SCOPE_FUNC;
    T v=0;
    for(PHASE_ID i(0);i<example.phases.m;i++)
        v=max(v,Max_Particle_Speed(example.phases(i)));
    return v;
}
//#####################################################################
// Function Grid_V_Upper_Bound
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_DRIVER<TV>::
Grid_V_Upper_Bound(const PHASE& ph) const
{
    T result=0;
    T xi=(T)6*sqrt((T)TV::m)*example.grid.one_over_dX.Min();
#pragma omp parallel for reduction(max:result)
    for(int k=0;k<ph.simulated_particles.m;k++){
        int p=ph.simulated_particles(k);
        T v=example.particles.V(p).Magnitude();
        if(example.particles.store_B) v+=example.particles.B(p).Frobenius_Norm()*xi;
        result=max(result,v);}
    return result;
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
    
    for(PHASE_ID i(0);i<example.phases.m;i++)
        result=max(result,Grid_V_Upper_Bound(example.phases(i)));
    return result;
}
//#####################################################################
// Function Update_Simulated_Particles
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Update_Simulated_Particles()
{
    TIMER_SCOPE_FUNC;

    if(!example.particles.store_phase){
        // If particles do not have phase attribute,
        // put them in the default(0) phase.
        PHASE& ph=example.phases(PHASE_ID(0));
        ph.simulated_particles.Remove_All();
        
        for(int p=0;p<example.particles.number;p++){
            if(example.particles.valid(p))
                ph.simulated_particles.Append(p);
        }
    }
    else{
        // Otherwise dispatch particles to their corresponding phase
        // based on their phase attribute value.
        for(PHASE_ID i(0);i<example.phases.m;i++)
            example.phases(i).simulated_particles.Remove_All();
        
        for(int p=0;p<example.particles.number;p++){
            int particle_phase=example.particles.phase(p);
            PHASE& ph=example.phases(PHASE_ID(particle_phase));
            if(example.particles.valid(p))
                ph.simulated_particles.Append(p);
        }
    }
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
