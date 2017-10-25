//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/FINE_TIMER.h>
#include <Core/Log/SCOPE.h>
#include <Core/Matrices/SPARSE_MATRIX_THREADED_CONSTRUCTION.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <Tools/Parallel_Computation/APPEND_HOLDER.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/CELL_ITERATOR_THREADED.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR_THREADED.h>
#include <Grid_Tools/Grids/FACE_RANGE_ITERATOR.h>
#include <Grid_PDE/Boundaries/BOUNDARY_MAC_GRID_PERIODIC.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_MAC.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Grid_PDE/Interpolation/WENO_INTERPOLATION.h>
#include <Geometry/Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MAC_DRIVER.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MAC_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_FACE_ITERATOR.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <Hybrid_Methods/Projection/MPM_PROJECTION_SYSTEM.h>
#include <Hybrid_Methods/Projection/MPM_PROJECTION_VECTOR.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_MAC_DRIVER<TV>::
MPM_MAC_DRIVER(MPM_MAC_EXAMPLE<TV>& example)
    :example(example)
{
    DEBUG_SUBSTEPS::write_substeps_level=example.write_substeps_level;
    DEBUG_SUBSTEPS::writer=[=](const std::string& title){Write_Substep(title);};
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_MAC_DRIVER<TV>::
~MPM_MAC_DRIVER()
{
    DEBUG_SUBSTEPS::writer=0;
}
//#####################################################################
// Execute_Main_Program
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Execute_Main_Program()
{
    Step([=](){Initialize();},"initialize");
    Simulate_To_Frame(example.last_frame);
    for(const auto& it:example.time_step_callbacks)
        if(!it.data.x)
            PHYSBAM_FATAL_ERROR("Unrecognized time step callback: "+it.key);
    FINE_TIMER::Dump_Timing_Info();
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Initialize()
{
    TIMER_SCOPE_FUNC;
    LOG::cout<<std::setprecision(16);
    DEBUG_SUBSTEPS::write_substeps_level=example.substeps_delay_frame<0?example.write_substeps_level:-1;

    // setup time
    output_number=current_frame=example.restart;

    example.Initialize();
    for(int i=0;i<TV::m;i++){
        bool a=example.bc_type(i*2)==example.BC_PERIODIC;
        bool b=example.bc_type(i*2+1)==example.BC_PERIODIC;
        PHYSBAM_ASSERT(a==b);
        example.periodic_boundary.is_periodic(i)=a;}

    example.radius_sphere=(T)0.36*example.grid.dX.Max();
    example.radius_escape=TV::m==2?(T)1.5*example.grid.dX.Max():(T)1.1*example.grid.dX.Max();

    // must have level sets in order to use multiphase projection
    if(example.use_multiphase_projection){
        PHYSBAM_ASSERT(example.use_phi);}

    PHYSBAM_ASSERT(example.grid.Is_MAC_Grid());
    if(example.restart)
        example.Read_Output_Files(example.restart);

    example.particles.Store_B(example.use_affine);
    if(example.particles.store_B)
        for(int i=0;i<TV::m;i++)
            example.Dp_inv(i).Resize(example.particles.X.m);

    for(PHASE_ID i(0);i<example.phases.m;i++){
        PHASE& ph=example.phases(i);
        ph.Initialize(example.grid,example.weights,example.ghost,example.threads);
        ph.phi.Resize(example.grid.Domain_Indices(example.ghost));
        if(example.use_phi){
            ph.levelset=new LEVELSET<TV>(example.grid,ph.phi,example.ghost);
            ph.levelset->boundary=&example.periodic_boundary;}}

    RANGE<TV_INT> range(example.grid.Cell_Indices(example.ghost));
    example.location.Resize(example.grid,example.ghost);
    for(FACE_ITERATOR<TV> it(example.grid,example.ghost);it.Valid();it.Next())
        example.location(it.Full_Index())=it.Location();

    Update_Simulated_Particles();
    // Need grid velocities for initial advection step
    if(example.flip || example.rk_particle_order){
        Update_Particle_Weights();
        Prepare_Scatter();
        Particle_To_Grid();}

    if(example.extrap_type=='p')
        example.force.Resize(example.grid.Domain_Indices(example.ghost));

    if(!example.restart) Write_Output_Files(0);
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Advance_One_Time_Step()
{
    TIMER_SCOPE_FUNC;
    Step([=](){Move_Particles();},"move-particles",true);
    Step([=](){Update_Simulated_Particles();},"simulated-particles",false);
    Step([=](){Update_Particle_Weights();},"update-weights",false);
    Step([=](){Prepare_Scatter();},"prepare-scatter",false);
    Step([=](){Particle_To_Grid();},"p2g");
    Step([=](){Build_Level_Sets();},"build-level-sets",false,example.use_phi);
    Step([=](){Bump_Particles();},"bump-particles",true,example.use_bump);
    Step([=](){Reseeding();},"reseeding",true,example.use_reseeding);
    Step([=](){Apply_Forces();},"forces");
    Step([=](){Compute_Boundary_Conditions();},"compute boundary conditions",false);
    Step([=](){Pressure_Projection();},"projection");
    Step([=](){Apply_Viscosity();},"viscosity",true,example.use_viscosity);
    Step([=](){
        if(example.extrap_type=='p')
            Reflect_Boundary(&MPM_MAC_DRIVER::Reflect_Boundary_Velocity_Copy_Only,RF::ghost|RF::delay_corners);
        else Extrapolate_Velocity(!example.flip);},"velocity-extrapolation",true);
    Step([=](){Grid_To_Particle();},"g2p");
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
        for(int i=0;i<example.begin_frame.m;i++)
            example.begin_frame(i)(current_frame);
        if(example.substeps_delay_frame==current_frame)
            DEBUG_SUBSTEPS::write_substeps_level=example.write_substeps_level;
        T time_at_frame=example.time+example.frame_dt;
        bool done=false;
        for(int substep=0;!done;substep++){
            LOG::SCOPE scope("SUBSTEP","substep %d",substep+1);
            example.dt=Compute_Dt();
            example.dt=clamp(example.dt,example.min_dt,example.max_dt);
            T next_time=example.time+example.dt;
            if(next_time+(T).0001*example.dt>time_at_frame){ // Account for roundoff errors
                next_time=time_at_frame;
                done=true;}
            else if(next_time+example.dt>time_at_frame) next_time=(example.time+time_at_frame)/2;
            example.dt=next_time-example.time;
            LOG::cout<<"substep dt: "<<example.dt<<std::endl;

            Step([=](){Advance_One_Time_Step();},"time-step");
            PHYSBAM_DEBUG_WRITE_SUBSTEP("end substep %i",0,substep);
            LOG::cout<<"max velocity: "<<Max_Particle_Speed()<<"  bound: "<<Grid_V_Upper_Bound()<<std::endl;
            example.time=next_time;}
        for(int i=0;i<example.end_frame.m;i++)
            example.end_frame(i)(current_frame);
        Write_Output_Files(++output_number);}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Write_Substep(const std::string& title)
{
    TIMER_SCOPE_FUNC;
    example.frame_title=title;
    LOG::printf("Writing substep [%s]: output_number=%i, time=%g, frame=%i\n",
        example.frame_title,output_number+1,example.time,current_frame);
    Write_Output_Files(++output_number);
    example.frame_title="";
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
    for(int i=0;i<TV::m;++i){
        example.weights(i)->Update(example.particles.X);
        if(example.particles.store_B)
            example.weights(i)->Dp_Inverse(example.particles.X,example.Dp_inv(i));}
}
//#####################################################################
// Function Extrapolate_Boundary
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Extrapolate_Boundary(PHASE& ph) const
{
    int d=1;
    auto extrapolate=[this,&ph,d](auto it){
        if(example.bc_type(it.side)==example.BC_PERIODIC) return;
        int side_axis=it.side/2;
        int sign_out=it.side%2?1:-1;
        auto further=[=](FACE_INDEX<TV::m> f){f.index(side_axis)-=sign_out;return f;};
        FACE_INDEX<TV::m> f0=it.face;
        for(int j=d;j>0;--j) f0=further(f0);
        FACE_INDEX<TV::m> f1=further(f0);
        ph.velocity(it.face)=(d+1)*ph.velocity(f0)-d*ph.velocity(f1);};

    for(int n=example.weights[0]->stencil_width-2;n>0;--n) // assume all dimensions use the same interpolation order
        for(FACE_RANGE_ITERATOR<TV::m> it(example.grid.Domain_Indices(),1-n,-n,
                RF::ghost|RF::skip_inner|RF::delay_corners);it.Valid();it.Next())
            extrapolate(it);
}
//#####################################################################
// Function Particle_To_Grid
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Particle_To_Grid(PHASE_ID pid) const
{
    TIMER_SCOPE_FUNC;
    const MPM_PARTICLES<TV>& particles=example.particles;
    PHASE& ph=example.phases(pid);

#pragma omp parallel for
    for(int i=0;i<ph.mass.array.m;i++){
        ph.mass.array(i)=0;
        ph.volume.array(i)=0;
        ph.velocity.array(i)=0;}

    ph.gather_scatter->template Scatter<int>(true,
        [this,&particles,&ph](int p,const PARTICLE_GRID_FACE_ITERATOR<TV>& it,int data)
        {
            T w=it.Weight();
            FACE_INDEX<TV::m> index=it.Index();
            ph.mass(index)+=w*particles.mass(p);
            ph.volume(index)+=w*particles.volume(p);
            T V=w*particles.V(p)(index.axis);
            if(example.use_affine){
                TV Y=example.grid.Face(index)-particles.X(p);
                if(!example.lag_Dp) Y=example.Dp_inv(index.axis)(p)*Y;
                V+=particles.B(p).Row(index.axis).Dot(Y)*w;}
            ph.velocity(index)+=particles.mass(p)*V;
        });
    Fix_Periodic_Accum(ph.mass);
    Fix_Periodic_Accum(ph.velocity);
    Fix_Periodic_Accum(ph.volume);

    if(example.flip && example.extrap_type=='p'){
        PHYSBAM_DEBUG_WRITE_SUBSTEP("before reflect",1);
        Reflect_Boundary(&MPM_MAC_DRIVER::Reflect_Boundary_Mass_Momentum,RF::ghost|RF::duplicate_corners);}

    // TODO: move mass, momentum, volume inside
    if(example.move_mass_inside)
        Move_Mass_Momentum_Inside(ph);
    else if(example.move_mass_inside_nearest)
        Move_Mass_Momentum_Inside_Nearest(ph);

    ph.valid_indices.Remove_All();
    ph.valid_flat_indices.Remove_All();

#pragma omp parallel
    for(FACE_ITERATOR_THREADED<TV> it(example.grid,example.ghost);it.Valid();it.Next()){
        int i=ph.mass.Standard_Index(it.Full_Index());
        if(ph.mass.array(i)){
            ph.velocity.array(i)/=ph.mass.array(i);}
        else ph.velocity.array(i)=0;}

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
        for(FACE_ITERATOR_THREADED<TV> it(example.grid);it.Valid();it.Next()){
            int i=ph.mass.Standard_Index(it.Full_Index());
            if(ph.mass.array(i)){
                flat_t.Append(i);
                indices_t.Append(it.Full_Index());}}
    }
    flat_h.Combine();
    indices_h.Combine();
    Fix_Periodic(ph.mass);
    Fix_Periodic(ph.velocity);
    Fix_Periodic(ph.volume);
    if(example.flip){
        if(example.extrap_type!='p'){
            Extrapolate_Boundary(ph);
            Extrapolate_Velocity(pid,false);}
        ph.velocity_save=ph.velocity;}
}
//#####################################################################
// Function Build_Level_Sets
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Build_Level_Sets()
{
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
        T influence_bound=example.radius_sphere+dx*(T)1.1;
        TV X=example.particles.X(p);
        RANGE<TV> bound(X-influence_bound,X+influence_bound);
        RANGE<TV_INT> range=example.grid.Clamp_To_Cell(bound,example.ghost).Intersect(grid_domain);
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            T d=(X-example.grid.Center(it.index)).Magnitude();
            ph.phi(it.index)=min(ph.phi(it.index),d-example.radius_sphere);}}

    ARRAY<TV_INT> seed_indices;
    for(FACE_ITERATOR<TV> it(example.grid,example.ghost,GRID<TV>::INTERIOR_REGION);it.Valid();it.Next()){
        TV_INT a=it.First_Cell_Index(),b=it.Second_Cell_Index();
        if(ph.phi(a)<0){if(ph.phi(b)>=0) seed_indices.Append(b);}
        else if(ph.phi(b)<0) seed_indices.Append(a);}
    Prune_Duplicates(seed_indices);

    FAST_MARCHING_METHOD_UNIFORM<TV> fmm;
    fmm.seed_indices=&seed_indices;
    fmm.correct_interface_phi=false;
    fmm.process_sign=1;
    fmm.Fast_Marching_Method(example.grid,example.ghost,ph.phi,3*dx);

    T shift=10*dx;
    for(int i=0;i<seed_indices.m;i++)
        ph.phi(seed_indices(i))-=shift;

    fmm.process_sign=-1;
    fmm.Fast_Marching_Method(example.grid,example.ghost,ph.phi,shift+3*dx);

    for(CELL_ITERATOR<TV> it(example.grid,example.ghost);it.Valid();it.Next())
        if(ph.phi(it.index)<=0) ph.phi(it.index)+=shift;
}
//#####################################################################
// Function Bump_Particles
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Bump_Particles()
{
    T weno_eps=(T)1e-6;
    T dx=example.grid.dX.Max();
    T one_over_dx=(T)1/dx;
    TV halfdx=(T).5*example.grid.dX;
    for(PHASE_ID i(0);i<example.phases.m;i++) {
        PHASE& ph=example.phases(i);

        ARRAY<T,TV_INT> curvature(example.grid.Cell_Indices(example.ghost-1));
        ARRAY<TV,TV_INT> gradient(example.grid.Cell_Indices(example.ghost-1));
        ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT> Hessian(example.grid.Cell_Indices(example.ghost-1));
        for(CELL_ITERATOR<TV> it(example.grid,example.ghost-1);it.Valid();it.Next()){
            curvature(it.index)=ph.levelset->Compute_Curvature(ph.phi,it.index);
            gradient(it.index)=ph.levelset->Gradient(ph.phi,it.index);
            Hessian(it.index)=ph.levelset->Hessian(ph.phi,it.index);}

        for(int j=0;j<ph.simulated_particles.m;j++){
            int p=ph.simulated_particles(j);
            TV X=example.particles.X(p);
            TV_INT index=example.grid.Cell(X-halfdx);
            TV c=(X-example.grid.Center(index))*example.grid.One_Over_DX();
            T phi=WENO_Interpolation(c,ph.phi,index,weno_eps);
            if(phi>-example.radius_sphere && phi<example.radius_escape){
                T cell_k[1<<TV::m];
                for(int l=0;l<(1<<TV::m);l++){
                    TV_INT cell(index);
                    for(int d=0;d<TV::m;d++)
                        if(l&(1<<d))
                            cell(d)++;
                    cell_k[l]=curvature(cell);}
                T absk=abs(LINEAR_INTERPOLATION<T,T>::Linear(cell_k,c));

                T target=(T)0;
                if(absk>(T)0.5*one_over_dx)
                    continue;
                if(absk<(T)0.25*one_over_dx)
                    target=example.radius_sphere;

                TV nearest=Nearest_Point_On_Surface(X,ph,gradient,Hessian);
                T norm_grad=WENO_Interpolation(c,gradient,index,weno_eps).Magnitude();
                T d=phi/norm_grad;

                TV nq=nearest-X;
                nq.Normalize();
                if(d>=0)
                    example.particles.X(p)=nearest+target*nq;
                else if(d>-target)
                    example.particles.X(p)=nearest-target*nq;}}}
}
//#####################################################################
// Function Nearest_Point_On_Surface
//#####################################################################
template<class TV> TV MPM_MAC_DRIVER<TV>::
Nearest_Point_On_Surface(const TV& p,const PHASE& ph,
    const ARRAY<TV,TV_INT>& gradient,
    const ARRAY<SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>,TV_INT>& Hessian) const
{
    typedef VECTOR<T,TV::m+1> ETV;

    int max_iterations=20;
    T phi_tolerance=(T)0.01*example.grid.dX.Min(),dir_tolerance=(T)0.02;
    auto Compute=[this,&gradient,&Hessian,&ph](const TV& X,TV& n,SYMMETRIC_MATRIX<T,TV::m>* h)
    {
        T weno_eps=(T)1e-6;
        TV_INT index=example.grid.Cell(X-(T)0.5*example.grid.dX);
        TV c=(X-example.grid.Center(index))*example.grid.one_over_dX;
        n=WENO_Interpolation(c,gradient,index,weno_eps).Normalized();
        if(h) *h=WENO_Interpolation(c,Hessian,index,weno_eps);
        return WENO_Interpolation(c,ph.phi,index,weno_eps);
    };

    TV n;
    T phi=Compute(p,n,0);
    TV x=p-phi*n;
    T lambda=-2*phi;
    for(int i=0;i<max_iterations;i++){
        SYMMETRIC_MATRIX<T,TV::m> h;
        T phi=Compute(x,n,&h);
        TV xp=x-p;
        if(abs(phi)<=phi_tolerance){
            T dist2=xp.Magnitude_Squared();
            if(xp.Projected_Orthogonal_To_Unit_Direction(n).Magnitude_Squared()<=sqr(dir_tolerance)*dist2) break;
            if(dist2<=sqr(phi_tolerance)){
                x=p-sign(xp.Dot(n))*phi*n;
                break;}}
        ETV g=(xp*2-lambda*n).Append(-phi);
        MATRIX<T,TV::m+1,TV::m+1> H;
        H.Set_Submatrix(0,0,(T)2-lambda*h);
        ETV en=(-n).Append(0);
        H.Set_Column(TV::m,en);
        H.Set_Row(TV::m,en);
        ETV dx=H.PLU_Solve(-g);
        x+=dx.Remove_Index(TV::m);
        lambda+=dx(TV::m);}
    return x;
}
//#####################################################################
// Function Reseeding
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Reseeding()
{
    LINEAR_INTERPOLATION_MAC<TV,T> li(example.grid);
    int target_ppc=1<<TV::m;
    ARRAY<ARRAY<int>,TV_INT> particle_counts(example.grid.Domain_Indices());
    for(int i=0;i<example.particles.X.m;i++){
        if(!example.particles.valid(i)) continue;
        TV X=example.particles.X(i);
        TV_INT index=example.grid.Cell(X);
        particle_counts(index).Append(i);}
    ARRAY<ARRAY<int>,PHASE_ID> particles_by_phase(example.phases.m);
    for(CELL_ITERATOR<TV> it(example.grid);it.Valid();it.Next()){
        ARRAY<int>& a=particle_counts(it.index);
        int target=target_ppc;
        PAIR<PHASE_ID,T> poc=Phase_Of_Cell(it.index);
        if(abs(poc.y)<example.grid.dX.Max()) target*=4;
        if(a.m<target){
            RANGE<TV> range=example.grid.Cell_Domain(it.index);
            for(int i=0;i<target-a.m;i++){
                TV X;
                example.random.Fill_Uniform(X,range);
                PHASE_ID pid;
                T phi=Phase_And_Phi(X,pid);
                if(abs(phi)>=example.radius_sphere){
                    int p=example.particles.Add_Element_From_Deletion_List();
                    example.particles.valid(p)=true;
                    example.particles.mass(p)=1;
                    example.particles.phase(p)=pid;
                    example.particles.V(p)=li.Clamped_To_Array(example.phases(pid).velocity,X);
                    example.particles.X(p)=X;}}}
        if(a.m>2*target){
            for(int i=0;i<a.m;i++)
                particles_by_phase(example.particles.phase(a(i))).Append(a(i));
            for(PHASE_ID i(0);i<particles_by_phase.m;i++){
                if(particles_by_phase(i).m>2*target){
                    example.random.Random_Shuffle(particles_by_phase(i));
                    int p=particles_by_phase(i).Pop();
                    Invalidate_Particle(p);}
                particles_by_phase(i).Remove_All();}}}
    Update_Simulated_Particles();
    Update_Particle_Weights();
    Prepare_Scatter();
}
//#####################################################################
// Function Phase_And_Phi
//#####################################################################
template<class TV> auto MPM_MAC_DRIVER<TV>::
Phase_And_Phi(const TV& X,PHASE_ID& phase) const->T
{
    T min_phi=FLT_MAX;
    for(PHASE_ID i(0);i<example.phases.m;i++){
        T phi=example.phases(i).levelset->Phi(X);
        if(phi<min_phi){
            min_phi=phi;
            phase=i;}}
    return min_phi;
}
//#####################################################################
// Function Particle_To_Grid
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Particle_To_Grid()
{
    for(PHASE_ID i(0);i<example.phases.m;i++)
        Particle_To_Grid(i);
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
            if(use_flip) particles.V(p)=(1-example.flip)*h.V+example.flip*(particles.V(p)+h.flip_V);
            else particles.V(p)=h.V;
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
            FACE_INDEX<TV::m> first_face_index=example.grid.Face(loc_reflect,fit.face.axis); // get the smallest face index
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
            PHYSBAM_DEBUG_WRITE_SUBSTEP("faces mass moved to",1);
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
    PHYSBAM_DEBUG_WRITE_SUBSTEP("purge",1);
    for(FACE_ITERATOR<TV> it(example.grid,2);it.Valid();it.Next())
        if(ph.mass(it.Full_Index()))
            Add_Debug_Particle(it.Location(),color_map(ph.mass(it.Full_Index())));
    PHYSBAM_DEBUG_WRITE_SUBSTEP("masses before ",1);

    for(FACE_ITERATOR<TV> fit(example.grid,2);fit.Valid();fit.Next()){
        TV loc=fit.Location();
        T mass=ph.mass(fit.Full_Index());
        if(!mass) continue;
        for(int i=0;i<example.collision_objects.m;i++){
            MPM_COLLISION_OBJECT<TV>* o=example.collision_objects(i);
            if(o->Phi(loc,example.time)>=0) continue; // outside collision object
            IMPLICIT_OBJECT<TV>* io=o->Get_Implicit_Object(example.time);
            TV loc_bd=io->Closest_Point_On_Boundary(loc);
            TV loc_reflect=loc+2*(loc_bd-loc);
            T closest_distance=FLT_MAX;
            FACE_INDEX<TV::m> first_face_index=example.grid.Face(loc_reflect,fit.face.axis); // get the smallest face index
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

    PHYSBAM_DEBUG_WRITE_SUBSTEP("purge",1);
    for(FACE_ITERATOR<TV> it(example.grid,2);it.Valid();it.Next())
        if(ph.mass(it.Full_Index()))
            Add_Debug_Particle(it.Location(),color_map(ph.mass(it.Full_Index())));
    PHYSBAM_DEBUG_WRITE_SUBSTEP("masses after",1);
}
//#####################################################################
// Function Face_Fraction
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_DRIVER<TV>::
Face_Fraction(const FACE_INDEX<TV::m>& face_index,const ARRAY<T,TV_INT>& phi) const
{
    T phi0=phi(face_index.First_Cell_Index());
    T phi1=phi(face_index.Second_Cell_Index());
    T dx=example.grid.dX(face_index.axis);
    T m=sqr(dx)-sqr(phi1-phi0);
    T s=phi0+phi1;
    if(s*s<m) return 0.5-s/(2*sqrt(m));
    else if(s>=0) return 0;
    else return 1;
}
//#####################################################################
// Function Face_Fraction
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Face_Fraction(const FACE_INDEX<TV::m>& face_index,ARRAY<PAIR<PHASE_ID,T> >& face_fractions) const
{
    face_fractions.Remove_All();
    if(!example.use_multiphase_projection){
        for(PHASE_ID p=PHASE_ID(0);p<example.phases.m;p++){
            if(example.phases(p).mass(face_index)){
                face_fractions={{p,1}};
                return;}}
        return;}

    LINEAR_INTERPOLATION_UNIFORM<TV,T> li;
    GRID<TV> center_grid=example.grid.Get_Center_Grid();

    VECTOR<T,GRID<TV>::number_of_nodes_per_face> phis;
    VECTOR<PHASE_ID,GRID<TV>::number_of_nodes_per_face> phases;
    for(int i=0;i<phis.m;i++){
        TV X=example.grid.Node(example.grid.Face_Node_Index(face_index.axis,face_index.index,i));
        for(PHASE_ID p(0);p<example.phases.m;p++){
            T x=li.Clamped_To_Array(center_grid,example.phases(p).phi,X);
            if(x<=0){
                phis(i)=-x;
                phases(i)=p;}}}
    face_fractions.Append({phases(0),1});
    for(int i=1;i<phases.m;i++)
        if(phases(i)!=phases(0)){
            face_fractions.Append({phases(i),0});
            for(int j=i+1;j<phases.m;j++)
                PHYSBAM_ASSERT(phases(j)==phases(0) || phases(j)==phases(i));
            break;}

    if(face_fractions.m==1) return;
    T fraction=Face_Fraction(face_index,example.phases(face_fractions(0).x).phi);
    face_fractions(0).y=fraction;
    face_fractions(1).y=1-fraction;
}
//#####################################################################
// Function Phase_Of_Cell
//#####################################################################
template<class TV> PAIR<PHASE_ID,typename TV::SCALAR> MPM_MAC_DRIVER<TV>::
Phase_Of_Cell(const TV_INT& cell_index) const
{
    for(PHASE_ID p(0);p<example.phases.m;p++){
        T x=example.phases(p).phi(cell_index);
        if(x<=0) return {p,x};}
    PHYSBAM_FATAL_ERROR("No phase!");
}
//#####################################################################
// Function Density
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_DRIVER<TV>::
Density(const FACE_INDEX<TV::m>& face_index) const
{
    if(!example.use_multiphase_projection){
        for(PHASE_ID p=PHASE_ID(0);p<example.phases.m;p++)
            if(T m=example.phases(p).mass(face_index)){
                if(example.use_constant_density)
                    return example.phases(p).density;
                if(example.use_particle_volumes)
                    m/=example.phases(p).volume(face_index);
                return m;}
        return 0;}

    auto p0=Phase_Of_Cell(face_index.First_Cell_Index());
    auto p1=Phase_Of_Cell(face_index.Second_Cell_Index());
    if(p0.x==p1.x) return example.phases(p0.x).density;
    T theta=p0.y/(p0.y+p1.y);
    return example.phases(p0.x).density*theta+example.phases(p1.x).density*(1-theta);
}
//#####################################################################
// Function Neumann_Boundary_Condition
//#####################################################################
template<class TV> bool MPM_MAC_DRIVER<TV>::
Neumann_Boundary_Condition(const FACE_INDEX<TV::m>& face,ARRAY<T,PHASE_ID>& bc) const
{
    VECTOR<RANGE<TV_INT>,TV::m> domains;
    for(int i=0;i<TV::m;i++){
        domains(i)=example.grid.Domain_Indices();
        domains(i).min_corner(i)++;}
    for(int i=0;i<TV::m;i++){
        bool wall=example.bc_type(2*i)==example.BC_SLIP||example.bc_type(2*i)==example.BC_NOSLIP;
        if(face.index(i)<domains(face.axis).min_corner(i) && wall){
            for(PHASE_ID p(0);p<example.phases.m;p++)
                if(example.phases(p).mass(face)){
                    bc(p)=0;
                    if(example.bc_velocity(2*i)){
                        TV X=example.grid.domain.Clamp(example.grid.Face(face));
                        bc(p)=example.bc_velocity(2*i)(X,face.axis,p,example.time);}}
            return true;}
        wall=example.bc_type(2*i+1)==example.BC_SLIP||example.bc_type(2*i+1)==example.BC_NOSLIP;
        if(face.index(i)>=domains(face.axis).max_corner(i) && wall){
            for(PHASE_ID p(0);p<example.phases.m;p++)
                if(example.phases(p).mass(face)){
                    bc(p)=0;
                    if(example.bc_velocity(2*i+1)){
                        TV X=example.grid.domain.Clamp(example.grid.Face(face));
                        bc(p)=example.bc_velocity(2*i+1)(X,face.axis,p,example.time);}}
            return true;}}
    TV X=example.grid.Face(face);
    for(int i=0;i<example.collision_objects.m;i++){
        MPM_COLLISION_OBJECT<TV>* o=example.collision_objects(i);
        if(o->Phi(X,example.time)<0){
            bc.Fill(o->Velocity(X,example.time)(face.axis));
            return true;}}
    return false;
}
//#####################################################################
// Function Compute_Boundary_Conditions
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Compute_Boundary_Conditions()
{
    example.psi_N.Resize(example.grid,example.ghost,false);
#pragma omp parallel
    for(FACE_ITERATOR_THREADED<TV> it(example.grid,example.ghost);it.Valid();it.Next()){
        ARRAY<T,PHASE_ID> u_bc(example.phases.m);
        bool N=Neumann_Boundary_Condition(it.Full_Index(),u_bc);
        example.psi_N(it.Full_Index())=N;
        if(!N) continue;
        for(PHASE_ID p(0);p<example.phases.m;p++)
            if(example.phases(p).mass(it.Full_Index()))
                example.phases(p).velocity(it.Full_Index())=u_bc(p);}
    Fix_Periodic(example.psi_N);
    example.dof=Allocate_Projection_System_Variable();
}
//#####################################################################
// Function Allocate_Projection_System_Variable
//#####################################################################
template<class TV> int MPM_MAC_DRIVER<TV>::
Allocate_Projection_System_Variable()
{
    int ghost=2;
    example.cell_index.Resize(example.grid.Domain_Indices(ghost),false,false);
    for(int s=0;s<2*TV::m;s++){
        int value=pressure_uninit;
        typename MPM_MAC_EXAMPLE<TV>::BC_TYPE bc_type=example.bc_type(s);
        if(bc_type==example.BC_SLIP || bc_type==example.BC_NOSLIP) value=pressure_N;
        else if(example.bc_type(s)==example.BC_FREE) value=pressure_D;
        for(CELL_ITERATOR<TV> it(example.grid,ghost,GRID<TV>::GHOST_REGION,s);it.Valid();it.Next())
            example.cell_index(it.index)=value;}

    int nvar=0;
    for(CELL_ITERATOR<TV> it(example.grid);it.Valid();it.Next()){
        bool dirichlet=false,all_N=true;
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> face(a,it.index);
            for(int s=0;s<2;s++){
                if(!example.psi_N(face)){
                    all_N=false;
                    if(!example.use_multiphase_projection){
                        bool all_phase_dirichlet=true;
                        for(PHASE_ID p=PHASE_ID(0);p<example.phases.m;p++){
                            if(example.phases(p).mass(face)){
                                all_phase_dirichlet=false;}}
                        if(all_phase_dirichlet)
                            dirichlet=true;}}
                face.index(a)++;}}
        if(all_N) example.cell_index(it.index)=pressure_N;
        else if(dirichlet) example.cell_index(it.index)=pressure_D;
        else example.cell_index(it.index)=nvar++;}
    Fix_Periodic(example.cell_index,1);
    example.rhs.v.Remove_All();
    example.rhs.v.Resize(nvar);
    return nvar;
}
//#####################################################################
// Function Compute_Laplacian
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Compute_Laplacian(int nvar)
{
    example.projection_system.A.Reset(nvar);
    ARRAY<int> tmp0,tmp1;
#pragma omp parallel
    {
        bool has_dirichlet=false;
        SPARSE_MATRIX_THREADED_CONSTRUCTION<T> helper(example.projection_system.A,tmp0,tmp1);
        for(CELL_ITERATOR_THREADED<TV> it(example.grid,0);it.Valid();it.Next()){
            int center_index=example.cell_index(it.index);
            if(center_index<0) continue;

            T diag=0;
            helper.Start_Row();
            for(int a=0;a<TV::m;a++){
                FACE_INDEX<TV::m> face(a,it.index);
                for(int s=0;s<2;s++){
                    TV_INT cell=it.index;
                    cell(a)+=2*s-1;
                    PHYSBAM_ASSERT(face.Cell_Index(s)==cell);

                    if(!example.psi_N(face)){
                        T entry=sqr(example.grid.one_over_dX(a))/Density(face);
                        diag+=entry;
                        int ci=example.cell_index(cell);
                        if(ci>=0) helper.Add_Entry(ci,-entry);
                        else{
                            has_dirichlet=true;
                            T rhs=0;
                            if(example.bc_pressure)
                                rhs=entry*example.dt*example.bc_pressure(cell,PHASE_ID(0),example.time);
                            example.rhs.v(center_index)+=rhs;}}
                    face.index(a)++;}}
            helper.Add_Entry(center_index,diag);}
        helper.Finish();
        if(has_dirichlet)
#pragma omp critical
            example.projection_system.dc_present=true;
    }

    if(example.projection_system.use_preconditioner)
        example.projection_system.A.Construct_Incomplete_Cholesky_Factorization();
}
//#####################################################################
// Function Compute_Gradient
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Compute_Gradient(int nvar)
{
    example.projection_system.mass.Remove_All();
    example.projection_system.faces.Remove_All();
    example.projection_system.gradient.Reset(nvar);
    example.projection_system.gradp_bc.Remove_All();
    APPEND_HOLDER<T> mass_h(example.projection_system.mass);
    APPEND_HOLDER<FACE_INDEX<TV::m> > faces_h(example.projection_system.faces);
    APPEND_HOLDER<T> gradp_bc_h(example.projection_system.gradp_bc);
    ARRAY<int> tmp0,tmp1,tmp2,tmp3;
#pragma omp parallel
    {
#pragma omp single
        {
            mass_h.Init();
            faces_h.Init();
            gradp_bc_h.Init();
        }
#pragma omp barrier
        ARRAY<T>& mass_t=mass_h.Array();
        ARRAY<FACE_INDEX<TV::m> >& faces_t=faces_h.Array();
        ARRAY<T>& gradp_bc_t=gradp_bc_h.Array();
        SPARSE_MATRIX_THREADED_CONSTRUCTION<T> G_helper(example.projection_system.gradient,tmp2,tmp3);
        ARRAY<PAIR<PHASE_ID,T> > face_fractions;
        for(FACE_ITERATOR_THREADED<TV> it(example.grid);it.Valid();it.Next()){
            FACE_INDEX<TV::m> face=it.Full_Index();
            int c0=example.cell_index(it.First_Cell_Index());
            int c1=example.cell_index(it.Second_Cell_Index());
            if(example.psi_N(face)){
                Face_Fraction(face,face_fractions);
                T u_star=0;
                for(int j=0;j<face_fractions.m;j++)
                    u_star+=face_fractions(j).y*example.phases(face_fractions(j).x).velocity(face);
                T rhs=example.grid.one_over_dX(it.axis)*u_star;
                if(c0>=0) example.rhs.v(c0)-=rhs;
                if(c1>=0) example.rhs.v(c1)+=rhs;
                continue;}
            if(example.bc_type(it.axis)==example.BC_PERIODIC)
                if(it.index(it.axis)==example.grid.numbers_of_cells(it.axis))
                    continue;

            T mass=Density(face);
            if(!mass) continue;

            G_helper.Start_Row();
            T gradp_bc=0;
            if(c0>=0) G_helper.Add_Entry(c0,-example.grid.one_over_dX(it.axis));
            else if(c0==pressure_D){
                T pr=example.bc_pressure?example.bc_pressure(it.First_Cell_Index(),PHASE_ID(0),example.time):0;
                gradp_bc+=-example.dt*example.grid.one_over_dX(it.axis)*pr;}
            if(c1>=0) G_helper.Add_Entry(c1,example.grid.one_over_dX(it.axis));
            else if(c1==pressure_D){
                T pr=example.bc_pressure?example.bc_pressure(it.Second_Cell_Index(),PHASE_ID(0),example.time):0;
                gradp_bc+=example.dt*example.grid.one_over_dX(it.axis)*pr;}
            gradp_bc_t.Append(gradp_bc);
            faces_t.Append(face);
            mass_t.Append(mass);}
        G_helper.Finish();
    }
    mass_h.Combine();
    faces_h.Combine();
    gradp_bc_h.Combine();
}
//#####################################################################
// Function Compute_Poisson_Matrix
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Compute_Poisson_Matrix()
{
    TIMER_SCOPE_FUNC;
    Compute_Laplacian(example.dof);
    Compute_Gradient(example.dof);
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

    ARRAY<T> tmp(example.projection_system.gradient.m);
    ARRAY<PAIR<PHASE_ID,T> > face_fractions;
    for(int i=0;i<tmp.m;i++){
        FACE_INDEX<TV::m> face(example.projection_system.faces(i));
        Face_Fraction(face,face_fractions);
        T u_star=0;
        for(int j=0;j<face_fractions.m;j++)
            u_star+=face_fractions(j).y*example.phases(face_fractions(j).x).velocity(face);
        tmp(i)=u_star;}
    example.projection_system.gradient.Transpose_Times_Add(tmp,example.rhs.v);

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
        OCTAVE_OUTPUT<T>(LOG::sprintf("C-%i.txt",solve_id).c_str()).Write_Preconditioner("C",system,*example.av(0),*example.av(1));
        OCTAVE_OUTPUT<T>(LOG::sprintf("P-%i.txt",solve_id).c_str()).Write_Projection("P",system,*example.av(0));
        OCTAVE_OUTPUT<T>(LOG::sprintf("b-%i.txt",solve_id).c_str()).Write("b",example.rhs);
        OCTAVE_OUTPUT<T>(LOG::sprintf("G-%i.txt",solve_id).c_str()).Write("G",example.projection_system.gradient);}

    CONJUGATE_GRADIENT<T> cg;
    cg.finish_before_indefiniteness=true;
    cg.relative_tolerance=true;
    bool converged=cg.Solve(example.projection_system,example.sol,example.rhs,
        example.av,example.solver_tolerance,0,example.solver_iterations);
    if(!converged) LOG::printf("SOLVER DID NOT CONVERGE.\n");

    if(example.print_matrix)
        OCTAVE_OUTPUT<T>(LOG::sprintf("x-%i.txt",solve_id).c_str()).Write("x",example.sol);

    if(example.test_system){
        example.projection_system.Multiply(example.sol,*example.av(0));
        *example.av(0)-=example.rhs;
        T r=example.projection_system.Convergence_Norm(*example.av(0));
        LOG::cout<<"residual: "<<r<<std::endl;}

    tmp.Resize(example.projection_system.gradient.m);
    example.projection_system.gradient.Times(example.sol.v,tmp);
    for(int i=0;i<tmp.m;i++){
        FACE_INDEX<TV::m> face(example.projection_system.faces(i));
        Face_Fraction(face,face_fractions);
        T dv=(tmp(i)+example.projection_system.gradp_bc(i))/example.projection_system.mass(i);
        for(int j=0;j<face_fractions.m;j++)
            example.phases(face_fractions(j).x).velocity(face)-=dv;}

    for(PHASE_ID p(0);p<example.phases.m;p++)
        Fix_Periodic(example.phases(p).velocity);
}
//#####################################################################
// Function Apply_Particle_Forces
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Apply_Particle_Forces()
{
#pragma omp parallel for
    for(int i=0;i<example.force.array.m;i++)
        example.force.array(i)=0;
    for(PHASE_ID pid(0);pid<example.phases.m;pid++)
        example.phases(pid).gather_scatter->template Scatter<int>(true,
            [this](int p,const PARTICLE_GRID_FACE_ITERATOR<TV>& it,int data)
            {
                T w=it.Weight();
                FACE_INDEX<TV::m> index=it.Index();
                example.force(index)+=w*example.gravity(index.axis)*example.particles.mass(p);
            });
    Reflect_Boundary(&MPM_MAC_DRIVER::Reflect_Boundary_Particle_Force,RF::ghost|RF::duplicate_corners);
    for(PHASE_ID p(0);p<example.phases.m;p++){
        PHASE& ph=example.phases(p);
#pragma omp parallel for
        for(int i=0;i<ph.velocity.array.m;i++)
            if(ph.mass.array(i))
                ph.velocity.array(i)+=example.force.array(i)/ph.mass.array(i)*example.dt;}
}
//#####################################################################
// Function Apply_Grid_Forces
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Apply_Grid_Forces()
{
#pragma omp parallel for
    for(int i=0;i<example.force.array.m;i++)
        example.force.array(i)=0;
    for(PHASE_ID p(0);p<example.phases.m;p++){
        PHASE& ph=example.phases(p);
        for(int i=0;i<ph.valid_flat_indices.m;i++){
            int k=ph.valid_flat_indices(i);
            FACE_INDEX<TV::m> f=ph.valid_indices(i);
            TV af=example.Compute_Analytic_Force(p,example.grid.Face(f),example.time);
            example.force.array(k)+=af(f.axis);}}
    Reflect_Boundary(&MPM_MAC_DRIVER::Reflect_Boundary_Grid_Force,RF::ghost|RF::delay_corners);
    for(PHASE_ID p(0);p<example.phases.m;p++){
        PHASE& ph=example.phases(p);
#pragma omp parallel for
        for(int i=0;i<ph.velocity.array.m;i++)
            ph.velocity.array(i)+=example.force.array(i)*example.dt;}
}
//#####################################################################
// Function Apply_Forces
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Apply_Forces()
{
    TIMER_SCOPE_FUNC;
    if(example.extrap_type=='p'){
        Apply_Particle_Forces();
        Apply_Grid_Forces();}
    else example.Apply_Forces(example.time);
    for(PHASE_ID p(0);p<example.phases.m;p++)
        Fix_Periodic(example.phases(p).velocity);
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
    T v2=0;
    for(PHASE_ID p(0);p<example.phases.m;p++){
        const PHASE& ph=example.phases(p);
#pragma omp parallel for reduction(max:v2)
        for(int k=0;k<ph.simulated_particles.m;k++){
            int p=ph.simulated_particles(k);
            v2=max(v2,example.particles.V(p).Magnitude_Squared());}}
    return sqrt(v2);
}
//#####################################################################
// Function Grid_V_Upper_Bound
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MAC_DRIVER<TV>::
Grid_V_Upper_Bound() const
{
    T result=0;
    T xi=(T)6*sqrt((T)TV::m)*example.grid.one_over_dX.Min();
    for(PHASE_ID i(0);i<example.phases.m;i++)
#pragma omp parallel for reduction(max:result)
    for(int k=0;k<example.phases(i).simulated_particles.m;k++){
        int p=example.phases(i).simulated_particles(k);
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

    if(!example.particles.store_phase){
        // If particles do not have phase attribute,
        // put them in the default(0) phase.
        PHASE& ph=example.phases(PHASE_ID(0));
        ph.simulated_particles.Remove_All();

        for(int p=0;p<example.particles.number;p++){
            if(example.particles.valid(p))
                ph.simulated_particles.Append(p);}
    }
    else{
        // Otherwise dispatch particles to their corresponding phase
        // based on their phase attribute value.
        for(PHASE_ID i(0);i<example.phases.m;i++)
            example.phases(i).simulated_particles.Remove_All();

        for(int p=0;p<example.particles.number;p++){
            PHASE& ph=example.phases(example.particles.phase(p));
            if(example.particles.valid(p))
                ph.simulated_particles.Append(p);}}
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
// Function Fix_Periodic
//#####################################################################
template<class TV> template<class T2> void MPM_MAC_DRIVER<TV>::
Fix_Periodic(ARRAY<T2,TV_INT>& u,int ghost) const
{
    if(!example.bc_type.Contains(example.BC_PERIODIC)) return;
    if(ghost==INT_MAX) ghost=example.ghost;
    Fill_Ghost_Cells_Periodic(example.grid,u,u,example.periodic_boundary.is_periodic,ghost);
}
//#####################################################################
// Function Fix_Periodic
//#####################################################################
template<class TV> template<class T2> void MPM_MAC_DRIVER<TV>::
Fix_Periodic(ARRAY<T2,FACE_INDEX<TV::m> >& u,int ghost) const
{
    if(!example.bc_type.Contains(example.BC_PERIODIC)) return;
    if(ghost==INT_MAX) ghost=example.ghost;
    Apply_Boundary_Condition_Face_Periodic(example.grid,u,example.periodic_boundary.is_periodic);
    Fill_Ghost_Faces_Periodic(example.grid,u,u,example.periodic_boundary.is_periodic,ghost);
}
//#####################################################################
// Function Fix_Periodic
//#####################################################################
template<class TV> template<class T2> void MPM_MAC_DRIVER<TV>::
Fix_Periodic_Accum(ARRAY<T2,TV_INT>& u,int ghost) const
{
    if(!example.bc_type.Contains(example.BC_PERIODIC)) return;
    if(ghost==INT_MAX) ghost=example.ghost;
    Fill_Ghost_Cells_Periodic_Accum(example.grid,u,u,example.periodic_boundary.is_periodic,ghost);
}
//#####################################################################
// Function Fix_Periodic
//#####################################################################
template<class TV> template<class T2> void MPM_MAC_DRIVER<TV>::
Fix_Periodic_Accum(ARRAY<T2,FACE_INDEX<TV::m> >& u,int ghost) const
{
    if(!example.bc_type.Contains(example.BC_PERIODIC)) return;
    if(ghost==INT_MAX) ghost=example.ghost;
    Fill_Ghost_Faces_Periodic_Accum(example.grid,u,u,example.periodic_boundary.is_periodic,ghost);
}
//#####################################################################
// Function Extrapolate_Velocity
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Extrapolate_Velocity(bool use_bc)
{
    for(PHASE_ID i(0);i<example.phases.m;i++)
        Extrapolate_Velocity(i,use_bc);
}
//#####################################################################
// Function Reflect_Boundary_Particle_Force
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Reflect_Boundary_Particle_Force(const FACE_INDEX<TV::m>& in,const FACE_INDEX<TV::m>& out,int side) const
{
    typename MPM_MAC_EXAMPLE<TV>::BC_TYPE bc_type=example.bc_type(side);
    if((bc_type==example.BC_SLIP && side/2==out.axis) || bc_type==example.BC_NOSLIP){
        example.force(in)-=example.force(out);
        example.force(out)=-example.force(in);}
    else{
        example.force(in)+=example.force(out);
        example.force(out)=example.force(in);}
}
//#####################################################################
// Function Reflect_Boundary_Grid_Force
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Reflect_Boundary_Grid_Force(const FACE_INDEX<TV::m>& in,const FACE_INDEX<TV::m>& out,int side) const
{
    typename MPM_MAC_EXAMPLE<TV>::BC_TYPE bc_type=example.bc_type(side);
    if((bc_type==example.BC_SLIP && side/2==out.axis) || bc_type==example.BC_NOSLIP)
        example.force(out)=-example.force(in);
    else
        example.force(out)=example.force(in);
}
//#####################################################################
// Function Reflect_Boundary_Mass_Momentum
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Reflect_Boundary_Mass_Momentum(const FACE_INDEX<TV::m>& in,const FACE_INDEX<TV::m>& out,int side) const
{
    for(PHASE_ID pid(0);pid<example.phases.m;pid++){
        PHASE& ph=example.phases(pid);
        T& m_in=ph.mass(in),&m_out=ph.mass(out);

        T& moment_in=ph.velocity(in),&moment_out=ph.velocity(out);
        typename MPM_MAC_EXAMPLE<TV>::BC_TYPE bc_type=example.bc_type(side);
        if((bc_type==example.BC_SLIP && side/2==out.axis) || bc_type==example.BC_NOSLIP){
            T bc=0;
            if(example.bc_velocity(side)){
                TV nearest_point=example.grid.Clamp(example.grid.Face(out));
                bc=example.bc_velocity(side)(nearest_point,out.axis,pid,example.time);}
            T mv_in_old=moment_in;
            T mv_out_old=moment_out;
            moment_in+=2*bc*m_out-moment_out;
            moment_out=mv_out_old+2*bc*m_in-mv_in_old;
        }
        else{
            moment_in+=moment_out;
            moment_out=moment_in;}
        m_in+=m_out;
        m_out=m_in;}
}
//#####################################################################
// Function Reflect_Boundary_Copy_Only
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Reflect_Boundary_Velocity_Copy_Only(const FACE_INDEX<TV::m>& in,const FACE_INDEX<TV::m>& out,int side) const
{
    for(PHASE_ID pid(0);pid<example.phases.m;pid++){
        PHASE& ph=example.phases(pid);

        const T& vel_in=ph.velocity(in);
        T& vel_out=ph.velocity(out);
        typename MPM_MAC_EXAMPLE<TV>::BC_TYPE bc_type=example.bc_type(side);
        if((bc_type==example.BC_SLIP && side/2==out.axis) || bc_type==example.BC_NOSLIP){
            T bc=0;
            if(example.bc_velocity(side)){
                TV nearest_point=example.grid.Clamp(example.grid.Face(out));
                bc=example.bc_velocity(side)(nearest_point,out.axis,pid,example.time);}
            vel_out=(2*bc-vel_in);
        }
        else vel_out=vel_in;}
}
//#####################################################################
// Function Reflect_Boundary
//#####################################################################
template<class TV> template<class F> void MPM_MAC_DRIVER<TV>::
Reflect_Boundary(F func,RF flag) const
{
    RANGE<TV_INT> domain=example.grid.Domain_Indices(0);
    TV_INT corner[2]={domain.min_corner,domain.max_corner};
    for(FACE_RANGE_ITERATOR<TV::m> it(domain,example.ghost,0,flag);it.Valid();it.Next()){
        typename MPM_MAC_EXAMPLE<TV>::BC_TYPE bc_type=example.bc_type(it.side);
        if(bc_type==example.BC_PERIODIC) continue;
        int side_axis=it.side/2;
        FACE_INDEX<TV::m> f=it.face;
        f.index(side_axis)=2*corner[it.side%2](side_axis)-it.face.index(side_axis);
        if(side_axis!=it.face.axis) f.index(side_axis)-=1;
        (this->*func)(f,it.face,it.side);}
}
//#####################################################################
// Function Extrapolate_Velocity
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Extrapolate_Velocity(PHASE_ID pid,bool use_bc) const
{
    PHASE& ph=example.phases(pid);
    RANGE<TV_INT> domain=example.grid.Domain_Indices(0);
    RANGE<TV_INT> ghost_domain=example.grid.Domain_Indices(example.ghost);
    TV_INT bound[TV::m];
    for(int i=0;i<TV::m;i++){
        bound[i]=domain.max_corner;
        for(int j=0;j<TV::m;j++) if(i!=j)
            bound[i](j)++;}
    for(FACE_RANGE_ITERATOR<TV::m> it(ghost_domain,domain,RF::ghost|RF::skip_inner|RF::delay_corners);it.Valid();it.Next()){
        typename MPM_MAC_EXAMPLE<TV>::BC_TYPE bc_type=example.bc_type(it.side);
        if(bc_type==example.BC_PERIODIC) continue;
        int side_axis=it.side/2;
        bool normal=side_axis==it.face.axis;
        bool has_bc=(bc_type==example.BC_SLIP && normal)||(bc_type==example.BC_NOSLIP);
        bool is_normal_bc=(bc_type==example.BC_SLIP||bc_type==example.BC_NOSLIP)&&normal;
        char extrap_type=example.extrap_type;
        TV_INT corner=it.side%2?domain.max_corner:domain.min_corner;
        if(bc_type==example.BC_FREE || !has_bc){
            if(extrap_type=='c') extrap_type='C';
            else if(extrap_type=='l') extrap_type='L';}
        if(!normal){
            if(it.face.index.Componentwise_Greater_Equal(TV_INT()).Count_Matches(false)+
               it.face.index.Componentwise_Greater_Equal(bound[side_axis]).Count_Matches(true)>=2){
                if(extrap_type=='c') extrap_type='C';
                else if(extrap_type=='l') extrap_type='L';}}
        int sign_out=it.side%2?1:-1;
        auto bc=[this,it](const TV& X,int axis,PHASE_ID p,T t){
            if(example.bc_velocity(it.side)) return example.bc_velocity(it.side)(X,axis,p,t);
            else return (T)0;};
        auto further=[=](FACE_INDEX<TV::m> f){f.index(side_axis)-=sign_out;return f;};
        FACE_INDEX<TV::m> nearest_face[2]={it.face,it.face};
        nearest_face[0].index(side_axis)=nearest_face[1].index(side_axis)=corner(side_axis);
        nearest_face[1].axis=side_axis;
        if(it.side%2) nearest_face[0]=further(nearest_face[0]);
        TV nearest_point=example.grid.Face(it.face);
        nearest_point(side_axis)=example.grid.Face(nearest_face[1])(side_axis);
        T& u=ph.velocity(it.face);
        switch(extrap_type){
            case 'r':{
                FACE_INDEX<TV::m> f=it.face;
                f.index(side_axis)=2*corner(side_axis)-it.face.index(side_axis);
                if(!normal) f.index(side_axis)-=1;
                u=ph.velocity(f);
                if(has_bc){
                    T v=0;
                    if(use_bc) v=bc(nearest_point,it.face.axis,pid,example.time);
                    else if(normal) v=ph.velocity(nearest_face[normal]);
                    else{
                        FACE_INDEX<TV::m> f0=nearest_face[false];
                        FACE_INDEX<TV::m> f1=further(f0);
                        v=1.5*ph.velocity(f0)-0.5*ph.velocity(f1);}
                    u=2*v-u;}
                break;}
            case 'a': u=bc(example.grid.Face(it.face),it.face.axis,pid,example.time);break;
            case '0': u=0;break;
            case 'c': u=bc(nearest_point,it.face.axis,pid,example.time);break;
            case 'C':{
                FACE_INDEX<TV::m> f=nearest_face[normal];
                if(is_normal_bc) f=further(f);
                u=ph.velocity(f);
                break;}
            case 'l':{
                T u0=bc(nearest_point,it.face.axis,pid,example.time);
                FACE_INDEX<TV::m> f=nearest_face[normal];
                if(is_normal_bc) f=further(f);
                T u1=ph.velocity(f);
                // compute how far from u to u0(given by bc), measured in number of |u0-u1|
                int n=normal?
                    abs(it.face.index(side_axis)-nearest_face[normal].index(side_axis)):
                    2*abs(it.face.index(side_axis)-nearest_face[normal].index(side_axis))-1;
                u=u0+n*(u0-u1);
                break;}
            case 'L':{
                FACE_INDEX<TV::m> f0=nearest_face[normal];
                FACE_INDEX<TV::m> f1=further(f0);
                int n=abs(it.face.index(side_axis)-f0.index(side_axis));
                if(is_normal_bc){
                    f0=f1;
                    f1=further(f0);
                    n+=1;}
                T u0=ph.velocity(f0);
                T u1=ph.velocity(f1);
                u=u0+n*(u0-u1);
                break;}
            default: PHYSBAM_FATAL_ERROR("Unrecognized extrapolation type");}}
}
//#####################################################################
// Function Apply_Viscosity
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Apply_Viscosity()
{
    TIMER_SCOPE_FUNC;
    PHYSBAM_ASSERT(example.phases.m==PHASE_ID(1)); // Single phase for now.
    typedef KRYLOV_VECTOR_WRAPPER<T,ARRAY<T> > VEC;
    typedef MATRIX_SYSTEM<SPARSE_MATRIX_FLAT_MXN<T>,T,VEC> MAT;
    PHASE& ph=example.phases(PHASE_ID(0));
    if(!ph.viscosity) return;

    auto p_psi_D=[this](const TV_INT& c){return example.cell_index(c)==pressure_D;};

    // Based on logic from IMPLICIT_VISCOSITY_UNIFORM::Setup_Boundary_Conditions
    for(int axis=0;axis<TV::m;axis++){
        GRID<TV> face_grid(example.grid.Get_Face_MAC_Grid(axis));
        ARRAY<int,TV_INT> velocity_index(face_grid.Domain_Indices(1),true,-1);
        VEC rhs,sol,tmp_vector;

        // Allocate velocities that will be corrected; copy over initial guess
        for(CELL_ITERATOR<TV> it(face_grid);it.Valid();it.Next()){
            FACE_INDEX<TV::m> face(axis,it.index);
            if(example.bc_type(axis)==example.BC_PERIODIC)
                if(it.index(axis)==example.grid.numbers_of_cells(axis))
                    continue;
            if(p_psi_D(face.First_Cell_Index()) && p_psi_D(face.Second_Cell_Index()))
                continue;
            if(!example.psi_N(face))
                velocity_index(it.index)=rhs.v.Append(ph.velocity(face));}
        Fix_Periodic(velocity_index,1);

        // Compute as needed rather than store.  face_grid index.
        auto psi_N=[axis,this,p_psi_D](const FACE_INDEX<TV::m>& face,int s)
        {
            TV_INT a=face.index,b=a;
            b(axis)--;
            if(face.axis==axis) return p_psi_D(b);
            int index=face.index(face.axis);
            bool boundary=index==example.grid.Domain_Indices().min_corner(face.axis)||
                index==example.grid.Domain_Indices().max_corner(face.axis);
            if(boundary && example.bc_type(2*face.axis+s)==example.BC_SLIP) return true;
            if(p_psi_D(a) && p_psi_D(b)) return true;
            a(face.axis)--;
            b(face.axis)--;
            return p_psi_D(a) && p_psi_D(b);
        };

        SPARSE_MATRIX_FLAT_MXN<T> A;
        MAT sys(A);
        TV scale=ph.viscosity/ph.density*example.dt*sqr(example.grid.one_over_dX);

        A.Reset(rhs.v.m);
        sol.v=rhs.v;
        tmp_vector.v.Resize(rhs.v.m);
        ARRAY<int> tmp0,tmp1;
#pragma omp parallel
        {
            ARRAY<T,PHASE_ID> u_bc(example.phases.m);
            SPARSE_MATRIX_THREADED_CONSTRUCTION<T> helper(A,tmp0,tmp1);
            for(CELL_ITERATOR_THREADED<TV> it(face_grid,0);it.Valid();it.Next()){
                int center_index=velocity_index(it.index);
                if(center_index<0) continue;
                if(example.bc_type(axis)==example.BC_PERIODIC)
                    if(it.index(axis)==example.grid.numbers_of_cells(axis))
                        continue;

                T diag=1;
                helper.Start_Row();
                for(int a=0;a<TV::m;a++){
                    FACE_INDEX<TV::m> face(a,it.index);
                    for(int s=0;s<2;s++){
                        TV_INT cell=it.index;
                        cell(a)+=2*s-1;
                        PHYSBAM_ASSERT(face.Cell_Index(s)==cell);
                        if(psi_N(face,s)){
                            // Assume for now that the boundary condition is zero.
                        }
                        else{
                            diag+=scale(a);
                            int ci=velocity_index(cell);
                            if(ci>=0) helper.Add_Entry(ci,-scale(a));
                            else rhs.v(center_index)+=scale(a)*ph.velocity(FACE_INDEX<TV::m>(axis,cell));}
                        face.index(a)++;}}
                helper.Add_Entry(center_index,diag);}
            helper.Finish();
        }
        if(example.projection_system.use_preconditioner){
            A.Construct_Incomplete_Cholesky_Factorization();
            sys.Set_Preconditioner(*A.C,tmp_vector);}

        ARRAY<KRYLOV_VECTOR_BASE<T>*> av;
        static int solve_id=-1;
        solve_id++;
        if(example.test_system) sys.Test_System(sol);
        if(example.print_matrix){
            LOG::cout<<"solve id: "<<solve_id<<std::endl;
            KRYLOV_SOLVER<T>::Ensure_Size(av,sol,2);
            OCTAVE_OUTPUT<T>(LOG::sprintf("visc-M-%i.txt",solve_id).c_str()).Write("M",sys,*av(0),*av(1));
            OCTAVE_OUTPUT<T>(LOG::sprintf("visc-C-%i.txt",solve_id).c_str()).Write_Preconditioner("C",sys,*av(0),*av(1));
            OCTAVE_OUTPUT<T>(LOG::sprintf("visc-b-%i.txt",solve_id).c_str()).Write("b",rhs);}
        CONJUGATE_GRADIENT<T> cg;
        cg.finish_before_indefiniteness=true;
        cg.relative_tolerance=true;
        bool converged=cg.Solve(sys,sol,rhs,
            av,example.solver_tolerance,0,example.solver_iterations);
        if(!converged) LOG::printf("SOLVER DID NOT CONVERGE.\n");

        if(example.print_matrix)
            OCTAVE_OUTPUT<T>(LOG::sprintf("visc-x-%i.txt",solve_id).c_str()).Write("x",sol);

        if(example.test_system){
            sys.Multiply(sol,*av(0));
            *av(0)-=rhs;
            T r=sys.Convergence_Norm(*av(0));
            LOG::cout<<"visc residual: "<<r<<std::endl;}

        for(CELL_ITERATOR<TV> it(face_grid);it.Valid();it.Next()){
            FACE_INDEX<TV::m> face(axis,it.index);
            int i=velocity_index(it.index);
            if(i>=0) ph.velocity(face)=sol.v(i);}}

    Fix_Periodic(ph.velocity);
}
//#####################################################################
// Function Step
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Step(std::function<void()> func,const char* name,bool dump_substep,bool do_step)
{
    auto& p=example.time_step_callbacks.Get_Or_Insert(name);
    p.x=true; // Flag the callback name as recognized, for sanity checking later
    if(!do_step) return;
    for(int i=0;i<p.y(1).m;i++) p.y(1)(i)();
    func();
    // Note: p may be invalidated by func(), so we must re-access it here.
    const auto& q=example.time_step_callbacks.Get_Or_Insert(name);
    for(int i=0;i<q.y(0).m;i++) q.y(0)(i)();
    if(dump_substep) PHYSBAM_DEBUG_WRITE_SUBSTEP(name,1);
}
//#####################################################################
// Function Invalidate_Particle
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Invalidate_Particle(int p)
{
    example.particles.valid(p)=false;
    example.particles.Add_To_Deletion_List(p);
}
//#####################################################################
// Function Move_Particles
//#####################################################################
template<class TV> void MPM_MAC_DRIVER<TV>::
Move_Particles()
{
    auto Clip=[this](int p,int data=0)
        {
            for(int i=0;i<TV::m;i++){
                T &x=example.particles.X(p)(i),a=example.grid.domain.min_corner(i),b=example.grid.domain.max_corner(i);
                if(x<a){
                    typename MPM_MAC_EXAMPLE<TV>::BC_TYPE bc_type=example.bc_type(2*i);
                    bool wall=bc_type==example.BC_SLIP||bc_type==example.BC_NOSLIP;
                    if(bc_type==example.BC_PERIODIC) x=wrap(x,a,b);
                    else if(bc_type==example.BC_FREE) Invalidate_Particle(p);
                    else if(example.clamp_particles && wall) x=a;}
                if(x>b){
                    typename MPM_MAC_EXAMPLE<TV>::BC_TYPE bc_type=example.bc_type(2*i+1);
                    bool wall=bc_type==example.BC_SLIP||bc_type==example.BC_NOSLIP;
                    if(bc_type==example.BC_PERIODIC) x=wrap(x,a,b);
                    else if(bc_type==example.BC_FREE) Invalidate_Particle(p);
                    else if(example.clamp_particles && wall) x=b;}}
        };

    if(example.rk_particle_order){
        LINEAR_INTERPOLATION_MAC<TV,T> li(example.grid);
#pragma omp parallel for
        for(int p=0;p<example.particles.X.m;p++)
            if(example.particles.valid(p)){
                const PHASE& ph=example.phases(example.particles.phase(p));
                for(RUNGEKUTTA<TV> rk(example.particles.X(p),example.rk_particle_order,example.dt,0);rk.Valid();rk.Next())
                    example.particles.X(p)+=example.dt*li.Clamped_To_Array(ph.velocity,example.particles.X(p));
                Clip(p);}}
    else if(example.flip){
        for(PHASE_ID i(0);i<example.phases.m;i++){
            PHASE& ph=example.phases(i);
            ph.gather_scatter->template Gather<int>(false,
                [this,&ph](int p,const PARTICLE_GRID_FACE_ITERATOR<TV>& it,int data)
                {example.particles.X(p)(it.Index().axis)+=example.dt*it.Weight()*ph.velocity(it.Index());},
                Clip);}}
    else{
#pragma omp parallel for
        for(int p=0;p<example.particles.X.m;p++)
            if(example.particles.valid(p)){
                example.particles.X(p)+=example.dt*example.particles.V(p);
                Clip(p);}}
}
//#####################################################################
namespace PhysBAM{
template class MPM_MAC_DRIVER<VECTOR<float,2> >;
template class MPM_MAC_DRIVER<VECTOR<float,3> >;
template class MPM_MAC_DRIVER<VECTOR<double,2> >;
template class MPM_MAC_DRIVER<VECTOR<double,3> >;
}
