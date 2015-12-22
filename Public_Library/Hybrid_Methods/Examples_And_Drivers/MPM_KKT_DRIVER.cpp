//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Log/DEBUG_SUBSTEPS.h>
#include <Tools/Log/SCOPE.h>
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_KKT_DRIVER.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_KKT_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <Hybrid_Methods/System/MPM_KKT_KRYLOV_SYSTEM.h>
#include <Hybrid_Methods/System/MPM_KKT_KRYLOV_VECTOR.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
using namespace PhysBAM;
namespace{
    template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
    {
        ((MPM_KKT_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
    }
};
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_KKT_DRIVER<TV>::
MPM_KKT_DRIVER(MPM_KKT_EXAMPLE<TV>& example)
    :example(example),
    kkt_sys(*new MPM_KKT_KRYLOV_SYSTEM<TV>(example)),
    kkt_lhs(*new MPM_KKT_KRYLOV_VECTOR<TV>(example.valid_velocity_indices,example.valid_pressure_indices)),
    kkt_rhs(*new MPM_KKT_KRYLOV_VECTOR<TV>(example.valid_velocity_indices,example.valid_pressure_indices))
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_KKT_DRIVER<TV>::
~MPM_KKT_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
    av.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Execute_Main_Program
//#####################################################################
template<class TV> void MPM_KKT_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void MPM_KKT_DRIVER<TV>::
Initialize()
{
    LOG::cout<<std::setprecision(16)<<std::endl;
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.substeps_delay_frame<0?example.write_substeps_level:-1);

    // setup time
    output_number=current_frame=example.restart;

    example.Initialize();
    if(example.restart)
        example.Read_Output_Files(example.restart);

    example.mass.Resize(example.grid.Domain_Indices(example.ghost));
    example.mass_coarse.Resize(example.coarse_grid.Domain_Indices(example.ghost));
    example.velocity.Resize(example.grid.Domain_Indices(example.ghost));
    example.velocity_new.Resize(example.grid.Domain_Indices(example.ghost));
    example.inv_valid_pressure_cell.Resize(example.coarse_grid.Domain_Indices(example.ghost));
    example.inv_valid_velocity_cell.Resize(example.grid.Domain_Indices(example.ghost));
    example.particles.Store_B(example.use_affine);
    example.particles.Store_S(example.use_oldroyd);
    example.current_velocity=&example.velocity;
    PHYSBAM_ASSERT(!example.particles.store_B || !example.particles.store_C);

    kkt_lhs.u.Resize(example.grid.Domain_Indices(example.ghost));
    kkt_lhs.p.Resize(example.coarse_grid.Domain_Indices(example.ghost));
    //kkt_lhs.lambda.Resize(kkt_lhs.p.array.m);
    kkt_rhs.u.Resize(example.grid.Domain_Indices(example.ghost));
    kkt_rhs.p.Resize(example.coarse_grid.Domain_Indices(example.ghost));
    //kkt_rhs.lambda.Resize(kkt_rhs.p.array.m);
    example.one_over_lambda.Resize(example.coarse_grid.Domain_Indices(example.ghost));
    example.J.Resize(example.coarse_grid.Domain_Indices(example.ghost));

    RANGE<TV_INT> range(example.grid.Cell_Indices(example.ghost));
    RANGE<TV_INT> coarse_range(example.coarse_grid.Cell_Indices(example.ghost));
    Initialize_Location(range,example.grid,example.location);
    Initialize_Location(coarse_range,example.coarse_grid,example.coarse_location);
    
    if(!example.restart) Write_Output_Files(0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after init",0,1);
}
//#####################################################################
// Function Initialize_Location 
//#####################################################################
template<class TV> void MPM_KKT_DRIVER<TV>::
Initialize_Location(const RANGE<TV_INT>& range,const GRID<TV>& grid,ARRAY<TV,TV_INT>& location)
{
    location.Resize(range,false,false);
#pragma omp parallel for
    for(int t=0;t<example.threads;t++){
        int a=(range.max_corner.x-range.min_corner.x)*t/example.threads+range.min_corner.x;
        int b=(range.max_corner.x-range.min_corner.x)*(t+1)/example.threads+range.min_corner.x;
        RANGE<TV_INT> local_range(range);
        local_range.max_corner.x=a;
        int i=local_range.Size();
        local_range.min_corner.x=a;
        local_range.max_corner.x=b;
        for(RANGE_ITERATOR<TV::m> it(local_range);it.Valid();it.Next())
            location.array(i++)=grid.Center(it.index);}
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class TV> void MPM_KKT_DRIVER<TV>::
Advance_One_Time_Step()
{
    example.Begin_Time_Step(example.time);
    Update_Simulated_Particles();
    Print_Particle_Stats("particle state",example.dt);
    Update_Particle_Weights();
    example.gather_scatter.Prepare_Scatter(example.particles);
    example.gather_scatter_coarse.Prepare_Scatter(example.particles);
    Particle_To_Grid();
    Print_Grid_Stats("after particle to grid",example.dt,example.velocity,0);
    Print_Energy_Stats("after particle to grid",example.velocity);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after particle to grid",0,1);
    Solve_KKT_System();
    Grid_To_Particle();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after grid to particle",0,1);

    example.End_Time_Step(example.time);
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void MPM_KKT_DRIVER<TV>::
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
            LOG::cout<<"substep dt: "<<example.dt<<std::endl;
            T next_time=example.time+example.dt;
            if(next_time>time_at_frame){
                next_time=time_at_frame;
                done=true;}
            else if(next_time+example.dt>time_at_frame) next_time=(example.time+time_at_frame)/2;
            example.dt=next_time-example.time;

            Advance_One_Time_Step();
            example.time=next_time;}
        example.End_Frame(current_frame);
        Write_Output_Files(++output_number);}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void MPM_KKT_DRIVER<TV>::
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
template<class TV> void MPM_KKT_DRIVER<TV>::
Write_Output_Files(const int frame)
{
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
template<class TV> void MPM_KKT_DRIVER<TV>::
Update_Particle_Weights()
{
    example.weights->Update(example.particles.X);
    example.coarse_weights->Update(example.particles.X);
}
//#####################################################################
// Function Particle_To_Grid
//#####################################################################
template<class TV> void MPM_KKT_DRIVER<TV>::
Particle_To_Grid()
{
    MPM_PARTICLES<TV>& particles=example.particles;
    example.current_velocity=&example.velocity;

#pragma omp parallel for
    for(int i=0;i<example.mass.array.m;i++){
        example.mass.array(i)=0;
        example.velocity.array(i)=TV();
        example.velocity_new.array(i)=TV();}

#pragma omp parallel for
    for(int i=0;i<example.mass_coarse.array.m;i++){
        example.mass_coarse.array(i)=0;
        example.one_over_lambda.array(i)=0;
        example.J.array(i)=0;}

    T scale=1;
    if(particles.store_B && !example.weights->use_gradient_transfer)
        scale=example.weights->Constant_Scalar_Inverse_Dp();
    bool use_gradient=example.weights->use_gradient_transfer;
    ARRAY_VIEW<MATRIX<T,TV::m> > dV(particles.B);
    example.gather_scatter.template Scatter<int>(
        [this,scale,&particles,use_gradient,dV](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
        {
            T w=it.Weight();
            TV_INT index=it.Index();
            example.mass(index)+=w*particles.mass(p);
            TV V=w*particles.V(p);
            if(example.use_affine){
                if(use_gradient) V+=particles.B(p)*it.Gradient();
                else V+=dV(p)*(w*scale*(example.grid.Center(index)-particles.X(p)));}
            example.velocity(index)+=particles.mass(p)*V;
        },true);

    example.gather_scatter_coarse.template Scatter<int>(
        [this,&particles](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
        {
            T w=it.Weight();
            TV_INT index=it.Index();
            example.mass_coarse(index)+=w*particles.mass(p);
            if(particles.lambda(p)!=FLT_MAX) example.one_over_lambda(index)+=w*particles.mass(p)/particles.lambda(p);
            example.J(index)+=w*particles.mass(p)*particles.F(p).Determinant();
        },true);

    example.valid_velocity_indices.Remove_All();
    example.valid_velocity_cell_indices.Remove_All();
    example.valid_pressure_indices.Remove_All();
    example.valid_pressure_cell_indices.Remove_All();
    example.inv_valid_pressure_cell.Fill(-1);
    example.inv_valid_velocity_cell.Fill(-1);
    for(RANGE_ITERATOR<TV::m> it(example.mass.domain);it.Valid();it.Next()){
        int i=example.mass.Standard_Index(it.index);
        if(example.mass.array(i)){
            example.inv_valid_velocity_cell.array(i)=example.valid_velocity_cell_indices.m;
            example.valid_velocity_indices.Append(i);
            example.valid_velocity_cell_indices.Append(it.index);
            example.velocity.array(i)/=example.mass.array(i);}
        else example.velocity.array(i)=TV();}
            
    for(RANGE_ITERATOR<TV::m> it(example.mass_coarse.domain);it.Valid();it.Next()){
        int i=example.mass_coarse.Standard_Index(it.index);
        if(example.mass_coarse.array(i)){
            example.inv_valid_pressure_cell.array(i)=example.valid_pressure_cell_indices.m;
            example.valid_pressure_indices.Append(i);
            example.valid_pressure_cell_indices.Append(it.index);
            example.one_over_lambda.array(i)/=example.mass_coarse.array(i);
            example.J.array(i)/=example.mass_coarse.array(i);}}
    
    // Add ghost pressure and velocity
    ARRAY<TV_INT> ghost_velocity,ghost_pressure;
    RANGE<TV_INT> local_range(TV_INT::Constant_Vector(-1),TV_INT::Constant_Vector(2));
    int initial_velocity_size=example.valid_velocity_indices.m;
    int initial_pressure_size=example.valid_pressure_indices.m;
    LOG::printf("initial_velocity_size=%P\n",initial_velocity_size);
    LOG::printf("initial_presure_size=%P\n",initial_pressure_size);
    for(int i=0;i<initial_velocity_size;i++){
        TV_INT vid=example.valid_velocity_cell_indices(i);
        for(RANGE_ITERATOR<TV::m> it(local_range);it.Valid();it.Next()){
            TV_INT index_possible_ghost_velocity=it.index+vid;
            if(example.inv_valid_velocity_cell(index_possible_ghost_velocity)==-1){
                example.inv_valid_velocity_cell(index_possible_ghost_velocity)=example.valid_velocity_cell_indices.m;
                example.valid_velocity_cell_indices.Append(index_possible_ghost_velocity);
                example.valid_velocity_indices.Append(example.velocity.Standard_Index(index_possible_ghost_velocity));}}
        TV_INT pressure_min_corner=vid/2;
        TV_INT pressure_max_corner=pressure_min_corner+TV_INT::Constant_Vector(2);
        for(int axis=0;axis<TV::m;axis++)
            if(vid(axis)%2==0){
                pressure_min_corner(axis)--;}
        RANGE<TV_INT> pressure_range(pressure_min_corner,pressure_max_corner);
        for(RANGE_ITERATOR<TV::m> pit(pressure_range);pit.Valid();pit.Next()){
            TV_INT pid=pit.index;
            if(example.inv_valid_pressure_cell(pid)==-1){
                example.inv_valid_pressure_cell(pid)=example.valid_pressure_cell_indices.m;
                example.valid_pressure_cell_indices.Append(pid);
                example.valid_pressure_indices.Append(example.mass_coarse.Standard_Index(pid));}}}
    example.valid_pressure_indices.Sort();
    example.valid_velocity_indices.Sort();
    auto compare=[](TV_INT a,TV_INT b)->bool{
        for(int i=0;i<TV::m;i++){
            if(a(i)<b(i)) return true;
            if(a(i)>b(i)) return false;}
        return false;};
    example.valid_velocity_cell_indices.Sort(compare);
    example.valid_pressure_cell_indices.Sort(compare);
    for(int ii=0;ii<example.valid_velocity_indices.m;ii++){
        int int_id=example.valid_velocity_indices(ii);
        example.inv_valid_velocity_cell.array(int_id)=ii;}
    for(int ii=0;ii<example.valid_pressure_indices.m;ii++){
        int int_id=example.valid_pressure_indices(ii);
        example.inv_valid_pressure_cell.array(int_id)=ii;}
    kkt_lhs.lambda.Resize(example.valid_pressure_indices.m);
    kkt_rhs.lambda.Resize(example.valid_pressure_indices.m);
}
//#####################################################################
// Function Grid_To_Particle
//#####################################################################
template<class TV> void MPM_KKT_DRIVER<TV>::
Grid_To_Particle()
{
    MPM_PARTICLES<TV>& particles=example.particles;
    T dt=example.dt;

#pragma omp parallel for
    for(int tid=0;tid<example.gather_scatter.threads;++tid){
        int a=example.simulated_particles.m*tid/example.gather_scatter.threads;
        int b=example.simulated_particles.m*(tid+1)/example.gather_scatter.threads;
        typename PARTICLE_GRID_ITERATOR<TV>::SCRATCH scratch;

        for(int k=a;k<b;k++){
            int p=example.simulated_particles(k);
            TV Vn_interpolate,V_pic,V_flip=particles.V(p);
            MATRIX<T,TV::m> B,grad_Vp,D;
            T J_pic=0;
            
            for(PARTICLE_GRID_ITERATOR<TV> it(example.coarse_weights,p,true,scratch);it.Valid();it.Next()){
                T w=it.Weight();
                TV_INT index=it.Index();
                PHYSBAM_ASSERT(example.inv_valid_pressure_cell(index)!=-1);
                T p_grid=kkt_lhs.p(index);
                J_pic+=w*(-p_grid*example.one_over_lambda(index)+1);}

            for(PARTICLE_GRID_ITERATOR<TV> it(example.weights,p,true,scratch);it.Valid();it.Next()){
                T w=it.Weight();
                TV_INT index=it.Index();
                TV V_grid=example.velocity_new(index);
                V_pic+=w*V_grid;
                V_flip+=w*(V_grid-example.velocity(index));
                Vn_interpolate+=w*example.velocity(index);
                if(example.use_midpoint)
                    V_grid=(T).5*(V_grid+example.velocity(index));
                grad_Vp+=MATRIX<T,TV::m>::Outer_Product(V_grid,it.Gradient());}
            MATRIX<T,TV::m> A=dt*grad_Vp+1;
            //particles.F(p)=A*particles.F(p);
            //T scale_J=J_pic/particles.F(p).Determinant();
            //particles.F(p)*=pow(abs(scale_J),1.0/TV::m);
            //if(scale_J<0) LOG::printf("WARNING: scale_J is negative!\n");
            if(particles.store_S){
                T k=example.dt*example.inv_Wi;
                particles.S(p)=(SYMMETRIC_MATRIX<T,TV::m>::Conjugate(A,particles.S(p))+k)/(1+k);}
            if(example.use_affine && example.use_early_gradient_transfer)
                B=grad_Vp/example.weights->Constant_Scalar_Inverse_Dp();
            else if(example.use_affine)
                for(PARTICLE_GRID_ITERATOR<TV> it(example.weights,p,false,scratch);it.Valid();it.Next()){
                    TV_INT index=it.Index();
                    TV V_grid=example.velocity_new(index);
                    TV Z=example.grid.Center(index);
                    TV xi_new,xp_new;
                    if(example.use_midpoint){
                        xi_new=Z+dt/2*(V_grid+example.velocity(index));
                        xp_new=particles.X(p)+dt/2*(Vn_interpolate+V_pic);}
                    else{
                        xi_new=Z+dt*V_grid;
                        xp_new=particles.X(p)+dt*V_pic;}
                    B+=it.Weight()/2*(MATRIX<T,TV::m>::Outer_Product(V_grid,Z-particles.X(p)+xi_new-xp_new)
                        +MATRIX<T,TV::m>::Outer_Product(Z-particles.X(p)-xi_new+xp_new,V_grid));
                    if(particles.store_C && example.weights->Order()>1) D+=it.Weight()*MATRIX<T,TV::m>::Outer_Product(Z-particles.X(p),Z-particles.X(p));}

            particles.V(p)=V_pic;
            Perform_Particle_Collision(p,example.time+example.dt);
            if(particles.store_B) particles.B(p)=B;
            if(particles.store_C) particles.C(p)=example.weights->Order()==1?grad_Vp:B*D.Inverse();
            if(example.use_midpoint) particles.X(p)+=(particles.V(p)+Vn_interpolate)*(dt/2);
            else particles.X(p)+=particles.V(p)*dt;
            particles.V(p)=V_flip*example.flip+V_pic*(1-example.flip);

            if(!example.grid.domain.Lazy_Inside(particles.X(p))) particles.valid(p)=false;}}
}
//#####################################################################
// Function Solve_KKT_System
//#####################################################################
template<class TV> void MPM_KKT_DRIVER<TV>::
Solve_KKT_System()
{
    bool have_non_zero=false;
    kkt_sys.Build_Div_Matrix();
    kkt_sys.Build_B_Matrix(have_non_zero);
    Build_FEM_Mass_Matrix();
    kkt_lhs.u.Fill(TV());kkt_rhs.u.Fill(TV());
    kkt_lhs.p.Fill((T)0);kkt_rhs.p.Fill((T)0);
    kkt_lhs.lambda.Fill((T)0);kkt_rhs.lambda.Fill((T)0);

    // Add gravity
    example.Capture_Stress();
    example.Precompute_Forces(example.time,example.dt,false);
    example.Add_Forces(kkt_rhs.u,example.time);
    const T pressure_volume_scale=example.coarse_grid.DX().Product();
    if(example.use_FEM_mass){
        int index=example.M.offsets(0);
        for(int row=0;row<example.M.m;row++){
            int end=example.M.offsets(row+1); TV sum=TV();
            for(;index<end;index++){
                int col=example.M.A(index).j;
                sum+=example.M.A(index).a*example.velocity.array(col);}
            kkt_rhs.u.array(example.valid_velocity_indices(row))+=sum/example.dt;}}
    else{
        for(int t=0;t<example.valid_velocity_indices.m;t++){
            int id=example.valid_velocity_indices(t);
            kkt_rhs.u.array(id)+=example.mass.array(id)*(example.velocity.array(id)/example.dt);}}
    for(int t=0;t<example.valid_pressure_indices.m;t++){
        int id=example.valid_pressure_indices(t);
        if(example.mass_coarse.array(id))
            kkt_rhs.p.array(id)=(example.J.array(id)-(T)1)*pressure_volume_scale/(example.dt*example.J.array(id));}
    // MINRES 
    MINRES<T> mr;
    mr.print_diagnostics=true;
    mr.print_residuals=false;
    mr.Solve(kkt_sys,kkt_lhs,kkt_rhs,av,1e-12,0,example.solver_iterations);
    for(int i=0;i<example.av.m;i++){
        (*av(i))*=0;}
    example.velocity_new=kkt_lhs.u;
    LOG::printf("example.velocity_new=%P\n",example.velocity_new);
    static int n=0;
    n++;
    KRYLOV_VECTOR_BASE<T> *a,*b;
    a=kkt_lhs.Clone_Default();
    b=kkt_lhs.Clone_Default();
    OCTAVE_OUTPUT<T> oo(LOG::sprintf("kmatrix-%i.dat",n).c_str());
    oo.Write("K",kkt_sys,*a,*b);
    oo.Write("lhs",kkt_lhs);
    oo.Write("rhs",kkt_rhs);
    exit(0);
//LOG::printf("kkt_lhs.lambda=%P\n",kkt_lhs.lambda);
//LOG::printf("kkt_lhs.u=%P\n",kkt_lhs.u);
//LOG::printf("kkt_lhs.p=%P\n",kkt_lhs.p);
//LOG::printf("example.velocity=%P\n",example.velocity);
//LOG::printf("kkt_rhs.u=%P\n",kkt_rhs.u);
//LOG::printf("example.masss=%P\n",example.mass);
//kkt_sys.Test_System(*a);
//exit(0); 
}
//#####################################################################
// Function Apply_Forces
//#####################################################################
template<class TV> void MPM_KKT_DRIVER<TV>::
Apply_Forces()
{
}
//#####################################################################
// Function Perform_Particle_Collision
//#####################################################################
template<class TV> void MPM_KKT_DRIVER<TV>::
Perform_Particle_Collision(int p,T time)
{
    if(!example.use_particle_collision) return;;
    for(int i=0;i<example.collision_objects.m;i++){
        TV X=example.particles.X(p);
        MPM_COLLISION_OBJECT<TV>* io=example.collision_objects(i);
        T phi=io->Phi(X,time);
        if(phi>=0) continue;
        if(example.collision_objects(i)->type==COLLISION_TYPE::stick) return;
        X-=phi*io->Normal(X,time);
        example.particles.X(p)=X;}
}
//#####################################################################
// Function Apply_Friction
//#####################################################################
template<class TV> void MPM_KKT_DRIVER<TV>::
Apply_Friction()
{
}
//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> typename TV::SCALAR MPM_KKT_DRIVER<TV>::
Compute_Dt() const
{
    T critical_speed=example.cfl*example.grid.DX().Min()/example.max_dt;
    T v=Grid_V_Upper_Bound();
    return (v>critical_speed)?(example.cfl*example.grid.DX().Min()/v):example.max_dt;
}
//#####################################################################
// Function Max_Particle_Speed
//#####################################################################
template<class TV> typename TV::SCALAR MPM_KKT_DRIVER<TV>::
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
template<class TV> typename TV::SCALAR MPM_KKT_DRIVER<TV>::
Grid_V_Upper_Bound() const
{
    if(!example.use_affine || !example.weights->constant_scalar_inertia_tensor) return Max_Particle_Speed();
    T result=0;
    T xi=(T)6*sqrt((T)TV::m)*example.grid.One_Over_DX().Min();
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
template<class TV> void MPM_KKT_DRIVER<TV>::
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

    for(int i=0;i<example.lagrangian_forces.m;i++)
        example.lagrangian_forces(i)->Update_Mpi(example.particle_is_simulated,0);
}
//#####################################################################
// Function Print_Grid_Stats
//#####################################################################
template<class TV> void MPM_KKT_DRIVER<TV>::
Print_Grid_Stats(const char* str,T dt,const ARRAY<TV,TV_INT>& u,const ARRAY<TV,TV_INT>* u0)
{
    if(!example.print_stats) return;
    typename TV::SPIN am=example.Total_Grid_Angular_Momentum(dt,u,u0);
    TV lm=example.Total_Grid_Linear_Momentum(u);
    T ke=example.Total_Grid_Kinetic_Energy(u);
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
template<class TV> void MPM_KKT_DRIVER<TV>::
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
template<class TV> void MPM_KKT_DRIVER<TV>::
Print_Energy_Stats(const char* str,const ARRAY<TV,TV_INT>& u)
{
    if(!example.print_stats) return;
    example.Capture_Stress();
    example.Precompute_Forces(example.time,example.dt,false);
    T ke=example.Total_Grid_Kinetic_Energy(u);
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
// Function Get_Appropriate_Mass_Element
//#####################################################################
template<class TV> ARRAY<ARRAY<typename TV::SCALAR,VECTOR<int,TV::m> >,VECTOR<int,TV::m> >
Get_Appropriate_Mass_Element(typename TV::SCALAR dx_scale)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    const T mass_element[][3]={{1.0/20,13.0/120,1.0/120},{13.0/120,9.0/20,13.0/120},{1.0/120,13.0/120,1.0/20}};
RANGE<TV_INT> velocity_range(TV_INT(),TV_INT::Constant_Vector(3));
    ARRAY<ARRAY<T,TV_INT>,TV_INT> wt(velocity_range);
    for(RANGE_ITERATOR<TV::m> vit1(velocity_range);vit1.Valid();vit1.Next()){
        wt(vit1.index).Resize(velocity_range);
        wt(vit1.index).Fill(dx_scale);
        for(RANGE_ITERATOR<TV::m> vit2(velocity_range);vit2.Valid();vit2.Next()){
            for(int paxis=0;paxis<TV::m;paxis++){
                int row=vit1.index(paxis);
                int col=vit2.index(paxis);
                wt(vit1.index)(vit2.index)*=mass_element[row][col];}}}
    return wt;
}
//#####################################################################
// Function Build_FEM_Mass_Matrix
//#####################################################################
template<class TV> void MPM_KKT_DRIVER<TV>::
Build_FEM_Mass_Matrix()
{
    T dx_scale=example.grid.DX().Product();
    ARRAY<ARRAY<T,TV_INT>,TV_INT> wt=Get_Appropriate_Mass_Element<TV>(dx_scale);
    example.M.Reset(example.velocity.array.m);
    LOG::printf("example.velocity.array.m=%P\n",example.velocity.array.m);
    RANGE<TV_INT> local_range(TV_INT::Constant_Vector(-2),TV_INT::Constant_Vector(3));
    for(int row=0;row<example.valid_velocity_indices.m;row++){
        TV_INT vid_cell=example.valid_velocity_cell_indices(row);
        for(RANGE_ITERATOR<TV::m> it(local_range);it.Valid();it.Next())
            example.M.Append_Entry_To_Current_Row(example.velocity.Standard_Index(vid_cell+it.index),0);
        example.M.Finish_Row();}
    example.M.Sort_Entries();
    
    for(int element=0;element<example.valid_velocity_indices.m;element++){
        int eid=example.valid_velocity_indices(element);
        if(example.mass.array(eid)==0) continue;
        TV v_loc=example.location.array(eid);
        bool skip=false;
        for(int fw=0;fw<example.fluid_walls.m;fw++)
            if(example.fluid_walls(fw)->Lazy_Inside(v_loc)) skip=true;
        if(skip) continue;
        TV_INT eid_cell=example.valid_velocity_cell_indices(element);
        for(RANGE_ITERATOR<TV::m> rit(wt.domain);rit.Valid();rit.Next()){
            TV_INT row_index=eid_cell+rit.index-1;
            int row=example.inv_valid_velocity_cell(row_index);
            if(row==-1) continue;
            for(RANGE_ITERATOR<TV::m> cit(wt(rit.index).domain);cit.Valid();cit.Next()){
                TV_INT col_index=eid_cell+cit.index-1;
                if(example.inv_valid_velocity_cell(col_index)!=-1){
                    int col=example.velocity.Standard_Index(col_index);
                    example.M(row,col)+=wt(rit.index)(cit.index);}}}}
    OCTAVE_OUTPUT<T> oo("mass.dat");
    oo.Write("M",example.M);
}
//#####################################################################
namespace PhysBAM{
template class MPM_KKT_DRIVER<VECTOR<float,1> >;
template class MPM_KKT_DRIVER<VECTOR<float,2> >;
template class MPM_KKT_DRIVER<VECTOR<float,3> >;
template class MPM_KKT_DRIVER<VECTOR<double,1> >;
template class MPM_KKT_DRIVER<VECTOR<double,2> >;
template class MPM_KKT_DRIVER<VECTOR<double,3> >;
}
