//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/GMRES.h>
#include <Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Log/DEBUG_SUBSTEPS.h>
#include <Tools/Log/LOG.h>
#include <Tools/Log/SCOPE.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_DRIVER.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <Hybrid_Methods/System/FLUID_KRYLOV_SYSTEM.h>
#include <Hybrid_Methods/System/FLUID_KRYLOV_VECTOR.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_VECTOR.h>
#include <Hybrid_Methods/System/MPM_OBJECTIVE.h>
#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Dynamics/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
#include <boost/function.hpp>
#ifdef USE_OPENMP
#include <omp.h>
#endif
using namespace PhysBAM;
namespace{
    template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
    {
        ((MPM_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
    }
};
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_DRIVER<TV>::
MPM_DRIVER(MPM_EXAMPLE<TV>& example)
    :example(example),objective(*new MPM_OBJECTIVE<TV>(example)),
    dv(*new MPM_KRYLOV_VECTOR<TV>(example.valid_grid_indices)),
    rhs(*new MPM_KRYLOV_VECTOR<TV>(example.valid_grid_indices)),
    fluid_sys(*new FLUID_KRYLOV_SYSTEM<TV>(example)),
    fluid_p(*new FLUID_KRYLOV_VECTOR<TV>(example.valid_pressure_indices)),
    fluid_rhs(*new FLUID_KRYLOV_VECTOR<TV>(example.valid_pressure_indices))
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_DRIVER<TV>::
~MPM_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
    delete &objective;
    delete &dv;
    delete &rhs;
    delete &fluid_sys;
    delete &fluid_p;
    delete &fluid_rhs;
    av.Delete_Pointers_And_Clean_Memory();
    bv.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Execute_Main_Program
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
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

    example.mass.Resize(example.grid.Domain_Indices(example.ghost));
    if(example.use_max_weight) example.max_weight.Resize(example.grid.Domain_Indices(example.ghost));
    example.volume.Resize(example.grid.Domain_Indices(example.ghost));
    example.velocity.Resize(example.grid.Domain_Indices(example.ghost));
    example.velocity_new.Resize(example.grid.Domain_Indices(example.ghost));
    example.velocity_check.Resize(example.grid.Domain_Indices(example.ghost));
    example.cell_C.Resize(example.grid.Domain_Indices(example.ghost));
    dv.u.Resize(example.grid.Domain_Indices(example.ghost));
    rhs.u.Resize(example.grid.Domain_Indices(example.ghost));
    objective.system.tmp.u.Resize(example.grid.Domain_Indices(example.ghost));

    if(example.use_fluid){
        fluid_p.p.Resize(example.grid.Domain_Indices(example.ghost));
        fluid_rhs.p.Resize(example.grid.Domain_Indices(example.ghost));
        example.mass_f.Resize(example.grid.Domain_Indices(example.ghost));
        example.volume_f.Resize(example.grid.Domain_Indices(example.ghost));
        example.density_f.Resize(example.grid.Domain_Indices(example.ghost));
        example.velocity_f.Resize(example.grid.Domain_Indices(example.ghost));
        example.velocity_new_f.Resize(example.grid.Domain_Indices(example.ghost));
        example.velocity_check_f.Resize(example.grid.Domain_Indices(example.ghost));
        example.cell_solid.Resize(example.grid.Domain_Indices(example.ghost));
        example.cell_pressure.Resize(example.grid.Domain_Indices(example.ghost));
        example.cell_C.Resize(example.grid.Domain_Indices(example.ghost));
        for(CELL_ITERATOR<TV> iterator(example.grid,example.ghost,GRID<TV>::GHOST_REGION);iterator.Valid();iterator.Next())
            example.cell_solid(iterator.Cell_Index())=true;}

    RANGE<TV_INT> range(example.grid.Cell_Indices(example.ghost));
    example.location.Resize(range,false,false);
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
            example.location.array(i++)=example.grid.Center(it.index);}

    if(!example.restart) Write_Output_Files(0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after init",0,1);
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Advance_One_Time_Step()
{
    example.Begin_Time_Step(example.time);

    Update_Simulated_Particles();
    Print_Particle_Stats("particle state",example.dt);
    Update_Particle_Weights();
    example.gather_scatter.Prepare_Scatter(example.particles);
    Particle_To_Grid();
    Print_Grid_Stats("after particle to grid",example.dt,example.velocity,0);
    Print_Energy_Stats("after particle to grid",example.velocity);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after particle to grid",0,1);
    Apply_Forces();
    Print_Grid_Stats("after forces",example.dt,example.velocity_new,&example.velocity);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after forces",0,1);
    if(example.use_fluid){
        // Energy and Momentum test
        T energyP=0,energyC=0,energyF=0,energyCAfter=0,energyPAfter=0,difVP=0,dif_check_old=0;
        TV momentumP,momentumC,momentumF;
        for(int k=0;k<example.simulated_particles.m;k++){
            int p=example.simulated_particles(k);
            energyP+=(T).5*example.particles.V(p).Magnitude_Squared()*example.particles.mass(p);
            momentumP+=example.particles.mass(p)*example.particles.V(p);}
        for(int i=0;i<example.valid_pressure_indices.m;i++){
            int j=example.valid_pressure_indices(i);
            energyC+=(T).5*example.mass.array(j)*example.velocity_new.array(j).Magnitude_Squared();
            momentumC+=example.mass.array(j)*example.velocity_new.array(j);}

        Compute_Cell_C();
        for(FACE_ITERATOR<TV> iterator(example.grid,example.ghost);iterator.Valid();iterator.Next()){
            FACE_INDEX<TV::m> face=iterator.Full_Index();
            energyF+=(T).5*example.mass_f(face)*example.velocity_f(face)*example.velocity_f(face);
            momentumF(face.axis)+=example.mass_f(face)*example.velocity_f(face);}
        //Cell_To_Face();
        Cell_To_Face_C();
        example.velocity_new_f=example.velocity_f;
        div=Pressure_Projection();
        if(example.use_f2p)
            Face_To_Particle();
        else{
            Face_To_Cell();
        for(int i=0;i<example.valid_pressure_indices.m;i++){
            int j=example.valid_pressure_indices(i);
            dif_check_old+=(example.velocity.array(j)-example.velocity_check.array(j)).Magnitude();}

            Grid_To_Particle();}
        for(int i=0;i<example.valid_pressure_indices.m;i++){
            int j=example.valid_pressure_indices(i);
            energyCAfter+=(T).5*example.mass.array(j)*example.velocity_check.array(j).Magnitude_Squared();}
        for(int k=0;k<example.simulated_particles.m;k++){
            int p=example.simulated_particles(k);
            energyPAfter+=(T).5*example.particles.V_check(p).Magnitude_Squared()*example.particles.mass(p);
            difVP+=(T)(example.particles.V_check(p)-example.particles.V(p)).Magnitude_Squared();}

        //LOG::printf("ENERGY [PARTICLE]:\t%e\n",energyP);
        //LOG::printf("ENERGY [CELL CENTER]:\t%e\n",energyC);
        //LOG::printf("ENERGY [FACE]:\t\t%e\n",energyF);
        //LOG::printf("ENERGY DIFF [CELL CENTER]:\t%e\n",energyCAfter-energyC);
        //LOG::printf("ENERGY DIFF [PARTICLE]:\t%e\n",energyPAfter-energyP);
        //LOG::printf("VELOCITY DIFF [PARTICLE]:\t%e\n",difVP);
        //LOG::printf("VELOCITY DIFF [CELL]:\t%e\n",dif_check_old);
        //LOG::printf("ENERGY RATIO [CELL CENTER/PARTICLE]:\t%e\n",energyC/energyP);
        //LOG::printf("ENERGY RATIO [FACE/CELL CENTER]:\t%e\n",energyF/energyC);
        //LOG::printf("ENERGY RATIO [FACE/PARTICLE]:\t\t%e\n",energyF/energyP);
        //LOG::printf("MOMENTUM [PARTICLE]:\t%P\n",momentumP);
        //LOG::printf("MOMENTUM [CELL CENTER]:\t%P\n",momentumC);
        //LOG::printf("MOMENTUM [FACE]:\t\t%P\n",momentumF);
    }
    else
        Grid_To_Particle();
    //TV v=(example.particles.X(0)-TV::Constant_Vector(0.5)).Normalized()*0.5;
    //std::swap(v(0),v(1)); v(0)*=-1;
    //example.particles.V(0)=v;
    //TV_INT index=example.grid.Index(example.particles.X(0));
    //LOG::printf("DIFFERENCE IN C: %e\n",(example.particles.C(0)-example.cell_C(index)).Frobenius_Norm_Squared());
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after grid to particle",0,1);
    /*if(div>1e-1){
        Write_Output_Files(++output_number);
        PHYSBAM_FATAL_ERROR("Div kills it");}*/

    example.End_Time_Step(example.time);
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
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
template<class TV> void MPM_DRIVER<TV>::
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
template<class TV> void MPM_DRIVER<TV>::
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
template<class TV> void MPM_DRIVER<TV>::
Update_Particle_Weights()
{
    example.weights->Update(example.particles.X);
    for(int i=0;i<TV::m;++i)
        example.face_weights(i)->Update(example.particles.X);
}
//#####################################################################
// Function Particle_To_Grid
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Particle_To_Grid()
{
    MPM_PARTICLES<TV>& particles=example.particles;

#pragma omp parallel for
    for(int i=0;i<example.mass.array.m;i++){
        example.mass.array(i)=0;
        example.volume.array(i)=0;
        example.velocity.array(i)=TV();
        example.velocity_new.array(i)=TV();
<<<<<<< HEAD
        if(example.use_max_weight) example.max_weight.array(i)=0;}
=======
        example.velocity_check.array(i)=TV();}
>>>>>>> density pressure projection working

    if(example.weights->use_gradient_transfer)
    {
        example.gather_scatter.template Scatter<int>(
            [this,&particles](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
            {
                example.mass(it.Index())+=it.Weight()*particles.mass(p);
<<<<<<< HEAD
                if(example.use_max_weight){T& w=example.max_weight(it.Index());w=max(w,it.Weight());}
=======
                example.volume(it.Index())+=it.Weight()*particles.volume(p);
>>>>>>> density pressure projection working
                TV V=particles.V(p);
                if(example.use_affine){
                    V+=particles.C(p)*((example.grid.Center(it.Index())-particles.X(p)));
                    example.velocity(it.Index())+=it.Weight()*particles.mass(p)*V;}
                else example.velocity(it.Index())+=particles.mass(p)*(it.Weight()*particles.V(p)+particles.B(p)*it.Gradient());
            },true);
    }
    else if(example.weights->constant_scalar_inertia_tensor)
    {
        T Dp_inverse=example.weights->Constant_Scalar_Inverse_Dp();
        example.gather_scatter.template Scatter<int>(
            [this,Dp_inverse,&particles](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
            {
                example.mass(it.Index())+=it.Weight()*particles.mass(p);
<<<<<<< HEAD
                if(example.use_max_weight) example.max_weight(it.Index())+=it.Weight();
=======
                example.volume(it.Index())+=it.Weight()*particles.volume(p);
>>>>>>> density pressure projection working
                example.cell_C(it.Index())+=MATRIX<T,TV::m>::Outer_Product(particles.V(p),it.Gradient());
                TV V=particles.V(p);
                if(example.use_affine) {
                    if(example.use_fluid)
                        V+=particles.C(p)*(example.grid.Center(it.Index())-particles.X(p));
                    else
                        V+=particles.B(p)*(Dp_inverse*(example.grid.Center(it.Index())-particles.X(p)));}
                example.velocity(it.Index())+=it.Weight()*particles.mass(p)*V;
            },true);
    }
    else PHYSBAM_FATAL_ERROR("General case for rasterization not implemented");

    example.valid_grid_indices.Remove_All();
    example.valid_grid_cell_indices.Remove_All();
    if(example.use_fluid){
        example.valid_pressure_indices.Remove_All();
        example.valid_pressure_cell_indices.Remove_All();
        example.cell_pressure.Fill(false);}
    for(RANGE_ITERATOR<TV::m> it(example.mass.domain);it.Valid();it.Next()){
        int i=example.mass.Standard_Index(it.index);
        if(example.mass.array(i)){
            if(example.use_fluid){
                example.cell_pressure.array(i)=true;
                if(!example.cell_solid.array(i)){
                    bool valid=false;
                    if(example.weights->Order()==1)
                        valid=true;
                    else {
                        int c=0;
                        for(int d=0;d<TV::m;++d)
                            c+=(example.cell_solid(it.index-TV_INT::Axis_Vector(d))||example.cell_solid(it.index+TV_INT::Axis_Vector(d)));
                        if(c<TV::m) valid=true;}
                    if(valid){
                        example.valid_pressure_indices.Append(i);
                        example.valid_pressure_cell_indices.Append(it.index);}}}
            example.valid_grid_indices.Append(i);
            example.valid_grid_cell_indices.Append(it.index);
            example.velocity.array(i)/=example.mass.array(i);}
        else example.velocity.array(i)=TV();}
        //LOG::printf("These are the flags:\n%P\n",example.cell_pressure);
}
//#####################################################################
// Function Grid_To_Particle
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
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
            TV Vn_interpolate,V_pic,V_flip=particles.V(p),V_pic_check,V_flip_check=particles.V(p);
            MATRIX<T,TV::m> B,C,grad_Vp;

            for(PARTICLE_GRID_ITERATOR<TV> it(example.weights,p,true,scratch);it.Valid();it.Next()){
                if(output_number==87)
                LOG::cout<<"index="<<it.Index()<<std::endl;
                T w=it.Weight();
                TV_INT index=it.Index();
                TV V_grid=example.velocity_new(index);
                TV V_grid_check=example.velocity_check(index);
                V_pic+=w*V_grid;
                V_pic_check+=w*V_grid_check;
                V_flip+=w*(V_grid-example.velocity(index));
                V_flip_check+=w*(V_grid-example.velocity_check(index));
                Vn_interpolate+=w*example.velocity(index);
                if(example.use_midpoint)
                    V_grid=(T).5*(V_grid+example.velocity(index));
                grad_Vp+=MATRIX<T,TV::m>::Outer_Product(V_grid,it.Gradient());}
            particles.F(p)+=dt*grad_Vp*particles.F(p);
<<<<<<< HEAD
            if(particles.store_S) particles.S(p)+=dt*(grad_Vp*particles.S(p)).Twice_Symmetric_Part();

            if(example.use_affine && example.use_early_gradient_transfer)
                B=grad_Vp/example.weights->Constant_Scalar_Inverse_Dp();
            else if(example.use_affine)
=======
            if(example.use_affine)
>>>>>>> density pressure projection working
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
                            +MATRIX<T,TV::m>::Outer_Product(Z-particles.X(p)-xi_new+xp_new,V_grid));}

<<<<<<< HEAD
            particles.V(p)=V_pic;
            Perform_Particle_Collision(p,example.time+example.dt);
            if(example.use_midpoint) particles.X(p)+=(particles.V(p)+Vn_interpolate)*(dt/2);
            else particles.X(p)+=particles.V(p)*dt;
            particles.V(p)=V_flip*example.flip+V_pic*(1-example.flip);
            particles.B(p)=B;
            SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> Dp_inverse_transpose=example.weights->Dp(particles.X(p)).Inverse();
            particles.C(p)=B*Dp_inverse_transpose;
=======
                if(output_number==87)
                LOG::cout<<"passes"<<std::endl;
                    particles.V(p)=V_pic;
                    Perform_Particle_Collision(p);
                    SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> Dp_inverse_transpose=example.weights->Dp(particles.X(p)).Inverse();
                    if(example.use_midpoint) particles.X(p)+=(particles.V(p)+Vn_interpolate)*(dt/2);
                    else particles.X(p)+=particles.V(p)*dt;
                    particles.V(p)=V_flip*example.flip+V_pic*(1-example.flip);
                    particles.V_check(p)=V_flip_check*example.flip+V_pic_check*(1-example.flip);
                    particles.B(p)=B;
                    particles.C(p)=B*Dp_inverse_transpose;
>>>>>>> density pressure projection working

            if(!example.grid.domain.Lazy_Inside(particles.X(p))) particles.valid(p)=false;}}
}
//#####################################################################
// Function Face_To_Particle
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Face_To_Particle()
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
            MATRIX<T,TV::m> C;
            for(int i=0;i<TV::m;++i)
                for(PARTICLE_GRID_ITERATOR<TV> it(example.face_weights(i),p,true,scratch);it.Valid();it.Next()){
                    T w=it.Weight();
                    FACE_INDEX<TV::m> face_index(i,it.Index());
                    T V_grid=example.velocity_new_f(face_index);
                    V_pic(i)+=w*V_grid;
                    V_flip(i)+=w*(V_grid-example.velocity_f(face_index));
                    if(example.use_midpoint){
                        PHYSBAM_NOT_IMPLEMENTED("Midpoint with Face_To_Particle is not supported");}

                    if(example.use_affine){
                        TV tmp=example.velocity_new_f(face_index)*it.Gradient();
                        for(int k=0;k<TV::m;++k)
                            C(i,k)+=tmp(k);}}
            particles.V(p)=V_pic;
            Perform_Particle_Collision(p);
            if(example.use_midpoint) PHYSBAM_NOT_IMPLEMENTED("Midpoint with Face_To_Particle is not supported");
            else particles.X(p)+=particles.V(p)*dt;
            particles.V(p)=V_flip*example.flip+V_pic*(1-example.flip);
            particles.C(p)=C;
            particles.F(p)+=dt*C*particles.F(p);
            if(!example.grid.domain.Lazy_Inside(particles.X(p))) particles.valid(p)=false;}}
}
//#####################################################################
// Function Compute_Cell_C
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Compute_Cell_C()
{
    example.cell_C.Fill(MATRIX<T,TV::m>());
    T one_over_dx=example.grid.one_over_dX(0);
    for(int t=0;t<example.valid_grid_cell_indices.m;t++){
        TV_INT index=example.valid_grid_cell_indices(t);
        MATRIX<T,TV::m> C;
        for(int i=0;i<TV::m;i++)
            for(int j=0;j<TV::m;j++){
                TV_INT e=TV_INT::Axis_Vector(j);
                if(example.mass(index+e)) 
                    C(i,j)=(example.velocity_new(index+e)(i)-example.velocity_new(index)(i))*one_over_dx;
                else if(example.mass(index-e)) 
                    C(i,j)=(example.velocity_new(index)(i)-example.velocity_new(index-e)(i))*one_over_dx;
                else
                    PHYSBAM_FATAL_ERROR("Failed to compute cell centered C");}
        example.cell_C(index)=C;}
}
//#####################################################################
// Function Cell_To_Face
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Cell_To_Face()
{
    // TODO: parallelzie
    example.mass_f.Fill(0);
    example.velocity_f.Fill(0);
    example.velocity_new_f.Fill(0);

#pragma omp parallel for
    for(int t=0;t<example.valid_grid_cell_indices.m;t++){
        TV_INT index=example.valid_grid_cell_indices(t);
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> face(a,example.grid.First_Face_Index_In_Cell(a,index));
            example.mass_f(face)+=example.mass(index);
            example.velocity_f(face)+=example.mass(index)*example.velocity_new(index)(a);}}
#pragma omp parallel for
    for(int t=0;t<example.valid_grid_cell_indices.m;t++){
        TV_INT index=example.valid_grid_cell_indices(t);
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> face(a,example.grid.Second_Face_Index_In_Cell(a,index));
            example.mass_f(face)+=example.mass(index);
            example.velocity_f(face)+=example.mass(index)*example.velocity_new(index)(a);}}

    // TODO: improve efficiency
    for(FACE_ITERATOR<TV> iterator(example.grid,example.ghost);iterator.Valid();iterator.Next()){
        FACE_INDEX<TV::m> face=iterator.Full_Index();
        if(example.mass_f(face)){
            example.velocity_f(face)/=example.mass_f(face);
            example.mass_f(face)*=(T).5;}}
}
//#####################################################################
// Function Cell_To_Face_C
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Cell_To_Face_C()
{
    // TODO: parallelzie
    example.mass_f.Fill(0);
    example.volume_f.Fill(0);
    example.velocity_f.Fill(0);
    example.velocity_new_f.Fill(0);
    example.density_f.Fill(0);

#pragma omp parallel for
    for(int t=0;t<example.valid_grid_cell_indices.m;t++){
        TV_INT index=example.valid_grid_cell_indices(t);
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> face(a,example.grid.First_Face_Index_In_Cell(a,index));
            example.mass_f(face)+=example.mass(index);
            example.volume_f(face)+=example.volume(index);
            example.velocity_f(face)+=example.mass(index)*(example.velocity_new(index)(a)-(example.cell_C(index)*TV::Axis_Vector(a)*example.grid.dX(a))(a)*(T).5);}}
#pragma omp parallel for
    for(int t=0;t<example.valid_grid_cell_indices.m;t++){
        TV_INT index=example.valid_grid_cell_indices(t);
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> face(a,example.grid.Second_Face_Index_In_Cell(a,index));
            example.mass_f(face)+=example.mass(index);
            example.volume_f(face)+=example.volume(index);
            example.velocity_f(face)+=example.mass(index)*(example.velocity_new(index)(a)+(example.cell_C(index)*TV::Axis_Vector(a)*example.grid.dX(a))(a)*(T).5);}}

    // TODO: improve efficiency
    for(FACE_ITERATOR<TV> iterator(example.grid,example.ghost);iterator.Valid();iterator.Next()){
        FACE_INDEX<TV::m> face=iterator.Full_Index();
        if(example.mass_f(face)){
            example.velocity_f(face)/=example.mass_f(face);
            example.mass_f(face)*=(T).5;
            example.volume_f(face)*=(T).5;
            if(abs(example.volume_f(face))>1e-23)
                example.density_f(face)=example.mass_f(face)/example.volume_f(face);
            else{
                LOG::printf("volume_f[%P]=%P\n",face,example.volume_f(face));
                PHYSBAM_FATAL_ERROR("THERE IS A CELL WITH VOLUME <= 1E-16");}
            if(abs(example.volume_f(face))<1e-10)
                LOG::printf("density_f[%P]=%f\n",face,example.density_f(face));
            if(abs(example.density_f(face)-(T)2)>1e-10)
                PHYSBAM_FATAL_ERROR("DENSITY IS WRONG!\n");
        }}
}
//#####################################################################
// Function Face_To_Cell
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Face_To_Cell()
{
#pragma omp parallel for
    for(int t=0;t<example.valid_grid_cell_indices.m;t++){
        TV_INT index=example.valid_grid_cell_indices(t);
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> faceF(a,example.grid.First_Face_Index_In_Cell(a,index));
            FACE_INDEX<TV::m> faceS(a,example.grid.Second_Face_Index_In_Cell(a,index));
            example.velocity_new(index)(a)=(example.velocity_new_f(faceF)+example.velocity_new_f(faceS))*(T).5;
            example.velocity_check(index)(a)=(example.velocity_check_f(faceF)+example.velocity_check_f(faceS))*(T).5; 
        }}
}
//#####################################################################
// Function Face_To_Cell_FLIP
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Face_To_Cell_FLIP()
{
    T flip_ratio=(T)0.75;
#pragma omp parallel for
    for(int t=0;t<example.valid_grid_cell_indices.m;t++){
        TV_INT index=example.valid_grid_cell_indices(t);
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> faceF(a,example.grid.First_Face_Index_In_Cell(a,index));
            FACE_INDEX<TV::m> faceS(a,example.grid.Second_Face_Index_In_Cell(a,index));
            T pic=(example.velocity_new_f(faceF)+example.velocity_new_f(faceS))*(T).5;
            T flip=example.velocity_new(index)(a)+pic-
                (example.velocity_f(faceF)+example.velocity_f(faceS))*(T).5;
            example.velocity_new(index)(a)=flip_ratio*flip+(1-flip_ratio)*pic;
        }}
}
//#####################################################################
// Function Pressure_Projection
//#####################################################################
//template<class TV> void MPM_DRIVER<TV>::
template<class TV> typename TV::SCALAR MPM_DRIVER<TV>::
Pressure_Projection()
{
    CONJUGATE_GRADIENT<T> cg; 
    T one_over_dx=example.grid.one_over_dX(0);
    fluid_p.p.Fill(0);
    fluid_rhs.p.Fill(0);
    for(int t=0;t<example.valid_pressure_cell_indices.m;t++){
        TV_INT index=example.valid_pressure_cell_indices(t);
        for(int a=0;a<TV::m;a++){
            const TV_INT axis=TV_INT::Axis_Vector(a);
            FACE_INDEX<TV::m> faceF(a,example.grid.First_Face_Index_In_Cell(a,index));
            FACE_INDEX<TV::m> faceS(a,example.grid.Second_Face_Index_In_Cell(a,index));
            fluid_rhs.p(index)+=-(example.velocity_f(faceS)-example.velocity_f(faceF));
            if(example.weights->Order()==1){
                if(example.cell_solid(index-axis)&&example.cell_solid(index+axis))
                    PHYSBAM_NOT_IMPLEMENTED("Fluid trapped between two solids (rhs)");
                else if(example.cell_solid(index-axis)&&!example.cell_solid(index+axis))
                    fluid_rhs.p(index)-=(example.velocity_f(faceS)+example.velocity_f(faceF));
                else if(!example.cell_solid(index-axis)&&example.cell_solid(index+axis))
                    fluid_rhs.p(index)+=(example.velocity_f(faceS)+example.velocity_f(faceF));}
            else PHYSBAM_NOT_IMPLEMENTED();
        }
        fluid_rhs.p(index)*=one_over_dx;}
    cg.finish_before_indefiniteness=true;
    cg.relative_tolerance=false;
    bool converged=cg.Solve(fluid_sys,fluid_p,fluid_rhs,bv,1e-8,0,1000);
    if(!converged)
        LOG::printf("CG DID NOT CONVERGE.\n");
    example.velocity_new_f=example.velocity_f;
    example.velocity_check_f=example.velocity_f;
    //UPDATE VELOCITY ON FACES
    for(int t=0;t<example.valid_pressure_cell_indices.m;t++){
        TV_INT index=example.valid_pressure_cell_indices(t);
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> faceF(a,example.grid.First_Face_Index_In_Cell(a,index));
            FACE_INDEX<TV::m> faceS(a,example.grid.Second_Face_Index_In_Cell(a,index));
            example.velocity_new_f(faceF)-=fluid_p.p(index)*one_over_dx/(T)example.density_f(faceF);
            example.velocity_new_f(faceS)+=fluid_p.p(index)*one_over_dx/(T)example.density_f(faceS);}}
    /*for(int t=0;t<example.valid_pressure_cell_indices.m;t++){
        TV_INT index=example.valid_pressure_cell_indices(t);
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> faceF(a,example.grid.First_Face_Index_In_Cell(a,index));
            FACE_INDEX<TV::m> faceS(a,example.grid.Second_Face_Index_In_Cell(a,index));
            example.velocity_new_f(faceF)-=fluid_p.p(index)*one_over_dx;
            example.velocity_new_f(faceS)+=fluid_p.p(index)*one_over_dx;}}*/
    //APPLY BC ON THE FACES
    for(int t=0;t<example.valid_pressure_cell_indices.m;t++){
        TV_INT index=example.valid_pressure_cell_indices(t);
        for(int a=0;a<TV::m;a++){
            const TV_INT axis=TV_INT::Axis_Vector(a);
            FACE_INDEX<TV::m> faceF(a,example.grid.First_Face_Index_In_Cell(a,index));
            FACE_INDEX<TV::m> faceS(a,example.grid.Second_Face_Index_In_Cell(a,index));
            if(example.cell_solid(index+axis))
                example.velocity_new_f(faceS)=-example.velocity_new_f(faceF);
            if(example.cell_solid(index-axis))
                example.velocity_new_f(faceF)=-example.velocity_new_f(faceS);}}

    // Diagnostics
    T divmax=0;
    TV_INT where;
    int cell_type;
    for(int t=0;t<example.valid_pressure_cell_indices.m;t++){
        TV_INT index=example.valid_pressure_cell_indices(t);
        T div=0;
        for(int a=0;a<TV::m;a++){
            const TV_INT axis=TV_INT::Axis_Vector(a);
            FACE_INDEX<TV::m> faceF(a,example.grid.First_Face_Index_In_Cell(a,index));
            FACE_INDEX<TV::m> faceS(a,example.grid.Second_Face_Index_In_Cell(a,index));
            if(!example.cell_solid(index-axis)||!example.cell_solid(index+axis)||!example.cell_solid(index-2*axis)||!example.cell_solid(index+2*axis))
                div+=example.velocity_new_f(faceS)-example.velocity_new_f(faceF);}
        if(abs(div)>divmax) {divmax=abs(div); where=index; cell_type=example.cell_pressure(index);}}
    divmax*=one_over_dx;

    LOG::cout<<"Divergence after projection\t"<<divmax<<std::endl;
    LOG::printf("where: %P\t cell_type: %d\n",where,cell_type);
    return divmax;
}
//#####################################################################
// Function Apply_Forces
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Apply_Forces()
{
    example.Capture_Stress();
    objective.Reset();
    if(example.use_symplectic_euler){
        objective.tmp2*=0;
        example.Precompute_Forces(example.time,false);
        example.Add_Forces(objective.tmp2.u,example.time);
#pragma omp parallel for
        for(int i=0;i<example.valid_grid_indices.m;i++){
            int p=example.valid_grid_indices(i);
            dv.u.array(p)=example.dt/example.mass.array(p)*objective.tmp2.u.array(p);}}
    else{
        NEWTONS_METHOD<T> newtons_method;
        newtons_method.tolerance=example.newton_tolerance*example.dt;
        newtons_method.progress_tolerance=1e-5;
        newtons_method.max_iterations=example.newton_iterations;
        newtons_method.krylov_tolerance=example.solver_tolerance;
        newtons_method.max_krylov_iterations=example.solver_iterations;
        newtons_method.use_cg=true;
        newtons_method.debug=true;

        newtons_method.require_one_iteration=!objective.Initial_Guess(dv,newtons_method.tolerance);
        LOG::printf("max velocity: %P\n",Max_Particle_Speed());
        if(example.test_diff) objective.Test_Diff(dv);

        objective.system.forced_collisions.Remove_All();
        bool converged=newtons_method.Newtons_Method(objective,objective.system,dv,av);
        if(!converged) LOG::cout<<"WARNING: Newton's method did not converge"<<std::endl;}

    Apply_Friction();
    objective.Restore_F();

#pragma omp parallel for
    for(int i=0;i<example.valid_grid_indices.m;i++){
        int j=example.valid_grid_indices(i);
        example.velocity_new.array(j)=dv.u.array(j)+objective.v0.u.array(j);}
    example.velocity_new.array.Subset(objective.system.stuck_nodes)=objective.system.stuck_velocity;
}
//#####################################################################
// Function Perform_Particle_Collision
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
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
template<class TV> void MPM_DRIVER<TV>::
Apply_Friction()
{
    if(!example.collision_objects.m) return;
    objective.Adjust_For_Collision(dv);
    objective.Compute_Unconstrained(dv,0,&objective.tmp0,0);
    objective.tmp1=objective.tmp0;
    objective.Project_Gradient_And_Prune_Constraints(objective.tmp1,true);

    objective.v1.u.array.Subset(objective.system.stuck_nodes).Fill(TV());
    for(int i=0;i<objective.system.collisions.m;i++){
        const typename MPM_KRYLOV_SYSTEM<TV>::COLLISION& c=objective.system.collisions(i);
        TV& v=objective.v1.u.array(c.p);
        T normal_force=TV::Dot_Product(c.n,objective.tmp0.u.array(c.p)-objective.tmp1.u.array(c.p));
        TV t=v.Projected_Orthogonal_To_Unit_Direction(c.n);
        T t_mag=t.Normalize();
        T coefficient_of_friction=example.collision_objects(c.object)->friction;
        T k=coefficient_of_friction*normal_force/example.mass.array(c.p);
        if(t_mag<=k)
            v.Project_On_Unit_Direction(c.n);
        else v-=k*t;
        dv.u.array(c.p)=v-objective.v0.u.array(c.p);}
}
//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> typename TV::SCALAR MPM_DRIVER<TV>::
Compute_Dt() const
{
    T critical_speed=example.cfl*example.grid.DX().Min()/example.max_dt;
    T v=Grid_V_Upper_Bound();
    return (v>critical_speed)?(example.cfl*example.grid.DX().Min()/v):example.max_dt;
}
//#####################################################################
// Function Max_Particle_Speed
//#####################################################################
template<class TV> typename TV::SCALAR MPM_DRIVER<TV>::
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
template<class TV> typename TV::SCALAR MPM_DRIVER<TV>::
Grid_V_Upper_Bound() const
{
    if(!example.use_affine || !example.weights->constant_scalar_inertia_tensor) return Max_Particle_Speed();
    T result=0;
    T xi=(T)6*sqrt((T)TV::m)*example.grid.One_Over_DX().Min();
#pragma omp parallel for reduction(max:result)
    for(int k=0;k<example.simulated_particles.m;k++){
        int p=example.simulated_particles(k);
        result=max(result,example.particles.V(p).Magnitude()+example.particles.B(p).Frobenius_Norm()*xi);}
    return result;
}
//#####################################################################
// Function Update_Simulated_Particles
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
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

    for(int i=0;i<example.lagrangian_forces.m;i++){
        example.lagrangian_forces(i)->use_implicit_velocity_independent_forces=true;
        example.lagrangian_forces(i)->Update_Mpi(example.particle_is_simulated,0);}
}
//#####################################################################
// Function Print_Grid_Stats
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
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
template<class TV> void MPM_DRIVER<TV>::
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
template<class TV> void MPM_DRIVER<TV>::
Print_Energy_Stats(const char* str,const ARRAY<TV,TV_INT>& u)
{
    if(!example.print_stats) return;
    example.Capture_Stress();
    example.Precompute_Forces(example.time,false);
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
namespace PhysBAM{
template class MPM_DRIVER<VECTOR<float,2> >;
template class MPM_DRIVER<VECTOR<float,3> >;
template class MPM_DRIVER<VECTOR<double,2> >;
template class MPM_DRIVER<VECTOR<double,3> >;
}

