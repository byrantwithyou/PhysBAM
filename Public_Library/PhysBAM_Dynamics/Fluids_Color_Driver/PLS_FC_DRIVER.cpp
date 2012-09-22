//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/QUADRATIC_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <PhysBAM_Tools/Krylov_Solvers/MINRES.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Vectors/VECTOR_UTILITIES.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_VECTOR_COLOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER.h>
#include <PhysBAM_Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER_POLY.h>
#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_PHI_WATER.h>
#include <PhysBAM_Dynamics/Fluids_Color_Driver/PLS_FC_DRIVER.h>
#include <PhysBAM_Dynamics/Fluids_Color_Driver/PLS_FC_EXAMPLE.h>
using namespace PhysBAM;
namespace{
template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
{
    ((PLS_FC_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
}
};
//#####################################################################
// Initialize
//#####################################################################
template<class TV> PLS_FC_DRIVER<TV>::
PLS_FC_DRIVER(PLS_FC_EXAMPLE<TV>& example)
    :example(example)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> PLS_FC_DRIVER<TV>::
~PLS_FC_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Initialize()
{
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);

    // setup time
    current_frame=example.restart;
    output_number=current_frame;
    time=example.time_steps_per_frame*current_frame*example.dt;

    example.levelset_color.phi.Resize(example.grid.Node_Indices(1));
    example.levelset_color.color.Resize(example.grid.Node_Indices(1));
    example.face_color.Resize(example.grid,3);
    example.prev_face_color.Resize(example.grid,3);
    example.face_velocities.Resize(example.number_of_colors);
    example.prev_face_velocities.Resize(example.number_of_colors);
    for(int i=0;i<example.number_of_colors;i++){
        example.face_velocities(i).Resize(example.grid,3);
        example.prev_face_velocities(i).Resize(example.grid,3);}

    {
        example.particle_levelset_evolution.Initialize_Domain(example.grid);
        example.particle_levelset_evolution.particle_levelset.Set_Band_Width(6);
        example.collision_bodies_affecting_fluid.Initialize_Grids();
    }

    example.particle_levelset_evolution.time=time;
    example.particle_levelset_evolution.Set_CFL_Number((T).9);
    example.particle_levelset_evolution.Use_Reinitialization();

    example.boundary=&example.boundary_scalar;
    example.phi_boundary=&example.cell_extrapolate;

    VECTOR<VECTOR<bool,2>,TV::dimension> domain_open_boundaries=VECTOR_UTILITIES::Complement(example.domain_boundary);
    example.phi_boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    example.boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    example.particle_levelset_evolution.Levelset_Advection(0).Set_Custom_Advection(example.advection_scalar);

    {
        example.particle_levelset_evolution.Initialize_Domain(example.grid);
        example.particle_levelset_evolution.particle_levelset.Set_Band_Width(6);
        example.collision_bodies_affecting_fluid.Initialize_Grids();
    }
    
    example.particle_levelset_evolution.Particle_Levelset(0).number_of_ghost_cells=example.number_of_ghost_cells;
    example.particle_levelset_evolution.Set_Number_Particles_Per_Cell(16);
    example.particle_levelset_evolution.Set_Levelset_Callbacks(example);
    example.particle_levelset_evolution.Initialize_FMM_Initialization_Iterative_Solver(true);

    example.particle_levelset_evolution.particle_levelset.levelset.Set_Custom_Boundary(*example.phi_boundary);
    example.particle_levelset_evolution.Bias_Towards_Negative_Particles(false);
    example.particle_levelset_evolution.Particle_Levelset(0).Store_Unique_Particle_Id();
    example.particle_levelset_evolution.use_particle_levelset=true;
    example.particle_levelset_evolution.particle_levelset.levelset.Set_Collision_Body_List(example.collision_bodies_affecting_fluid);
// TODO    example.particle_levelset_evolution.particle_levelset.levelset.Set_Face_Velocities_Valid_Mask(&example.incompressible.valid_mask);
    example.particle_levelset_evolution.particle_levelset.Set_Collision_Distance_Factors(.1,1);

    {
        example.particle_levelset_evolution.Initialize_Domain(example.grid);
        example.particle_levelset_evolution.particle_levelset.Set_Band_Width(6);
        example.collision_bodies_affecting_fluid.Initialize_Grids();
    }

    if(example.restart){
        example.Read_Output_Files(example.restart);
        example.collision_bodies_affecting_fluid.Rasterize_Objects();
        example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*example.grid.Minimum_Edge_Length(),5);} // compute grid visibility (for advection later)
    else{
        example.collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false);
        example.collision_bodies_affecting_fluid.Rasterize_Objects();
        example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*example.grid.Minimum_Edge_Length(),5);
        example.Initialize();
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(example.grid,3);it.Valid();it.Next())
            example.face_color(it.Full_Index())=example.levelset_color.Color(it.Location());
        example.particle_levelset_evolution.Make_Signed_Distance();
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time);}

    example.collision_bodies_affecting_fluid.Compute_Grid_Visibility();
    example.particle_levelset_evolution.Set_Seed(2606);
//    Put these back in for PLS
//    if(!example.restart) example.particle_levelset_evolution.Seed_Particles(time);
//    example.particle_levelset_evolution.Delete_Particles_Outside_Grid();

    int extrapolation_cells=2*example.number_of_ghost_cells+2;
    ARRAY<T,TV_INT> exchanged_phi_ghost(example.grid.Domain_Indices(extrapolation_cells));
    example.particle_levelset_evolution.particle_levelset.levelset.boundary->Fill_Ghost_Cells(example.grid,example.particle_levelset_evolution.phi,exchanged_phi_ghost,0,time,extrapolation_cells);
// TODO    example.incompressible.Extrapolate_Velocity_Across_Interface(example.face_velocities,exchanged_phi_ghost,false,example.number_of_ghost_cells,0,TV());

    if(!example.restart) Write_Output_Files(0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after init",0,1);
}
//#####################################################################
// Function Update_Pls
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Update_Pls(T dt)
{
    T maximum_fluid_speed=0;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(example.grid);it.Valid();it.Next()){
        int c=example.face_color(it.Full_Index());
        if(c<0) continue;
        maximum_fluid_speed=max(maximum_fluid_speed,abs(example.face_velocities(c)(it.Full_Index())));}

    T max_particle_collision_distance=example.particle_levelset_evolution.Particle_Levelset(0).max_collision_distance_factor*example.grid.dX.Max();
    example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(true,dt*maximum_fluid_speed+2*max_particle_collision_distance+(T).5*example.grid.dX.Max(),10);

    LOG::Time("Fill ghost");
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before advection",0,1);
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities_ghost;
    face_velocities_ghost.Resize(example.grid,example.number_of_ghost_cells,false);
    // TODO example.incompressible.boundary->Fill_Ghost_Cells_Face(example.grid,example.face_velocities,face_velocities_ghost,time+dt,example.number_of_ghost_cells);

    ARRAY<T,TV_INT> phi_back(example.grid.Domain_Indices(example.number_of_ghost_cells));
    example.phi_boundary->Fill_Ghost_Cells(example.grid,example.particle_levelset_evolution.Particle_Levelset(0).levelset.phi,phi_back,dt,time,example.number_of_ghost_cells);
    LOG::Time("Advect Levelset");
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before phi",0,1);
    example.particle_levelset_evolution.Advance_Levelset(dt);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after phi",0,1);
    LOG::Time("advecting particles");
    example.particle_levelset_evolution.Particle_Levelset(0).Euler_Step_Particles(face_velocities_ghost,dt,time,true,true,false,false);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after advection",0,1);

    LOG::Time("updating removed particle velocities");
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after particles",0,1);

    example.particle_levelset_evolution.Make_Signed_Distance();
    // TODO example.boundary->Fill_Ghost_Cells_Face(example.grid,example.face_velocities,face_velocities_ghost,time+dt,example.number_of_ghost_cells);

    LOG::Time("modifying levelset");
    example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time+dt);
    example.particle_levelset_evolution.Particle_Levelset(0).Exchange_Overlap_Particles();
    example.particle_levelset_evolution.Modify_Levelset_And_Particles(&face_velocities_ghost);
    example.particle_levelset_evolution.Make_Signed_Distance();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after particles",0,1);

    LOG::Time("deleting particles"); // needs to be after adding sources, since it does a reseed
    example.particle_levelset_evolution.Delete_Particles_Outside_Grid();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete 1",0,1);
    LOG::Time("deleting particles in local maxima");
    example.particle_levelset_evolution.particle_levelset.Delete_Particles_In_Local_Maximum_Phi_Cells(0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete 2",0,1);
    LOG::Time("deleting particles far from interface");
    example.particle_levelset_evolution.Particle_Levelset(0).Delete_Particles_Far_From_Interface(); // uses visibility
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete 3",0,1);

    LOG::Time("re-incorporating removed particles");
    example.particle_levelset_evolution.Particle_Levelset(0).Identify_And_Remove_Escaped_Particles(face_velocities_ghost,1.5,time+dt);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after remove",0,1);
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Advance_One_Time_Step(bool first_step)
{
    Extrapolate_Velocity(example.face_velocities,example.face_color);
    example.Begin_Time_Step(time);
    T dt=example.dt;
//    Update_Pls(dt);

    Advection_And_BDF(dt,first_step);
    Apply_Pressure_And_Viscosity(dt,first_step);
    time+=dt;
    Extrapolate_Velocity(example.face_velocities,example.face_color);
    example.End_Time_Step(time);
}
//#####################################################################
// Function Advection_And_BDF
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Advection_And_BDF(T dt,bool first_step)
{
    Extrapolate_Velocity(example.face_velocities,example.face_color);
    Extrapolate_Velocity(example.prev_face_velocities,example.prev_face_color);
    if(!example.use_advection)
        for(int c=0;c<example.number_of_colors;c++)
            No_Advection_And_BDF(dt,first_step,c);
    else if(example.use_reduced_advection)
        for(int c=0;c<example.number_of_colors;c++)
            Reduced_Advection_And_BDF(dt,first_step,c);
    else
        for(int c=0;c<example.number_of_colors;c++)
            RK2_Advection_And_BDF(dt,first_step,c);
    Swap_State();
    Extrapolate_Velocity(example.face_velocities,example.face_color);
}
//#####################################################################
// Function Advection_And_BDF
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
No_Advection_And_BDF(T dt,bool first_step,int c)
{
    if(!first_step) example.prev_face_velocities(c).Copy((T)2/(T)1.5,example.face_velocities(c),-(T).5/(T)1.5,example.prev_face_velocities(c));
    else example.prev_face_velocities(c)=example.face_velocities(c);
}
//#####################################################################
// Function Advection_And_BDF
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Reduced_Advection_And_BDF(T dt,bool first_step,int c)
{
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T,AVERAGING_UNIFORM<GRID<TV> >,QUADRATIC_INTERPOLATION_UNIFORM<GRID<TV>,T> > quadratic_advection;
    BOUNDARY_UNIFORM<GRID<TV>,T> boundary;
    FACE_LOOKUP_UNIFORM<GRID<TV> > lookup_face_velocities(example.face_velocities(c)),lookup_prev_face_velocities(example.prev_face_velocities(c));
    if(!first_step){
        ARRAY<T,FACE_INDEX<TV::dimension> > temp(example.grid,3);
        quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,temp,lookup_prev_face_velocities,lookup_prev_face_velocities,boundary,2*dt,time+dt);
        quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,example.prev_face_velocities(c),lookup_face_velocities,lookup_face_velocities,boundary,dt,time+dt);
        example.prev_face_velocities(c).Copy((T)2/(T)1.5,example.prev_face_velocities(c),-(T).5/(T)1.5,temp);}
    else{
        quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,example.prev_face_velocities(c),lookup_face_velocities,lookup_face_velocities,boundary,dt,time+dt);}
}
//#####################################################################
// Function Advection_And_BDF
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
RK2_Advection_And_BDF(T dt,bool first_step,int c)
{
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T,AVERAGING_UNIFORM<GRID<TV> >,QUADRATIC_INTERPOLATION_UNIFORM<GRID<TV>,T> > quadratic_advection;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T,AVERAGING_UNIFORM<GRID<TV> >,LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> > linear_advection;
    BOUNDARY_UNIFORM<GRID<TV>,T> boundary;
    ARRAY<T,FACE_INDEX<TV::dimension> > temp(example.grid,3),temp2(example.grid,3),temp3(example.grid,3),temp4(example.grid,3),temp5(example.grid,3);
    FACE_LOOKUP_UNIFORM<GRID<TV> > lookup_temp(temp),lookup_temp2(temp2),lookup_temp3(temp3),lookup_temp4(temp4),lookup_temp5(temp5);
    FACE_LOOKUP_UNIFORM<GRID<TV> > lookup_face_velocities(example.face_velocities(c)),lookup_prev_face_velocities(example.prev_face_velocities(c));
    example.prev_face_color=example.face_color;
    if(!first_step){
        temp.Copy((T)1.5,example.face_velocities(c),-(T).5,example.prev_face_velocities(c));
        linear_advection.Update_Advection_Equation_Face_Lookup(example.grid,temp2,lookup_temp,lookup_face_velocities,boundary,dt/2,time+dt);
        linear_advection.Update_Advection_Equation_Face_Lookup(example.grid,temp3,lookup_face_velocities,lookup_face_velocities,boundary,dt,time+dt);
        Extrapolate_Velocity(temp2,example.prev_face_color,c);
        Extrapolate_Velocity(temp3,example.prev_face_color,c);
        quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,temp4,lookup_face_velocities,lookup_temp2,boundary,dt,time+dt);
        quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,temp5,lookup_prev_face_velocities,lookup_temp3,boundary,2*dt,time+dt);
        Extrapolate_Velocity(temp4,example.prev_face_color,c);
        Extrapolate_Velocity(temp5,example.prev_face_color,c);
        example.prev_face_velocities(c).Copy((T)2/(T)1.5,temp4,-(T).5/(T)1.5,temp5);}
    else{
        quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,example.prev_face_velocities(c),lookup_face_velocities,lookup_face_velocities,boundary,dt,time+dt);}
}
//#####################################################################
// Function Apply_Pressure_And_Viscosity
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Apply_Pressure_And_Viscosity(T dt,bool first_step)
{
    if(example.omit_solve) return;
    static int solve_id=-1;solve_id++;
    struct BOUNDARY_CONDITIONS_COLOR_LOCAL:public BOUNDARY_CONDITIONS_COLOR<TV>
    {
        PLS_FC_EXAMPLE<TV>* example;
        T time,dt;
        virtual TV j_surface(const TV& X,int color0,int color1) {return TV();}
        virtual TV d_surface(const TV& X,int color0,int color1) {return example->Dirichlet_Boundary_Condition(X,color0,color1,time);}
        virtual TV n_surface(const TV& X,int color0,int color1) {return example->Neumann_Boundary_Condition(X,color0,color1,time)*dt;}
    } bccl;
    bccl.example=&example;
    bccl.time=time+dt;
    bccl.dt=dt;

    INTERFACE_STOKES_SYSTEM_COLOR<TV> iss(example.grid,example.levelset_color.phi,example.levelset_color.color,true);
    iss.use_preconditioner=example.use_preconditioner;
    ARRAY<T> system_inertia=example.rho,dt_mu(example.mu*dt);
    if(!first_step) system_inertia*=(T)1.5;
    iss.Set_Matrix(dt_mu,example.wrap,&bccl,&system_inertia,&system_inertia);

    printf("\n");
    for(int i=0;i<TV::m;i++){for(int c=0;c<iss.cdi->colors;c++) printf("%c%d [%i]\t","uvw"[i],c,iss.cm_u(i)->dofs(c));printf("\n");}
    for(int c=0;c<iss.cdi->colors;c++) printf("p%d [%i]\t",c,iss.cm_p->dofs(c));printf("\n");
    printf("qn [%i]\t",iss.cdi->constraint_base_n);
    printf("qt [%i] ",iss.cdi->constraint_base_t);
    printf("\n");

    Extrapolate_Velocity(example.face_velocities,example.face_color);
    INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV> rhs,sol;
    {
        ARRAY<ARRAY<TV,TV_INT> > f_volume;
        f_volume.Resize(iss.cdi->colors);
        for(int c=0;c<iss.cdi->colors;c++)
            f_volume(c).Resize(example.grid.Domain_Indices());

        iss.Set_RHS(rhs,f_volume,&example.face_velocities);
        iss.Resize_Vector(sol);
    }

    MINRES<T> mr;
    KRYLOV_SOLVER<T>* solver=&mr;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;

    if(example.dump_matrix){
        KRYLOV_SOLVER<T>::Ensure_Size(vectors,rhs,2);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("M-%d.txt",solve_id).c_str()).Write("M",iss,*vectors(0),*vectors(1));
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("P-%d.txt",solve_id).c_str()).Write_Preconditioner("P",iss,*vectors(0),*vectors(1));
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("b-%d.txt",solve_id).c_str()).Write("b",rhs);}
    solver->Solve(iss,sol,rhs,vectors,1e-10,0,example.max_iter);

    iss.Multiply(sol,*vectors(0));
    *vectors(0)-=rhs;
    LOG::cout<<"Residual: "<<iss.Convergence_Norm(*vectors(0))<<std::endl;

    for(int i=0;i<iss.null_modes.m;i++){
        iss.Multiply(*iss.null_modes(i),*vectors(0));
        LOG::cout<<"null mode["<<i<<"] "<<iss.Convergence_Norm(*vectors(0))<<std::endl;}

    // TODO: Handle ghost region
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(example.grid);it.Valid();it.Next()){
        int c=example.levelset_color.Color(it.Location());
        example.face_color(it.Full_Index())=c;
        if(c<0) continue;
        int k=iss.cm_u(it.Axis())->Get_Index(it.index,c);
        assert(k>=0);
        example.face_velocities(c)(it.Full_Index())=sol.u(it.Axis())(c)(k);}
    vectors.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Extrapolate_Velocity
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Extrapolate_Velocity(ARRAY<T,FACE_INDEX<TV::dimension> >& u,const ARRAY<int,FACE_INDEX<TV::dimension> >& color,int c)
{
    struct MASK_FACE:public EXTRAPOLATION_HIGHER_ORDER_POLY<TV,T>::MASK_FACE
    {
        const ARRAY<int,FACE_INDEX<TV::dimension> >* face_color;
        int c;
        virtual bool Inside(const FACE_INDEX<TV::m>& index)
        {return (*face_color)(index)==c;}
    } mask;
    mask.face_color=&color;
    mask.c=c;

    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(example.grid);it.Valid();it.Next())
        if(color(it.Full_Index())!=c)
            u(it.Full_Index())=0;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(example.grid,3,GRID<TV>::GHOST_REGION);it.Valid();it.Next())
        u(it.Full_Index())=0;
    EXTRAPOLATION_HIGHER_ORDER_POLY<TV,T>::Extrapolate_Face(example.grid,mask,3,u,3,3);
}
//#####################################################################
// Function Extrapolate_Velocity
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Extrapolate_Velocity(ARRAY<ARRAY<T,FACE_INDEX<TV::dimension> > >& u,const ARRAY<int,FACE_INDEX<TV::dimension> >& color)
{
    for(int c=0;c<example.number_of_colors;c++)
        Extrapolate_Velocity(u(c),color,c);
}
//#####################################################################
// Function Swap_State
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Swap_State()
{
    example.prev_face_velocities.Exchange(example.face_velocities);
    example.prev_face_color.Exchange(example.face_color);
}
//#####################################################################
// Simulate_To_Frame
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    for(;current_frame<frame;current_frame++){
        LOG::SCOPE scope("FRAME","Frame %d",current_frame+1);
        for(int substep=0;substep<example.time_steps_per_frame;substep++){
            LOG::SCOPE scope("SUBSTEP","substep %d",substep+1);
            Advance_One_Time_Step(current_frame==0 && substep==0);}
        if(current_frame%1==0 && 0){
            example.particle_levelset_evolution.Reseed_Particles(time);
            example.particle_levelset_evolution.Delete_Particles_Outside_Grid();}
        Write_Output_Files(++output_number);}
} 
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    if(level<=example.write_substeps_level){
        example.frame_title=title;
        LOG::cout<<"Writing substep ["<<example.frame_title<<"]: output_number="<<output_number+1<<", time="<<time<<", frame="<<current_frame<<", substep="<<substep<<std::endl;
        Write_Output_Files(++output_number);
        example.frame_title="";}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Write_Output_Files(const int frame)
{
    FILE_UTILITIES::Create_Directory(example.output_directory);
    FILE_UTILITIES::Create_Directory(example.output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+STRING_UTILITIES::string_sprintf("/%d/frame_title",frame),example.frame_title);
    if(frame==0) 
        FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/first_frame",frame,"\n");
    example.Write_Output_Files(frame);
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
template class PLS_FC_DRIVER<VECTOR<float,2> >;
template class PLS_FC_DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PLS_FC_DRIVER<VECTOR<double,2> >;
template class PLS_FC_DRIVER<VECTOR<double,3> >;
#endif
