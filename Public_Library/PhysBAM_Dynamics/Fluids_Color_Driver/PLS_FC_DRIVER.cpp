//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/EXTRAPOLATION_HIGHER_ORDER_POLY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_HAMILTON_JACOBI_ENO.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/QUADRATIC_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <PhysBAM_Tools/Krylov_Solvers/MINRES.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Vectors/VECTOR_UTILITIES.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_VECTOR_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/VOLUME_FORCE_COLOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/REINITIALIZATION.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/BOUNDARY_PHI_WATER.h>
#include <PhysBAM_Dynamics/Fluids_Color_Driver/PLS_FC_DRIVER.h>
#include <PhysBAM_Dynamics/Fluids_Color_Driver/PLS_FC_EXAMPLE.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_MULTIPLE.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
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
    example.time=time=example.time_steps_per_frame*current_frame*example.dt;

    example.levelset_color.phi.Resize(example.grid.Domain_Indices(example.number_of_ghost_cells));
    example.levelset_color.color.Resize(example.grid.Domain_Indices(example.number_of_ghost_cells));
    example.face_color.Resize(example.grid,example.number_of_ghost_cells);
    example.prev_face_color.Resize(example.grid,example.number_of_ghost_cells);
    example.face_velocities.Resize(example.number_of_colors);
    example.prev_face_velocities.Resize(example.number_of_colors);
    for(int i=0;i<example.number_of_colors;i++){
        example.face_velocities(i).Resize(example.grid,example.number_of_ghost_cells);
        example.prev_face_velocities(i).Resize(example.grid,example.number_of_ghost_cells);}
    for(int i=0;i<example.bc_phis.m;i++)
        example.bc_phis(i).Resize(example.grid.Domain_Indices(example.number_of_ghost_cells));

    {
        example.particle_levelset_evolution_multiple.Initialize_Domain(example.grid,example.number_of_colors);
        example.particle_levelset_evolution_multiple.particle_levelset.Set_Band_Width(2*example.number_of_ghost_cells);
        example.collision_bodies_affecting_fluid.Initialize_Grids();
    }

    example.particle_levelset_evolution_multiple.time=time;
    example.particle_levelset_evolution_multiple.Set_CFL_Number((T).9);
    example.particle_levelset_evolution_multiple.Use_Reinitialization();

    VECTOR<VECTOR<bool,2>,TV::dimension> domain_open_boundaries=VECTOR_UTILITIES::Complement(example.domain_boundary);
    example.particle_levelset_evolution_multiple.Levelset_Advection(0).Set_Custom_Advection(example.advection_scalar);

    {
        example.particle_levelset_evolution_multiple.Initialize_Domain(example.grid,example.number_of_colors);
        example.particle_levelset_evolution_multiple.particle_levelset.Set_Band_Width(2*example.number_of_ghost_cells);
        example.collision_bodies_affecting_fluid.Initialize_Grids();
    }

    example.particle_levelset_evolution_multiple.Particle_Levelset(0).number_of_ghost_cells=example.number_of_ghost_cells;
    example.particle_levelset_evolution_multiple.Set_Number_Particles_Per_Cell(16);
    example.particle_levelset_evolution_multiple.Set_Levelset_Callbacks(example);
    example.particle_levelset_evolution_multiple.Initialize_FMM_Initialization_Iterative_Solver(true);
    example.particle_levelset_evolution_multiple.runge_kutta_order_levelset=3;
    example.particle_levelset_evolution_multiple.runge_kutta_order_particles=3;
    example.particle_levelset_evolution_multiple.Use_Hamilton_Jacobi_Weno_Advection();

    example.particle_levelset_evolution_multiple.particle_levelset.levelset.Set_Custom_Boundary(example.boundary);
    example.particle_levelset_evolution_multiple.Bias_Towards_Negative_Particles(false);
    example.particle_levelset_evolution_multiple.Particle_Levelset(0).Store_Unique_Particle_Id();
    example.particle_levelset_evolution_multiple.use_particle_levelset=true;
    example.particle_levelset_evolution_multiple.particle_levelset.levelset.Set_Collision_Body_List(example.collision_bodies_affecting_fluid);
// TODO    example.particle_levelset_evolution_multiple.particle_levelset.levelset.Set_Face_Velocities_Valid_Mask(&example.incompressible.valid_mask);
    example.particle_levelset_evolution_multiple.particle_levelset.Set_Collision_Distance_Factors(.1,1);

    {
        example.particle_levelset_evolution_multiple.Initialize_Domain(example.grid,example.number_of_colors);
        example.particle_levelset_evolution_multiple.particle_levelset.Set_Band_Width(2*example.number_of_ghost_cells);
        example.collision_bodies_affecting_fluid.Initialize_Grids();
    }
    example.particle_levelset_evolution_multiple.levelset_advection_multiple.Set_Custom_Advection(*new ADVECTION_HAMILTON_JACOBI_ENO<GRID<TV>,T>);

    if(example.restart){
        example.Read_Output_Files(example.restart);
        example.collision_bodies_affecting_fluid.Rasterize_Objects();
        example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*example.grid.dX.Min(),5);} // compute grid visibility (for advection later)
    else{
        example.collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false);
        example.collision_bodies_affecting_fluid.Rasterize_Objects();
        example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*example.grid.dX.Min(),5);
        example.Initialize();
        example.Rebuild_Levelset_Color();
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(example.grid,example.number_of_ghost_cells);it.Valid();it.Next())
            example.face_color(it.Full_Index())=example.levelset_color.Color(it.Location());
        example.prev_face_color.Fill(-9);
        example.particle_levelset_evolution_multiple.Make_Signed_Distance();
        example.particle_levelset_evolution_multiple.Fill_Levelset_Ghost_Cells(time);
        example.Get_Initial_Velocities();}

    example.collision_bodies_affecting_fluid.Compute_Grid_Visibility();
    example.particle_levelset_evolution_multiple.Set_Seed(2606);
//    Put these back in for PLS
//    if(!example.restart) example.particle_levelset_evolution_multiple.Seed_Particles(time);
//    example.particle_levelset_evolution_multiple.Delete_Particles_Outside_Grid();

    int extrapolation_cells=2*example.number_of_ghost_cells+2;
    ARRAY<T,TV_INT> exchanged_phi_ghost(example.grid.Domain_Indices(extrapolation_cells));
    example.particle_levelset_evolution_multiple.particle_levelset.levelset.boundary->Fill_Ghost_Cells(example.grid,example.particle_levelset_evolution_multiple.phi,exchanged_phi_ghost,0,time,extrapolation_cells);
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

    T max_particle_collision_distance=example.particle_levelset_evolution_multiple.Particle_Levelset(0).max_collision_distance_factor*example.grid.dX.Max();
    example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(true,dt*maximum_fluid_speed+2*max_particle_collision_distance+(T).5*example.grid.dX.Max(),10);

    LOG::Time("Fill ghost");
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before advection",0,1);
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities_ghost;
    face_velocities_ghost.Resize(example.grid,example.number_of_ghost_cells,false);
    // TODO example.incompressible.boundary->Fill_Ghost_Faces(example.grid,example.face_velocities,face_velocities_ghost,time+dt,example.number_of_ghost_cells);

    ARRAY<T,TV_INT> phi_back(example.grid.Domain_Indices(example.number_of_ghost_cells));
    example.boundary.Fill_Ghost_Cells(example.grid,example.particle_levelset_evolution_multiple.Particle_Levelset(0).levelset.phi,phi_back,dt,time,example.number_of_ghost_cells);
    LOG::Time("Advect Levelset");
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before phi",0,1);
    example.particle_levelset_evolution_multiple.Advance_Levelset(dt);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after phi",0,1);
    LOG::Time("advecting particles");
    example.particle_levelset_evolution_multiple.Particle_Levelset(0).Euler_Step_Particles(face_velocities_ghost,dt,time,true,true,false,false);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after advection",0,1);

    LOG::Time("updating removed particle velocities");
        example.particle_levelset_evolution_multiple.Fill_Levelset_Ghost_Cells(time);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after particles",0,1);

    example.particle_levelset_evolution_multiple.Make_Signed_Distance();
    // TODO example.boundary->Fill_Ghost_Faces(example.grid,example.face_velocities,face_velocities_ghost,time+dt,example.number_of_ghost_cells);

    LOG::Time("modifying levelset");
    example.particle_levelset_evolution_multiple.Fill_Levelset_Ghost_Cells(time+dt);
    example.particle_levelset_evolution_multiple.Particle_Levelset(0).Exchange_Overlap_Particles();
    example.particle_levelset_evolution_multiple.Modify_Levelset_And_Particles(&face_velocities_ghost);
    example.particle_levelset_evolution_multiple.Make_Signed_Distance();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after particles",0,1);

    LOG::Time("deleting particles"); // needs to be after adding sources, since it does a reseed
    example.particle_levelset_evolution_multiple.Delete_Particles_Outside_Grid();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete 1",0,1);
    LOG::Time("deleting particles in local maxima");
    example.particle_levelset_evolution_multiple.particle_levelset.Delete_Particles_In_Local_Maximum_Phi_Cells(0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete 2",0,1);
    LOG::Time("deleting particles far from interface");
    example.particle_levelset_evolution_multiple.Particle_Levelset(0).Delete_Particles_Far_From_Interface(); // uses visibility
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete 3",0,1);

    LOG::Time("re-incorporating removed particles");
    example.particle_levelset_evolution_multiple.Particle_Levelset(0).Identify_And_Remove_Escaped_Particles(face_velocities_ghost,1.5,time+dt);
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
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before levelset evolution",0,1);

    if(example.use_level_set_method) Update_Level_Set(dt,first_step);
    else if(example.use_pls) Update_Pls(dt);

    PHYSBAM_DEBUG_WRITE_SUBSTEP("before advection",0,1);
    Extrapolate_Velocity(example.face_velocities,example.face_color);
    Advection_And_BDF(dt,first_step);
    Extrapolate_Velocity(example.face_velocities,example.face_color);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before solve",0,1);
    Apply_Pressure_And_Viscosity(dt,first_step);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after solve",0,1);
    example.time=time+=dt;
    Extrapolate_Velocity(example.face_velocities,example.face_color);
    example.End_Time_Step(time);
}
//#####################################################################
// Function Update_Level_Set
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Update_Level_Set(T dt,bool first_step)
{
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before phi advection",0,1);
    example.particle_levelset_evolution_multiple.Advance_Levelset(dt);
    Enforce_Phi_Boundary_Conditions();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before reinitialization",0,1);
    example.particle_levelset_evolution_multiple.Make_Signed_Distance();
    Enforce_Phi_Boundary_Conditions();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after level set update",0,1);
    example.Rebuild_Levelset_Color();
}
//#####################################################################
// Function Advection_And_BDF
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Advection_And_BDF(T dt,bool first_step)
{
    if(!example.use_advection)
        for(int c=0;c<example.number_of_colors;c++)
            No_Advection_And_BDF(dt,first_step,c);
    else if(example.use_reduced_advection)
        for(int c=0;c<example.number_of_colors;c++)
            Reduced_Advection_And_BDF(dt,first_step,c);
    else
        for(int c=0;c<example.number_of_colors;c++)
            RK2_Advection_And_BDF(dt,first_step,c);
    example.prev_face_velocities.Exchange(example.face_velocities);
    example.prev_face_color=example.face_color;
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
// Function Assert_Advection_CFL
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Assert_Advection_CFL(const ARRAY<T,FACE_INDEX<TV::m> >& u,const ARRAY<int,FACE_INDEX<TV::m> >& color,int c,T dt) const
{
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(example.grid);it.Valid();it.Next())
        if(color(it.Full_Index())==c){
            PHYSBAM_ASSERT(u(it.Full_Index())*dt<example.grid.dX(it.Axis())*example.number_of_ghost_cells);}
}
//#####################################################################
// Function Advection_And_BDF
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Reduced_Advection_And_BDF(T dt,bool first_step,int c)
{
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T,AVERAGING_UNIFORM<GRID<TV> >,QUADRATIC_INTERPOLATION_UNIFORM<GRID<TV>,T> > quadratic_advection;
    BOUNDARY<TV,T> boundary;
    FACE_LOOKUP_UNIFORM<GRID<TV> > lookup_face_velocities(example.face_velocities(c)),lookup_prev_face_velocities(example.prev_face_velocities(c));
    PHYSBAM_DEBUG_ONLY(Assert_Advection_CFL(example.face_velocities(c),example.face_color,c,dt));
    if(!first_step){
        PHYSBAM_DEBUG_ONLY(Assert_Advection_CFL(example.prev_face_velocities(c),example.prev_face_color,c,dt));
        ARRAY<T,FACE_INDEX<TV::dimension> > temp(example.grid,example.number_of_ghost_cells);
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
    BOUNDARY<TV,T> boundary;
    ARRAY<T,FACE_INDEX<TV::dimension> > temp(example.grid,example.number_of_ghost_cells),temp2(example.grid,example.number_of_ghost_cells),
        temp3(example.grid,example.number_of_ghost_cells),temp4(example.grid,example.number_of_ghost_cells),temp5(example.grid,example.number_of_ghost_cells);
    FACE_LOOKUP_UNIFORM<GRID<TV> > lookup_temp(temp),lookup_temp2(temp2),lookup_temp3(temp3),lookup_temp4(temp4),lookup_temp5(temp5);
    FACE_LOOKUP_UNIFORM<GRID<TV> > lookup_face_velocities(example.face_velocities(c)),lookup_prev_face_velocities(example.prev_face_velocities(c));
    PHYSBAM_DEBUG_ONLY(Assert_Advection_CFL(example.face_velocities(c),example.face_color,c,dt));
    if(!first_step){
        PHYSBAM_DEBUG_ONLY(Assert_Advection_CFL(example.prev_face_velocities(c),example.prev_face_color,c,dt));
        temp.Copy((T)1.5,example.face_velocities(c),-(T).5,example.prev_face_velocities(c));
        linear_advection.Update_Advection_Equation_Face_Lookup(example.grid,temp2,lookup_temp,lookup_face_velocities,boundary,dt/2,time+dt);
        linear_advection.Update_Advection_Equation_Face_Lookup(example.grid,temp3,lookup_face_velocities,lookup_face_velocities,boundary,dt,time+dt);
        quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,temp4,lookup_face_velocities,lookup_temp2,boundary,dt,time+dt);
        quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,temp5,lookup_prev_face_velocities,lookup_temp3,boundary,2*dt,time+dt);
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
        virtual TV u_jump(const TV& X,int color0,int color1) {return example->Velocity_Jump(X,color0,color1,time);}
        virtual TV j_surface(const TV& X,int color0,int color1) {return example->Jump_Interface_Condition(X,color0,color1,time)*dt;}
    } bccl;
    bccl.example=&example;
    bccl.time=time+dt;
    bccl.dt=dt;
    bccl.use_discontinuous_velocity=example.use_discontinuous_velocity;

    INTERFACE_STOKES_SYSTEM_COLOR<TV> iss(example.grid,example.levelset_color.phi,example.levelset_color.color,true);
    iss.use_preconditioner=example.use_preconditioner;
    iss.use_p_null_mode=example.use_p_null_mode;
    ARRAY<T> system_inertia=example.rho,dt_mu(example.mu*dt);
    if(!first_step) system_inertia*=(T)1.5;
    iss.Set_Matrix(dt_mu,example.wrap,&bccl,&system_inertia,&system_inertia);

    printf("\n");
    for(int i=0;i<TV::m;i++){for(int c=0;c<iss.cdi->colors;c++) printf("%c%d [%i]\t","uvw"[i],c,iss.cm_u(i)->dofs(c));printf("\n");}
    for(int c=0;c<iss.cdi->colors;c++) printf("p%d [%i]\t",c,iss.cm_p->dofs(c));printf("\n");
    printf("qn [%i]\t",iss.cdi->constraint_base_n);
    printf("qt [%i] ",iss.cdi->constraint_base_t);
    printf("\n");

    INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV> rhs,sol;

    struct VOLUME_FORCE_COLOR_LOCAL:public VOLUME_FORCE_COLOR<TV>
    {
        PLS_FC_EXAMPLE<TV>* example;
        T time,dt;
        virtual TV F(const TV& X,int color) {return example->Volume_Force(X,color,time)*dt;}
    } vfcl;
    vfcl.example=&example;
    vfcl.time=time+dt;
    vfcl.dt=dt;

    iss.Set_RHS(rhs,&vfcl,&example.face_velocities,false);
    iss.Resize_Vector(sol);

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

    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(example.grid);it.Valid();it.Next()){
        int c=example.levelset_color.Color(it.Location());
        example.face_color(it.Full_Index())=c;
        if(c<0) continue;
        int k=iss.cm_u(it.Axis())->Get_Index(it.index,c);
        assert(k>=0);
        example.face_velocities(c)(it.Full_Index())=sol.u(it.Axis())(c)(k);}
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(example.grid,example.number_of_ghost_cells,GRID<TV>::GHOST_REGION);it.Valid();it.Next())
        example.face_color(it.Full_Index())=example.levelset_color.Color(it.Location());
    vectors.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Extrapolate_Velocity
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Extrapolate_Velocity(ARRAY<T,FACE_INDEX<TV::dimension> >& u,const ARRAY<int,FACE_INDEX<TV::dimension> >& color,int c)
{
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(example.grid);it.Valid();it.Next())
        if(color(it.Full_Index())!=c)
            u(it.Full_Index())=0;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(example.grid,example.number_of_ghost_cells,GRID<TV>::GHOST_REGION);it.Valid();it.Next())
        u(it.Full_Index())=0;

    const LEVELSET<TV>& phi=*example.particle_levelset_evolution_multiple.particle_levelset_multiple.levelset_multiple.levelsets(c);
    EXTRAPOLATION_HIGHER_ORDER<TV,T> eho(example.grid,phi,example.number_of_ghost_cells*4,3,example.number_of_ghost_cells);
    eho.periodic=true;
    eho.Extrapolate_Face([&](const FACE_INDEX<TV::m>& index){return color(index)==c;},u);
    example.boundary.Fill_Ghost_Faces(example.grid,u,u,0,example.number_of_ghost_cells);
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
            example.particle_levelset_evolution_multiple.Reseed_Particles(time);
            example.particle_levelset_evolution_multiple.Delete_Particles_Outside_Grid();}
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
// Function Enforce_Phi_Boundary_Conditions
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Enforce_Phi_Boundary_Conditions()
{
    ARRAY<ARRAY<T,TV_INT> >& phis=example.particle_levelset_evolution_multiple.particle_levelset_multiple.levelset_multiple.phis;
    for(int i=0;i<phis.m;i++){
        example.boundary.Fill_Ghost_Cells(example.grid,phis(i),phis(i),0,0,example.number_of_ghost_cells);
        if(phis(i).array.Min()>=0) phis(i).array.Fill(example.number_of_ghost_cells*example.grid.dX.Max());}
}
//#####################################################################
namespace PhysBAM{
template class PLS_FC_DRIVER<VECTOR<float,2> >;
template class PLS_FC_DRIVER<VECTOR<float,3> >;
template class PLS_FC_DRIVER<VECTOR<double,2> >;
template class PLS_FC_DRIVER<VECTOR<double,3> >;
}
