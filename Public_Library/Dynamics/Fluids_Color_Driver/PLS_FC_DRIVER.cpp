//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Extrapolation/EXTRAPOLATION_HIGHER_ORDER_POLY.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_HAMILTON_JACOBI_ENO.h>
#include <Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/QUADRATIC_INTERPOLATION_UNIFORM.h>
#include <Tools/Krylov_Solvers/GMRES.h>
#include <Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Log/DEBUG_SUBSTEPS.h>
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Tools/Vectors/VECTOR_UTILITIES.h>
#include <Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <Geometry/Finite_Elements/INTERFACE_STOKES_MULTIGRID.h>
#include <Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_COLOR.h>
#include <Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_VECTOR_COLOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Grids_Uniform_Computations/REINITIALIZATION.h>
#include <Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER.h>
#include <Incompressible/Boundaries/BOUNDARY_PHI_WATER.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM.h>
#include <Dynamics/Fluids_Color_Driver/PLS_FC_DRIVER.h>
#include <Dynamics/Fluids_Color_Driver/PLS_FC_EXAMPLE.h>
#include <Dynamics/Level_Sets/LEVELSET_ADVECTION.h>
#include <Dynamics/Level_Sets/LEVELSET_ADVECTION_MULTIPLE.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#include <boost/function.hpp>
using namespace PhysBAM;
namespace{
    template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
    {
        ((PLS_FC_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
    }
};
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PLS_FC_DRIVER<TV>::
PLS_FC_DRIVER(PLS_FC_EXAMPLE<TV>& example)
    :example(example)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PLS_FC_DRIVER<TV>::
~PLS_FC_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Execute_Main_Program
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
    LOG::cout<<std::setprecision(16)<<std::endl;
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.substeps_delay_frame<0?example.write_substeps_level:-1);

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
    if(example.save_pressure){
        example.pressure_color.Resize(example.grid.Domain_Indices(example.number_of_ghost_cells));
        example.pressure.Resize(example.grid.Domain_Indices(example.number_of_ghost_cells));}
    if(example.use_polymer_stress){
        example.polymer_stress.Resize(example.number_of_colors);
        example.prev_polymer_stress.Resize(example.number_of_colors);
        for(int i=0;i<example.number_of_colors;i++){
            example.polymer_stress(i).Resize(example.grid.Domain_Indices(example.number_of_ghost_cells));
            example.prev_polymer_stress(i).Resize(example.grid.Domain_Indices(example.number_of_ghost_cells));}}

    example.particle_levelset_evolution_multiple.Initialize_Domain(example.grid,example.collision_bodies_affecting_fluid,example.number_of_colors,false);//false= we use positive and negative particles, not just negative
    example.particle_levelset_evolution_multiple.particle_levelset_multiple.Set_Band_Width(2*example.number_of_ghost_cells);
    example.collision_bodies_affecting_fluid.Initialize_Grids();

    example.particle_levelset_evolution_multiple.time=time;
    example.particle_levelset_evolution_multiple.Set_CFL_Number((T).9);
    example.particle_levelset_evolution_multiple.Use_Reinitialization();

    example.particle_levelset_evolution_multiple.levelset_advection_multiple.Set_Custom_Advection(example.advection_scalar);
    for(int i=0;i<example.number_of_colors;i++){
        example.particle_levelset_evolution_multiple.Particle_Levelset(i).number_of_ghost_cells=example.number_of_ghost_cells;}

    example.particle_levelset_evolution_multiple.Set_Number_Particles_Per_Cell(pow<TV::m>(4));
    example.particle_levelset_evolution_multiple.Set_Levelset_Callbacks(example);
    example.particle_levelset_evolution_multiple.Initialize_FMM_Initialization_Iterative_Solver(true);
    example.particle_levelset_evolution_multiple.runge_kutta_order_levelset=3;
    example.particle_levelset_evolution_multiple.runge_kutta_order_particles=3;
    example.particle_levelset_evolution_multiple.Use_Hamilton_Jacobi_Weno_Advection();

    for(int i=0;i<example.number_of_colors;i++){
        example.particle_levelset_evolution_multiple.particle_levelset_multiple.particle_levelsets(i)->levelset.Set_Custom_Boundary(example.boundary);
        example.particle_levelset_evolution_multiple.Particle_Levelset(i).Store_Unique_Particle_Id();
        example.particle_levelset_evolution_multiple.Particle_Levelset(i).Set_Collision_Distance_Factors(.1,1);}
    example.particle_levelset_evolution_multiple.Bias_Towards_Negative_Particles(true);
    example.particle_levelset_evolution_multiple.use_particle_levelset=true;

    example.Initialize();
    if(example.restart){
        example.Read_Output_Files(example.restart);}
    else{
        example.Rebuild_Levelset_Color();
        for(FACE_ITERATOR<TV> it(example.grid,example.number_of_ghost_cells);it.Valid();it.Next())
            example.face_color(it.Full_Index())=example.levelset_color.Color(it.Location());
        example.prev_face_color.Fill(-9);
        example.particle_levelset_evolution_multiple.Make_Signed_Distance();
        example.Make_Levelsets_Consistent();
        example.Get_Initial_Velocities();
        if(example.use_polymer_stress) example.Get_Initial_Polymer_Stresses();}

    example.collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false);
    example.collision_bodies_affecting_fluid.Rasterize_Objects();
    example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*example.grid.dX.Min(),5);
    example.collision_bodies_affecting_fluid.Compute_Grid_Visibility();
    example.particle_levelset_evolution_multiple.Set_Seed(2606);
    if(!example.restart && example.use_pls) example.particle_levelset_evolution_multiple.Seed_Particles(time);
    example.particle_levelset_evolution_multiple.Delete_Particles_Outside_Grid();

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
    for(FACE_ITERATOR<TV> it(example.grid);it.Valid();it.Next()){
        int c=example.face_color(it.Full_Index());
        if(c<0) continue;
        maximum_fluid_speed=max(maximum_fluid_speed,abs(example.face_velocities(c)(it.Full_Index())));}

    T max_particle_collision_distance=example.particle_levelset_evolution_multiple.Particle_Levelset(0).max_collision_distance_factor*example.grid.dX.Max();
    example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(true,dt*maximum_fluid_speed+2*max_particle_collision_distance+(T).5*example.grid.dX.Max(),10);

    LOG::Time("Fill ghost");
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before PLS phi advection",0,1);
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities_ghost(example.grid,example.number_of_ghost_cells);
    example.Merge_Velocities(face_velocities_ghost,example.face_velocities,example.face_color);
    example.boundary.Apply_Boundary_Condition_Face(example.grid,face_velocities_ghost,time+example.dt);
    example.boundary.Fill_Ghost_Faces(example.grid,face_velocities_ghost,face_velocities_ghost,time+dt,example.number_of_ghost_cells);

    LOG::Time("Advect Levelset");
    PHYSBAM_DEBUG_WRITE_SUBSTEP("phi BCs prior to advection",0,1);
    example.particle_levelset_evolution_multiple.Advance_Levelset(dt);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after phi",0,1);
    LOG::Time("advecting particles");
    example.particle_levelset_evolution_multiple.particle_levelset_multiple.Euler_Step_Particles(face_velocities_ghost,dt,time,true,true,false);

    PHYSBAM_DEBUG_WRITE_SUBSTEP("after particle advection",0,1);

    LOG::Time("updating removed particle velocities");
    example.Make_Levelsets_Consistent();

    LOG::Time("modifying levelset");
    example.particle_levelset_evolution_multiple.particle_levelset_multiple.Modify_Levelset_Using_Escaped_Particles(&face_velocities_ghost);
    example.Make_Levelsets_Consistent();
    example.particle_levelset_evolution_multiple.Make_Signed_Distance();
    example.Make_Levelsets_Consistent();
    example.particle_levelset_evolution_multiple.particle_levelset_multiple.Modify_Levelset_Using_Escaped_Particles(&face_velocities_ghost);
    example.Make_Levelsets_Consistent();

    LOG::Time("adjusting particle radii");
    example.particle_levelset_evolution_multiple.particle_levelset_multiple.Adjust_Particle_Radii();
    LOG::Stop_Time();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after particles",0,1);

    LOG::Time("deleting particles"); // needs to be after adding sources, since it does a reseed
    example.particle_levelset_evolution_multiple.Delete_Particles_Outside_Grid();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete 1",0,1);
    LOG::Time("deleting particles in local maxima");
    for(int i=0;i<example.number_of_colors;i++)
        example.particle_levelset_evolution_multiple.Particle_Levelset(i).Delete_Particles_In_Local_Maximum_Phi_Cells(0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete 2",0,1);
    LOG::Time("deleting particles far from interface");
    for(int i=0;i<example.number_of_colors;i++)
        example.particle_levelset_evolution_multiple.Particle_Levelset(i).Delete_Particles_Far_From_Interface(); // uses visibility
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete 3",0,1);

    LOG::Time("re-incorporating removed particles");
    for(int i=0;i<example.number_of_colors;i++)
        example.particle_levelset_evolution_multiple.Particle_Levelset(i).Identify_And_Remove_Escaped_Particles(face_velocities_ghost,5,time+dt);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after remove",0,1);
    example.Rebuild_Levelset_Color();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after rebuilding level set again",0,1);
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

    if(example.use_polymer_stress) Update_Polymer_Stress(dt);

    PHYSBAM_DEBUG_WRITE_SUBSTEP("before velocity advection",0,1);
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
    example.Make_Levelsets_Consistent();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before reinitialization",0,1);
    example.particle_levelset_evolution_multiple.Make_Signed_Distance();
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
    for(FACE_ITERATOR<TV> it(example.grid);it.Valid();it.Next())
        if(color(it.Full_Index())==c){
            PHYSBAM_ASSERT(u(it.Full_Index())*dt<example.grid.dX(it.Axis())*example.number_of_ghost_cells);}
}
//#####################################################################
// Function Advection_And_BDF
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Reduced_Advection_And_BDF(T dt,bool first_step,int c)
{
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T,AVERAGING_UNIFORM<TV>,QUADRATIC_INTERPOLATION_UNIFORM<TV,T> > quadratic_advection;
    BOUNDARY_MAC_GRID_PERIODIC<TV,T> boundary;
    FACE_LOOKUP_UNIFORM<TV> lookup_face_velocities(example.face_velocities(c)),lookup_prev_face_velocities(example.prev_face_velocities(c));
    PHYSBAM_DEBUG_ONLY(Assert_Advection_CFL(example.face_velocities(c),example.face_color,c,dt));
    if(!first_step){
        PHYSBAM_DEBUG_ONLY(Assert_Advection_CFL(example.prev_face_velocities(c),example.prev_face_color,c,dt));
        ARRAY<T,FACE_INDEX<TV::dimension> > temp(example.grid,example.number_of_ghost_cells);
        quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,temp,lookup_prev_face_velocities,lookup_prev_face_velocities,boundary,2*dt,time+dt);
        quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,example.prev_face_velocities(c),lookup_face_velocities,lookup_face_velocities,boundary,dt,time+dt);
        example.prev_face_velocities(c).Copy((T)2/(T)1.5,example.prev_face_velocities(c),-(T).5/(T)1.5,temp);}
    else quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,example.prev_face_velocities(c),lookup_face_velocities,lookup_face_velocities,boundary,dt,time+dt);
}
//#####################################################################
// Function Advection_And_BDF
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
RK2_Advection_And_BDF(T dt,bool first_step,int c)
{
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T,AVERAGING_UNIFORM<TV>,QUADRATIC_INTERPOLATION_UNIFORM<TV,T> > quadratic_advection;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T,AVERAGING_UNIFORM<TV>,LINEAR_INTERPOLATION_UNIFORM<TV,T> > linear_advection;
    BOUNDARY_MAC_GRID_PERIODIC<TV,T> boundary;
    ARRAY<T,FACE_INDEX<TV::dimension> > temp(example.grid,example.number_of_ghost_cells),temp2(example.grid,example.number_of_ghost_cells),
    temp3(example.grid,example.number_of_ghost_cells),temp4(example.grid,example.number_of_ghost_cells),temp5(example.grid,example.number_of_ghost_cells);
    FACE_LOOKUP_UNIFORM<TV> lookup_temp(temp),lookup_temp2(temp2),lookup_temp3(temp3),lookup_temp4(temp4),lookup_temp5(temp5);
    FACE_LOOKUP_UNIFORM<TV> lookup_face_velocities(example.face_velocities(c)),lookup_prev_face_velocities(example.prev_face_velocities(c));
    PHYSBAM_DEBUG_ONLY(Assert_Advection_CFL(example.face_velocities(c),example.face_color,c,dt));
    if(!first_step){
        PHYSBAM_DEBUG_ONLY(Assert_Advection_CFL(example.prev_face_velocities(c),example.prev_face_color,c,dt));
        temp.Copy((T)1.5,example.face_velocities(c),-(T).5,example.prev_face_velocities(c));
        linear_advection.Update_Advection_Equation_Face_Lookup(example.grid,temp2,lookup_temp,lookup_face_velocities,boundary,dt/2,time+dt);
        linear_advection.Update_Advection_Equation_Face_Lookup(example.grid,temp3,lookup_face_velocities,lookup_face_velocities,boundary,dt,time+dt);
        quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,temp4,lookup_face_velocities,lookup_temp2,boundary,dt,time+dt);
        quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,temp5,lookup_prev_face_velocities,lookup_temp3,boundary,2*dt,time+dt);
        example.prev_face_velocities(c).Copy((T)2/(T)1.5,temp4,-(T).5/(T)1.5,temp5);}
    else quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,example.prev_face_velocities(c),lookup_face_velocities,lookup_face_velocities,boundary,dt,time+dt);
}
//#####################################################################
// Function Apply_Pressure_And_Viscosity
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Apply_Pressure_And_Viscosity(T dt,bool first_step)
{
    if(example.omit_solve) return;
    static int solve_id=-1;solve_id++;
    INTERFACE_STOKES_SYSTEM_COLOR<TV>* issp=0;

    ARRAY<ARRAY<T,TV_INT> >& phis=example.particle_levelset_evolution_multiple.particle_levelset_multiple.levelset_multiple.phis;
    if(example.use_multigrid)
        issp=new INTERFACE_STOKES_MULTIGRID<TV>(example.num_multigrid_levels,example.grid,example.levelset_color.phi,
            example.levelset_color.color,true,phis,example.bc_phis,example.number_of_ghost_cells);
    else issp=new INTERFACE_STOKES_SYSTEM_COLOR<TV>(example.grid,example.levelset_color.phi,example.levelset_color.color,true);
    INTERFACE_STOKES_SYSTEM_COLOR<TV>& iss=*issp;

    iss.use_preconditioner=example.use_preconditioner;
    iss.use_p_null_mode=example.use_p_null_mode;
    iss.use_polymer_stress=example.use_polymer_stress;
    ARRAY<T> inertia=example.rho,dt_mu(example.mu*dt);
    if(!first_step) inertia*=(T)1.5;
    iss.Set_Matrix(dt_mu,example.use_discontinuous_velocity,
        [=](const TV& X,int color0,int color1){return example.Velocity_Jump(X,color0,color1,time);},
        [=](const TV& X,int color0,int color1){return example.Jump_Interface_Condition(X,color0,color1,time)*dt;},
        &inertia,true);

    printf("\n");
    for(int i=0;i<TV::m;i++){for(int c=0;c<iss.cdi->colors;c++) printf("%c%d [%i]\t","uvw"[i],c,iss.cm_u(i)->dofs(c));printf("\n");}
    for(int c=0;c<iss.cdi->colors;c++) printf("p%d [%i]\t",c,iss.cm_p->dofs(c));printf("\n");
    printf("qn [%i]\t",iss.cdi->constraint_base_n);
    printf("qt [%i] ",iss.cdi->constraint_base_t);
    printf("\n");

    INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV> rhs,sol;
    iss.Set_RHS(rhs,[=](const TV& X,int color){return example.Volume_Force(X,color,time)*dt;},&example.face_velocities,false);
    if(example.use_polymer_stress) iss.Add_Polymer_Stress_RHS(rhs,example.polymer_stress); //if needed
    iss.Resize_Vector(sol);

    MINRES<T> mr;
    KRYLOV_SOLVER<T>* solver=&mr;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;

    if(example.dump_matrix){
        KRYLOV_SOLVER<T>::Ensure_Size(vectors,rhs,2);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("M-%d.txt",solve_id).c_str()).Write("M",iss,*vectors(0),*vectors(1));
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("Z-%d.txt",solve_id).c_str()).Write_Preconditioner("Z",iss,*vectors(0),*vectors(1));
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("P-%d.txt",solve_id).c_str()).Write_Projection("P",iss,*vectors(0));
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("b-%d.txt",solve_id).c_str()).Write("b",rhs);}
    if(example.sparse_dump_matrix){
        SPARSE_MATRIX_FLAT_MXN<T> M;
        iss.Get_Sparse_Matrix(M);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("M-%d.txt",solve_id).c_str()).Write("M",M);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("b-%d.txt",solve_id).c_str()).Write("b",rhs);}
    solver->Solve(iss,sol,rhs,vectors,1e-10,0,example.max_iter);

    if(example.dump_matrix){
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("x-%d.txt",solve_id).c_str()).Write("x",sol);}
    if(example.dump_largest_eigenvector) Dump_Largest_Eigenvector(iss,vectors);

    iss.Multiply(sol,*vectors(0));
    *vectors(0)-=rhs;
    LOG::cout<<"Residual: "<<iss.Convergence_Norm(*vectors(0))<<std::endl;

    for(int i=0;i<iss.null_modes.m;i++){
        iss.Multiply(*iss.null_modes(i),*vectors(0));
        LOG::cout<<"null mode["<<i<<"] "<<iss.Convergence_Norm(*vectors(0))<<std::endl;}

    for(FACE_ITERATOR<TV> it(example.grid);it.Valid();it.Next()){
        int c=example.levelset_color.Color(it.Location());
        example.face_color(it.Full_Index())=c;
        if(c<0) continue;
        int k=iss.cm_u(it.Axis())->Get_Index(it.index,c);
        assert(k>=0);
        example.face_velocities(c)(it.Full_Index())=sol.u(it.Axis())(c)(k);}

    for(FACE_ITERATOR<TV> it(example.grid,example.number_of_ghost_cells,GRID<TV>::GHOST_REGION);it.Valid();it.Next())
        example.face_color(it.Full_Index())=example.levelset_color.Color(it.Location());
    vectors.Delete_Pointers_And_Clean_Memory();

    for(int c=0;c<example.face_velocities.m;c++){
        example.boundary.Apply_Boundary_Condition_Face(example.grid,example.face_velocities(c),time+example.dt);
        example.boundary.Fill_Ghost_Faces(example.grid,example.face_velocities(c),example.face_velocities(c),0,example.number_of_ghost_cells);}
    example.boundary_int.Apply_Boundary_Condition_Face(example.grid,example.face_color,time+example.dt);
    example.boundary_int.Fill_Ghost_Faces(example.grid,example.face_color,example.face_color,0,example.number_of_ghost_cells);

    if(example.save_pressure){
        for(CELL_ITERATOR<TV> it(example.grid);it.Valid();it.Next()){
            int c=example.levelset_color.Color(it.Location());
            if(c<0) continue;
            int k=iss.cm_p->Get_Index(it.index,c);
            assert(k>=0);
            example.pressure(it.index)=sol.p(c)(k)/dt;
            example.pressure_color(it.index)=c;}}
}
//#####################################################################
// Function Update_Polymer_Stress
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Update_Polymer_Stress(T dt)
{

    example.prev_polymer_stress.Exchange(example.polymer_stress);
        //Right now we are filling in for all colors. We may not want to do that in the future.

            for(CELL_ITERATOR<TV> it(example.grid,1);it.Valid();it.Next()){
                for(int c=0;c<example.polymer_stress.m;c++)
                    example.polymer_stress(c)(it.index)=example.Polymer_Stress(it.Location(),c,time+dt);}

}
//#####################################################################
// Function Extrapolate_Velocity
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Extrapolate_Velocity(ARRAY<T,FACE_INDEX<TV::dimension> >& u,const ARRAY<int,FACE_INDEX<TV::dimension> >& color,int c)
{
    for(FACE_ITERATOR<TV> it(example.grid);it.Valid();it.Next())
        if(color(it.Full_Index())!=c)
            u(it.Full_Index())=1e20;
    for(FACE_ITERATOR<TV> it(example.grid,example.number_of_ghost_cells,GRID<TV>::GHOST_REGION);it.Valid();it.Next())
        u(it.Full_Index())=1e20;

    const LEVELSET<TV>& phi=*example.particle_levelset_evolution_multiple.particle_levelset_multiple.levelset_multiple.levelsets(c);
    EXTRAPOLATION_HIGHER_ORDER<TV,T> eho(example.grid,phi,example.number_of_ghost_cells*10,3,example.number_of_ghost_cells);
    eho.periodic=true;

    eho.Extrapolate_Face([&](const FACE_INDEX<TV::m>& index){return color(index)==c;},u);

    example.boundary.Apply_Boundary_Condition_Face(example.grid,u,time+example.dt);
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
        LOG::SCOPE scope("FRAME","frame %d",current_frame+1);
        if(example.substeps_delay_frame==current_frame)
            DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);
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
// Function Dump_Largest_Eigenvector
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Dump_Largest_Eigenvector(const INTERFACE_STOKES_SYSTEM_COLOR<TV>& iss,ARRAY<KRYLOV_VECTOR_BASE<T>*>& vectors) const
{
    static int solve_id=-1;solve_id++;
    INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV> sol,rhs,tmp;
    sol.Resize(*vectors(0));
    rhs.Resize(*vectors(0));
    tmp.Resize(*vectors(0));

    MINRES<T> mr;
    KRYLOV_SOLVER<T>* solver=&mr;

    ARRAY<INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>*> evs;
    RANDOM_NUMBERS<T> random;

    OCTAVE_OUTPUT<T> oo(STRING_UTILITIES::string_sprintf("ev-%d.txt",solve_id).c_str());

    for(int j=0;j<5;j++){
        random.Fill_Uniform(sol.u,-1,1);
        for(int i=0;i<evs.m;i++) sol.Copy(-sol.Dot(*evs(i)),*evs(i),sol);
        sol*=1/sqrt((T)iss.Inner_Product(sol,sol));
        T a;
        for(int i=0;i<20;i++){
            for(int d=0;d<TV::m;d++)
                for(int c=0;c<example.number_of_colors;c++)
                    iss.matrix_inertial_rhs(d)(c).Times(sol.u(d)(c),rhs.u(d)(c));
            iss.Project(rhs);

            tmp=sol;
            sol*=0;
            solver->Solve(iss,sol,rhs,vectors,1e-10,0,example.max_iter);
            iss.Project(sol);
            for(int c=0;c<example.number_of_colors;c++) sol.p(c).Fill(0);
            sol.q.Fill(0);
            for(int i=0;i<evs.m;i++) sol.Copy(-sol.Dot(*evs(i)),*evs(i),sol);

            a=sqrt((T)iss.Inner_Product(sol,sol));
            sol*=1/a;
            LOG::cout<<"eig approx "<<a<<"   dp "<<iss.Inner_Product(sol,tmp)<<std::endl;
            tmp=sol;}

        Dump_Vector(iss,sol,STRING_UTILITIES::string_sprintf("eigenvector (%g)",a).c_str());
        evs.Append(new INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>(sol));
        oo.Write(STRING_UTILITIES::string_sprintf("EV%d",j).c_str(),sol);
        oo.Write(STRING_UTILITIES::string_sprintf("ev%d",j).c_str(),a);}

    evs.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Dump_Vector
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Dump_Vector(const INTERFACE_STOKES_SYSTEM_COLOR<TV>& iss,const KRYLOV_VECTOR_BASE<T>& u,const char* str) const
{
    const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>& v=static_cast<const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>&>(u);
    ARRAY<ARRAY<T,FACE_INDEX<TV::dimension> > > fv=example.face_velocities;
    for(int i=0;i<example.face_velocities.m;i++)
        for(int j=0;j<TV::m;j++)
            example.face_velocities(i).Component(j).Fill(0);

    for(FACE_ITERATOR<TV> it(example.grid);it.Valid();it.Next()){
        int c=example.face_color(it.Full_Index());
        if(c<0) continue;
        int k=iss.cm_u(it.Axis())->Get_Index(it.index,c);
        assert(k>=0);
        example.face_velocities(c)(it.Full_Index())=v.u(it.Axis())(c)(k);}

    for(CELL_ITERATOR<TV> it(example.grid);it.Valid();it.Next()){
        int c=example.levelset_color.Color(it.Location());
        if(c<0) continue;
        int k=iss.cm_p->Get_Index(it.index,c);
        assert(k>=0);
        T p=v.p(c)(k);
        Add_Debug_Particle(it.Location(),VECTOR<T,3>(p<0,p>=0,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(p));}

    PHYSBAM_DEBUG_WRITE_SUBSTEP(str,0,1);
    example.face_velocities=fv;
}
//#####################################################################
namespace PhysBAM{
    template class PLS_FC_DRIVER<VECTOR<float,2> >;
    template class PLS_FC_DRIVER<VECTOR<float,3> >;
    template class PLS_FC_DRIVER<VECTOR<double,2> >;
    template class PLS_FC_DRIVER<VECTOR<double,3> >;
}
