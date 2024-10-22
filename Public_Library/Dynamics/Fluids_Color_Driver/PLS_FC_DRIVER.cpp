//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Core/Math_Tools/pow.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/VECTOR_UTILITIES.h>
#include <Tools/Krylov_Solvers/GMRES.h>
#include <Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Grid_Tools/Extrapolation/EXTRAPOLATION_HIGHER_ORDER_POLY.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Grid_PDE/Advection/ADVECTION_HAMILTON_JACOBI_ENO.h>
#include <Grid_PDE/Boundaries/BOUNDARY_MPI.h>
#include <Grid_PDE/Interpolation/AVERAGING_UNIFORM.h>
#include <Grid_PDE/Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <Grid_PDE/Interpolation/QUADRATIC_INTERPOLATION_UNIFORM.h>
#include <Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <Geometry/Finite_Elements/INTERFACE_STOKES_MULTIGRID.h>
#include <Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_COLOR.h>
#include <Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_VECTOR_COLOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER.h>
#include <Geometry/Level_Sets/REINITIALIZATION.h>
#include <Incompressible/Boundaries/BOUNDARY_PHI_WATER.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM.h>
#include <Dynamics/Fluids_Color_Driver/PLS_FC_DRIVER.h>
#include <Dynamics/Fluids_Color_Driver/PLS_FC_EXAMPLE.h>
#include <Dynamics/Level_Sets/LEVELSET_ADVECTION.h>
#include <Dynamics/Level_Sets/LEVELSET_ADVECTION_MULTIPLE.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#include <functional>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PLS_FC_DRIVER<TV>::
PLS_FC_DRIVER(PLS_FC_EXAMPLE<TV>& example)
    :example(example)
{
    DEBUG_SUBSTEPS::write_substeps_level=example.write_substeps_level;
    DEBUG_SUBSTEPS::writer=[=](const std::string& title){Write_Substep(title);};
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PLS_FC_DRIVER<TV>::
~PLS_FC_DRIVER()
{
    DEBUG_SUBSTEPS::writer=0;
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
    DEBUG_SUBSTEPS::write_substeps_level=example.substeps_delay_frame<0?example.write_substeps_level:-1;

    if(example.auto_restart){
        example.viewer_dir.Read_Last_Frame(0);
        example.restart=example.viewer_dir.frame_stack(0);}
    else if(example.restart)
        example.viewer_dir.Set(example.restart);
    if(example.restart) example.viewer_dir.Make_Common_Directory(true);
    current_frame=example.restart;
    example.time=time=current_frame*example.dt;

    example.levelset_color.phi.Resize(example.grid.Domain_Indices(example.number_of_ghost_cells));
    example.levelset_color.color.Resize(example.grid.Domain_Indices(example.number_of_ghost_cells));
    example.face_color.Resize(example.grid,example.number_of_ghost_cells);
    example.prev_face_color.Resize(example.grid,example.number_of_ghost_cells);
    example.cell_color.Resize(example.grid.Cell_Indices(example.number_of_ghost_cells));
    example.prev_cell_color.Resize(example.grid.Cell_Indices(example.number_of_ghost_cells));
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
        example.Read_Output_Files();}
    else{
        example.Rebuild_Levelset_Color();
        for(FACE_ITERATOR<TV> it(example.grid,example.number_of_ghost_cells);it.Valid();it.Next())
            example.face_color(it.Full_Index())=example.levelset_color.Color(it.Location());
        for(CELL_ITERATOR<TV> it(example.grid,example.number_of_ghost_cells);it.Valid();it.Next())
            example.cell_color(it.index)=example.levelset_color.Color(it.Location());
        example.prev_face_color.Fill(-9);
        example.prev_cell_color.Fill(-9);
        example.particle_levelset_evolution_multiple.Make_Signed_Distance();
        example.Make_Levelsets_Consistent();
        example.Get_Initial_Velocities(0);
        if(example.use_polymer_stress) example.Get_Initial_Polymer_Stresses(0);}

    example.collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false);
    example.collision_bodies_affecting_fluid.Rasterize_Objects();
    example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*example.grid.dX.Min(),5);
    example.collision_bodies_affecting_fluid.Compute_Grid_Visibility();
    example.particle_levelset_evolution_multiple.Set_Seed(2606);
    if(!example.restart && example.use_pls) example.particle_levelset_evolution_multiple.Seed_Particles(time);
    example.particle_levelset_evolution_multiple.Delete_Particles_Outside_Grid();

    if(!example.restart) Write_Output_Files();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after init",1);
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
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before PLS phi advection",1);
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities_ghost(example.grid,example.number_of_ghost_cells);
    example.Merge_Velocities(face_velocities_ghost,example.face_velocities,example.face_color);
    example.boundary.Apply_Boundary_Condition_Face(example.grid,face_velocities_ghost,time+example.dt);
    example.boundary.Fill_Ghost_Faces(example.grid,face_velocities_ghost,face_velocities_ghost,time+dt,example.number_of_ghost_cells);

    LOG::Time("Advect Levelset");
    PHYSBAM_DEBUG_WRITE_SUBSTEP("phi BCs prior to advection",1);
    example.particle_levelset_evolution_multiple.Advance_Levelset(dt);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after phi",1);
    LOG::Time("advecting particles");
    example.particle_levelset_evolution_multiple.particle_levelset_multiple.Euler_Step_Particles(face_velocities_ghost,dt,time,true,true,false);

    PHYSBAM_DEBUG_WRITE_SUBSTEP("after particle advection",1);

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
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after particles",1);

    LOG::Time("deleting particles"); // needs to be after adding sources, since it does a reseed
    example.particle_levelset_evolution_multiple.Delete_Particles_Outside_Grid();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete 1",1);
    LOG::Time("deleting particles in local maxima");
    for(int i=0;i<example.number_of_colors;i++)
        example.particle_levelset_evolution_multiple.Particle_Levelset(i).Delete_Particles_In_Local_Maximum_Phi_Cells(0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete 2",1);
    LOG::Time("deleting particles far from interface");
    for(int i=0;i<example.number_of_colors;i++)
        example.particle_levelset_evolution_multiple.Particle_Levelset(i).Delete_Particles_Far_From_Interface(); // uses visibility
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete 3",1);

    LOG::Time("re-incorporating removed particles");
    for(int i=0;i<example.number_of_colors;i++)
        example.particle_levelset_evolution_multiple.Particle_Levelset(i).Identify_And_Remove_Escaped_Particles(face_velocities_ghost,5,time+dt);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after remove",1);
    example.Rebuild_Levelset_Color();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after rebuilding level set again",1);
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Advance_One_Time_Step(bool first_step)
{
    Extrapolate_Velocity(example.face_velocities,example.face_color);
    if(example.use_polymer_stress) Extrapolate_Polymer_Stress(example.polymer_stress,example.cell_color);
    example.Begin_Time_Step(time);
    T dt=example.dt;
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before levelset evolution",1);

    if(example.use_level_set_method) Update_Level_Set(dt,first_step);
    else if(example.use_pls) Update_Pls(dt);

    PHYSBAM_DEBUG_WRITE_SUBSTEP("before velocity advection",1);
    Extrapolate_Velocity(example.face_velocities,example.face_color);
    if(example.use_polymer_stress) Extrapolate_Polymer_Stress(example.polymer_stress,example.cell_color);
    Advection_And_BDF(dt,first_step);
    Extrapolate_Velocity(example.face_velocities,example.face_color);
    if(example.use_polymer_stress) Extrapolate_Polymer_Stress(example.polymer_stress,example.cell_color);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before solve",1);
    Apply_Pressure_And_Viscosity(dt,first_step);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after solve",1);

    example.time=time+=dt;
    Extrapolate_Velocity(example.face_velocities,example.face_color);
    if(example.use_polymer_stress){
        Update_Polymer_Stress(dt);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after update polymer stress",1);}
    if(example.use_polymer_stress) Extrapolate_Polymer_Stress(example.polymer_stress,example.cell_color);
    example.End_Time_Step(time);
}
//#####################################################################
// Function Update_Level_Set
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Update_Level_Set(T dt,bool first_step)
{
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before phi advection",1);
    example.particle_levelset_evolution_multiple.Advance_Levelset(dt);
    example.Make_Levelsets_Consistent();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before reinitialization",1);
    example.particle_levelset_evolution_multiple.Make_Signed_Distance();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after level set update",1);
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
    if(example.use_polymer_stress) example.prev_polymer_stress.Exchange(example.polymer_stress);
    example.prev_face_color=example.face_color;
    example.prev_cell_color=example.cell_color;
}
//#####################################################################
// Function Advection_And_BDF
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
No_Advection_And_BDF(T dt,bool first_step,int c)
{
    if(!first_step) example.prev_face_velocities(c).Copy((T)2/(T)1.5,example.face_velocities(c),-(T).5/(T)1.5,example.prev_face_velocities(c));
    else example.prev_face_velocities(c)=example.face_velocities(c);

    if(example.use_polymer_stress) example.prev_polymer_stress(c)=example.polymer_stress(c);
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
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,SYMMETRIC_MATRIX<T,TV::m>,AVERAGING_UNIFORM<TV>,QUADRATIC_INTERPOLATION_UNIFORM<TV,SYMMETRIC_MATRIX<T,TV::m> > > quadratic_advection_polymer_stress;
    BOUNDARY_MAC_GRID_PERIODIC<TV,T> boundary;
    BOUNDARY_MAC_GRID_PERIODIC<TV,SYMMETRIC_MATRIX<T,TV::m> > boundary_polymer_stress;
    FACE_LOOKUP_UNIFORM<TV> lookup_face_velocities(example.face_velocities(c)),lookup_prev_face_velocities(example.prev_face_velocities(c));
    PHYSBAM_DEBUG_ONLY(Assert_Advection_CFL(example.face_velocities(c),example.face_color,c,dt));
    if(!first_step){
        PHYSBAM_DEBUG_ONLY(Assert_Advection_CFL(example.prev_face_velocities(c),example.prev_face_color,c,dt));
        ARRAY<T,FACE_INDEX<TV::m> > temp(example.grid,example.number_of_ghost_cells);
        quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,temp,lookup_prev_face_velocities,lookup_prev_face_velocities,boundary,2*dt,time+dt);
        quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,example.prev_face_velocities(c),lookup_face_velocities,lookup_face_velocities,boundary,dt,time+dt);
        example.prev_face_velocities(c).Copy((T)2/(T)1.5,example.prev_face_velocities(c),-(T).5/(T)1.5,temp);}
    else quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,example.prev_face_velocities(c),lookup_face_velocities,lookup_face_velocities,boundary,dt,time+dt);

    if(example.use_polymer_stress) quadratic_advection_polymer_stress.Update_Advection_Equation_Cell_Lookup(example.grid,example.prev_polymer_stress(c),example.polymer_stress(c),lookup_face_velocities,boundary_polymer_stress,dt,time+dt);
}
//#####################################################################
// Function Advection_And_BDF
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
RK2_Advection_And_BDF(T dt,bool first_step,int c)
{
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T,AVERAGING_UNIFORM<TV>,QUADRATIC_INTERPOLATION_UNIFORM<TV,T> > quadratic_advection;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,SYMMETRIC_MATRIX<T,TV::m>,AVERAGING_UNIFORM<TV>,QUADRATIC_INTERPOLATION_UNIFORM<TV,SYMMETRIC_MATRIX<T,TV::m> > > quadratic_advection_polymer_stress;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T,AVERAGING_UNIFORM<TV>,LINEAR_INTERPOLATION_UNIFORM<TV,T> > linear_advection;
    BOUNDARY_MAC_GRID_PERIODIC<TV,T> boundary;
    BOUNDARY_MAC_GRID_PERIODIC<TV,SYMMETRIC_MATRIX<T,TV::m> > boundary_polymer_stress;
    ARRAY<T,FACE_INDEX<TV::m> > temp(example.grid,example.number_of_ghost_cells),temp2(example.grid,example.number_of_ghost_cells),
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

    if(example.use_polymer_stress) quadratic_advection_polymer_stress.Update_Advection_Equation_Cell_Lookup(example.grid,example.prev_polymer_stress(c),example.polymer_stress(c),lookup_face_velocities,boundary_polymer_stress,dt,time+dt);
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
    if(example.use_polymer_stress){
        iss.polymer_stress_coefficient=example.polymer_stress_coefficient;
        iss.inv_Wi=example.inv_Wi;
        iss.stored_polymer_stress=&example.polymer_stress;}

    ARRAY<T> inertia=example.rho,dt_mu(example.mu*dt);
    if(!first_step) inertia*=(T)1.5;
    iss.Set_Matrix(dt_mu,example.use_discontinuous_velocity,
        [=](const TV& X,int color0,int color1){return example.Velocity_Jump(X,color0,color1,time+dt);},
        [=](const TV& X,int color0,int color1){return example.Jump_Interface_Condition(X,color0,color1,time+dt)*dt;},
        &inertia,true,dt);

    printf("\n");
    for(int i=0;i<TV::m;i++){for(int c=0;c<iss.cdi->colors;c++) printf("%c%d [%i]\t","uvw"[i],c,iss.cm_u(i)->dofs(c));printf("\n");}
    for(int c=0;c<iss.cdi->colors;c++) printf("p%d [%i]\t",c,iss.cm_p->dofs(c));
    printf("\n");
    printf("qn [%i]\t",iss.cdi->constraint_base_n);
    printf("qt [%i] ",iss.cdi->constraint_base_t);
    printf("\n");

    INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV> rhs,sol;
    iss.Set_RHS(rhs,[=](const TV& X,int color){return example.Volume_Force(X,color,time+dt)*dt;},&example.face_velocities,false);
    if(example.use_polymer_stress) iss.Add_Polymer_Stress_RHS(rhs,example.polymer_stress,dt);
    iss.Resize_Vector(sol);

    MINRES<T> mr;
    GMRES<T> gr;
    KRYLOV_SOLVER<T>* solver=&mr;
    if(example.use_multigrid) solver=&gr;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;

    if(example.test_system) iss.Test_System(sol);
    if(example.dump_matrix){
        OCTAVE_OUTPUT<T>(LOG::sprintf("M-%d.txt",solve_id).c_str()).Write("M",iss,rhs);
        OCTAVE_OUTPUT<T>(LOG::sprintf("Z-%d.txt",solve_id).c_str()).Write_Preconditioner("Z",iss,rhs);
        OCTAVE_OUTPUT<T>(LOG::sprintf("P-%d.txt",solve_id).c_str()).Write_Projection("P",iss,rhs);
        OCTAVE_OUTPUT<T>(LOG::sprintf("b-%d.txt",solve_id).c_str()).Write("b",rhs);}
    if(example.sparse_dump_matrix){
        SPARSE_MATRIX_FLAT_MXN<T> M;
        iss.Get_Sparse_Matrix(M);
        OCTAVE_OUTPUT<T>(LOG::sprintf("M-%d.txt",solve_id).c_str()).Write("M",M);
        OCTAVE_OUTPUT<T>(LOG::sprintf("b-%d.txt",solve_id).c_str()).Write("b",rhs);}
    solver->Solve(iss,sol,rhs,vectors,example.solver_tolerance,0,example.max_iter);

    if(example.dump_matrix){
        OCTAVE_OUTPUT<T>(LOG::sprintf("x-%d.txt",solve_id).c_str()).Write("x",sol);}
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
        int k=iss.cm_u(it.Axis())->Get_Index(it.face.index,c);
        assert(k>=0);
        example.face_velocities(c)(it.Full_Index())=sol.u(it.Axis())(c)(k);}

    for(CELL_ITERATOR<TV> it(example.grid);it.Valid();it.Next()){
        int c=example.levelset_color.Color(it.Location());
        example.cell_color(it.index)=c;}

    for(FACE_ITERATOR<TV> it(example.grid,example.number_of_ghost_cells,GRID<TV>::GHOST_REGION);it.Valid();it.Next())
        example.face_color(it.Full_Index())=example.levelset_color.Color(it.Location());
    for(CELL_ITERATOR<TV> it(example.grid,example.number_of_ghost_cells,GRID<TV>::GHOST_REGION);it.Valid();it.Next())
        example.cell_color(it.index)=example.levelset_color.Color(it.Location());
    vectors.Delete_Pointers_And_Clean_Memory();

    for(int c=0;c<example.face_velocities.m;c++){
        example.boundary.Apply_Boundary_Condition_Face(example.grid,example.face_velocities(c),time+example.dt);
        example.boundary.Fill_Ghost_Faces(example.grid,example.face_velocities(c),example.face_velocities(c),0,example.number_of_ghost_cells);}
    example.boundary_int.Apply_Boundary_Condition_Face(example.grid,example.face_color,time+example.dt);
    example.boundary_int.Fill_Ghost_Faces(example.grid,example.face_color,example.face_color,0,example.number_of_ghost_cells);
    example.boundary_int.Apply_Boundary_Condition(example.grid,example.cell_color,time+example.dt);
    example.boundary_int.Fill_Ghost_Cells(example.grid,example.cell_color,example.cell_color,0,example.number_of_ghost_cells);

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
    for(CELL_ITERATOR<TV> it(example.grid);it.Valid();it.Next()){
        int c=example.cell_color(it.index);
        if(c<0) continue;
        SYMMETRIC_MATRIX<T,TV::m>& S=example.polymer_stress(c)(it.index);
        MATRIX<T,TV::m> du=AVERAGING_UNIFORM<TV>::Cell_Centered_Gradient(example.grid,example.face_velocities(c),it.index);
        MATRIX<T,TV::m> A=du*dt+1;
        SYMMETRIC_MATRIX<T,TV::m> f=example.Polymer_Stress_Forcing_Term(it.Location(),c,example.time+example.dt/2);
        S=(SYMMETRIC_MATRIX<T,TV::m>::Conjugate(A,S)+dt*example.inv_Wi(c)+dt*f)/(1+dt*example.inv_Wi(c));}
}
//#####################################################################
// Function Extrapolate_Velocity
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Extrapolate_Velocity(ARRAY<T,FACE_INDEX<TV::m> >& u,const ARRAY<int,FACE_INDEX<TV::m> >& color,int c)
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
Extrapolate_Velocity(ARRAY<ARRAY<T,FACE_INDEX<TV::m> > >& S,const ARRAY<int,FACE_INDEX<TV::m> >& color)
{
    for(int c=0;c<example.number_of_colors;c++)
        Extrapolate_Velocity(S(c),color,c);
}
//#####################################################################
// Function Extrapolate_Polymer_Stress
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Extrapolate_Polymer_Stress(ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT>& S,const ARRAY<int,TV_INT>& color,int c)
{
    for(CELL_ITERATOR<TV> it(example.grid);it.Valid();it.Next())
        if(color(it.index)!=c)
            S(it.index)=SYMMETRIC_MATRIX<T,TV::m>()+1e20;
    for(CELL_ITERATOR<TV> it(example.grid,example.number_of_ghost_cells,GRID<TV>::GHOST_REGION);it.Valid();it.Next())
        S(it.index)=SYMMETRIC_MATRIX<T,TV::m>()+1e20;

    const LEVELSET<TV>& phi=*example.particle_levelset_evolution_multiple.particle_levelset_multiple.levelset_multiple.levelsets(c);
    EXTRAPOLATION_HIGHER_ORDER<TV,SYMMETRIC_MATRIX<T,TV::m> > eho(example.grid,phi,example.number_of_ghost_cells*10,3,example.number_of_ghost_cells);
    eho.periodic=true;

    eho.Extrapolate_Cell([&](const TV_INT& index){return color(index)==c;},S);

    example.boundary_symmetric.Apply_Boundary_Condition(example.grid,S,time+example.dt);
    example.boundary_symmetric.Fill_Ghost_Cells(example.grid,S,S,0,example.number_of_ghost_cells);
}
//#####################################################################
// Function Extrapolate_Polymer_Stress
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Extrapolate_Polymer_Stress(ARRAY<ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT> >& S,const ARRAY<int,TV_INT>& color)
{
    for(int c=0;c<example.number_of_colors;c++)
        Extrapolate_Polymer_Stress(S(c),color,c);
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
            DEBUG_SUBSTEPS::write_substeps_level=example.write_substeps_level;
        for(int substep=0;substep<example.time_steps_per_frame;substep++){
            LOG::SCOPE scope("SUBSTEP","substep %d",substep+1);
            Advance_One_Time_Step(current_frame==0 && substep==0);}
        if(current_frame%1==0 && 0){
            example.particle_levelset_evolution_multiple.Reseed_Particles(time);
            example.particle_levelset_evolution_multiple.Delete_Particles_Outside_Grid();}
        Write_Output_Files();}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Write_Substep(const std::string& title)
{
    example.frame_title=title;
    LOG::cout<<"Writing substep ["<<example.frame_title<<"]: output_number="<<example.viewer_dir.frame_stack<<", time="<<time<<", frame="<<current_frame<<std::endl;
    Write_Output_Files();
    example.frame_title="";
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Write_Output_Files()
{
    example.viewer_dir.Start_Directory(0,example.frame_title);
    example.frame_title="";
    example.Write_Output_Files();
    example.viewer_dir.Finish_Directory();
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

    OCTAVE_OUTPUT<T> oo(LOG::sprintf("ev-%d.txt",solve_id).c_str());

    for(int j=0;j<5;j++){
        random.Fill_Uniform(sol.u,-1,1);
        for(int i=0;i<evs.m;i++) sol.Copy(-iss.Inner_Product(sol,*evs(i)),*evs(i),sol);
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
            for(int i=0;i<evs.m;i++) sol.Copy(-iss.Inner_Product(sol,*evs(i)),*evs(i),sol);

            a=sqrt((T)iss.Inner_Product(sol,sol));
            sol*=1/a;
            LOG::cout<<"eig approx "<<a<<"   dp "<<iss.Inner_Product(sol,tmp)<<std::endl;
            tmp=sol;}

        Dump_Vector(iss,sol,LOG::sprintf("eigenvector (%g)",a).c_str());
        evs.Append(new INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>(sol));
        oo.Write(LOG::sprintf("EV%d",j).c_str(),sol);
        oo.Write(LOG::sprintf("ev%d",j).c_str(),a);}

    evs.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Dump_Vector
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Dump_Vector(const INTERFACE_STOKES_SYSTEM_COLOR<TV>& iss,const KRYLOV_VECTOR_BASE<T>& u,const char* str) const
{
    const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>& v=static_cast<const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>&>(u);
    ARRAY<ARRAY<T,FACE_INDEX<TV::m> > > fv=example.face_velocities;
    for(int i=0;i<example.face_velocities.m;i++)
        for(int j=0;j<TV::m;j++)
            example.face_velocities(i).Component(j).Fill(0);

    for(FACE_ITERATOR<TV> it(example.grid);it.Valid();it.Next()){
        int c=example.face_color(it.Full_Index());
        if(c<0) continue;
        int k=iss.cm_u(it.Axis())->Get_Index(it.face.index,c);
        assert(k>=0);
        example.face_velocities(c)(it.Full_Index())=v.u(it.Axis())(c)(k);}

    for(CELL_ITERATOR<TV> it(example.grid);it.Valid();it.Next()){
        int c=example.levelset_color.Color(it.Location());
        if(c<0) continue;
        int k=iss.cm_p->Get_Index(it.index,c);
        assert(k>=0);
        T p=v.p(c)(k);
        Add_Debug_Particle(it.Location(),VECTOR<T,3>(p<0,p>=0,0));
        Debug_Particle_Set_Attribute<TV>("display_size",abs(p));}

    PHYSBAM_DEBUG_WRITE_SUBSTEP(str,1);
    example.face_velocities=fv;
}
//#####################################################################
namespace PhysBAM{
    template class PLS_FC_DRIVER<VECTOR<float,2> >;
    template class PLS_FC_DRIVER<VECTOR<float,3> >;
    template class PLS_FC_DRIVER<VECTOR<double,2> >;
    template class PLS_FC_DRIVER<VECTOR<double,3> >;
}
