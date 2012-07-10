//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
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

    example.levelset_color.phi.Resize(example.grid.Node_Indices());
    example.levelset_color.color.Resize(example.grid.Node_Indices());
    example.face_velocities.Resize(example.grid);
    example.face_color.Resize(example.grid);
    {
        example.particle_levelset_evolution.Initialize_Domain(example.grid);
        example.particle_levelset_evolution.particle_levelset.Set_Band_Width(6);
        example.collision_bodies_affecting_fluid.Initialize_Grids();
    }

    example.particle_levelset_evolution.Set_Time(time);
    example.particle_levelset_evolution.Set_CFL_Number((T).9);

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
    example.particle_levelset_evolution.Use_Particle_Levelset(true);
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
        example.particle_levelset_evolution.Make_Signed_Distance();
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time);}

    example.collision_bodies_affecting_fluid.Compute_Grid_Visibility();
    example.particle_levelset_evolution.Set_Seed(2606);
    if(!example.restart) example.particle_levelset_evolution.Seed_Particles(time);
    example.particle_levelset_evolution.Delete_Particles_Outside_Grid();

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
    T maximum_fluid_speed=example.face_velocities.Maxabs().Max();
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
    example.boundary->Fill_Ghost_Cells_Face(example.grid,example.face_velocities,face_velocities_ghost,time+dt,example.number_of_ghost_cells);

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
    T dt=example.dt;
//    Update_Pls(dt);

    Apply_Pressure_And_Viscosity(dt);
    time+=dt;
}
//#####################################################################
// Function Apply_Pressure_And_Viscosity
//#####################################################################
template<class TV> void PLS_FC_DRIVER<TV>::
Apply_Pressure_And_Viscosity(T dt)
{
    struct BOUNDARY_CONDITIONS_COLOR_LOCAL:public BOUNDARY_CONDITIONS_COLOR<TV>
    {
        virtual TV f_surface(const TV& X,int color0,int color1) {return TV();}
    } bccl;

    ARRAY<T,TV_INT> phi_value(example.grid.Node_Indices());
    ARRAY<int,TV_INT> phi_color(example.grid.Node_Indices());

    INTERFACE_STOKES_SYSTEM_COLOR<TV> iss(example.grid,example.levelset_color.phi,example.levelset_color.color);
    iss.use_preconditioner=example.use_preconditioner;
    iss.Set_Matrix(ARRAY<T>(example.mu*dt),example.wrap,&bccl,0);

    printf("\n");
    for(int i=0;i<TV::m;i++){for(int c=0;c<iss.cdi->colors;c++) printf("%c%d [%i]\t","uvw"[i],c,iss.cm_u(i)->dofs(c));printf("\n");}
    for(int c=0;c<iss.cdi->colors;c++) printf("p%d [%i]\t",c,iss.cm_p->dofs(c));printf("\n");
    printf("qn [%i]\t",iss.cdi->constraint_base_n);
    printf("qt [%i] ",iss.cdi->constraint_base_t);
    printf("\n");

    INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV> rhs,sol;
    {
        ARRAY<ARRAY<TV,TV_INT> > f_volume;
        ARRAY<ARRAY<T,FACE_INDEX<TV::m> > > u;
        ARRAY<ARRAY<bool,FACE_INDEX<TV::m> > > inside;
        
        f_volume.Resize(iss.cdi->colors);
        u.Resize(iss.cdi->colors);
        inside.Resize(iss.cdi->colors);

        for(int c=0;c<iss.cdi->colors;c++){
            f_volume(c).Resize(example.grid.Domain_Indices());
            u(c).Resize(example.grid,2);
            inside(c).Resize(example.grid);}

        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(example.grid);it.Valid();it.Next()){
            int c=example.face_color(it.Full_Index());
            u(c)(it.Full_Index())=example.face_velocities(it.Full_Index());
            inside(c)(it.Full_Index())=true;}

        for(int c=0;c<iss.cdi->colors;c++){
            EXTRAPOLATION_HIGHER_ORDER_POLY<TV,T>::Extrapolate_Face(example.grid,inside(c),2,u(c),3,2);
            f_volume(c).Resize(example.grid.Domain_Indices());
            u(c).Resize(example.grid);
            inside(c).Resize(example.grid);}

        iss.Set_RHS(rhs,f_volume,&u);
        iss.Resize_Vector(sol);
    }

    MINRES<T> mr;
    KRYLOV_SOLVER<T>* solver=&mr;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;

    if(example.dump_matrix) OCTAVE_OUTPUT<T>("M.txt").Write("M",iss,*vectors(0),*vectors(1));
    solver->Solve(iss,sol,rhs,vectors,1e-10,0,example.max_iter);

    iss.Multiply(sol,*vectors(0));
    *vectors(0)-=rhs;
    LOG::cout<<"Residual: "<<iss.Convergence_Norm(*vectors(0))<<std::endl;

    for(int i=0;i<iss.null_modes.m;i++){
        iss.Multiply(iss.null_modes(i),*vectors(0));
        LOG::cout<<"null mode["<<i<<"] "<<iss.Convergence_Norm(*vectors(0))<<std::endl;}

    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(example.grid);it.Valid();it.Next()){
        int c=example.levelset_color.Color(it.Location());
        example.face_color(it.Full_Index())=c;
        int k=iss.cm_u(it.Axis())->Get_Index(it.index,c);
        assert(k>=0);
        example.face_velocities(it.Full_Index())=sol.u(it.Axis())(c)(k);}
    vectors.Delete_Pointers_And_Clean_Memory();
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
        if(current_frame%1==0){
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
