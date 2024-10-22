//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Interpolation/AVERAGING_UNIFORM.h>
#include <Incompressible/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <Dynamics/Heat_Flows/HEAT_LAPLACE.h>
#include <Dynamics/Incompressible_Flows/IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM.h>
#include <Dynamics/Incompressible_Flows/INCOMPRESSIBLE_MULTIPHASE_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<TV>::
IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM(PROJECTION_DYNAMICS_UNIFORM<TV>& projection_input,const ARRAY<T,TV_INT>& variable_viscosity_input,const ARRAY<T>& densities_input,const ARRAY<T>& viscosities_input,T_MPI_GRID* mpi_grid_input,const int axis_input,bool use_variable_viscosity_input)
    :IMPLICIT_VISCOSITY_UNIFORM<TV>(*projection_input.elliptic_solver,variable_viscosity_input,(T)0,(T)0,mpi_grid_input,axis_input,use_variable_viscosity_input,false),projection(projection_input),densities(densities_input),viscosities(viscosities_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<TV>::
~IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM()
{}
//#####################################################################
// Function Allocate_Heat_Solver
//#####################################################################
template<class TV> void IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<TV>::
Allocate_Heat_Solver()
{
    heat_solver=new HEAT_LAPLACE<POISSON_COLLIDABLE_UNIFORM<TV> >(face_grid,u);
}
//#####################################################################
// Function Setup_Viscosity
//#####################################################################
template<class TV> void IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<TV>::
Setup_Viscosity(const T dt)
{
    if(use_variable_viscosity) PHYSBAM_NOT_IMPLEMENTED();

    POISSON_COLLIDABLE_UNIFORM<TV>& heat_poisson=dynamic_cast<POISSON_COLLIDABLE_UNIFORM<TV>&>(*heat_solver);
    PROJECTION_DYNAMICS_UNIFORM<TV>& projection_dynamics=dynamic_cast<PROJECTION_DYNAMICS_UNIFORM<TV>&>(projection);
    heat_poisson.multiphase=true;
    const int number_of_regions=densities.m;

    // set viscosity coefficients
    ARRAY<T> dt_times_kinematic_viscosity(number_of_regions);
    for(int i=0;i<number_of_regions;i++) dt_times_kinematic_viscosity(i)=dt*viscosities(i)/densities(i);
    heat_poisson.Set_Constant_beta(dt_times_kinematic_viscosity);

    // set up internal levelset
    heat_poisson.Use_Internal_Level_Set(number_of_regions);
    T_AVERAGING averaging;const LEVELSET_MULTIPLE<TV>& cell_centered_levelset_multiple=*projection_dynamics.poisson_collidable->levelset_multiple;
    for(CELL_ITERATOR<TV> iterator(face_grid,2);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index(),p_face_index=cell_index;
        for(int i=0;i<number_of_regions;i++) heat_poisson.levelset_multiple->phis(i)(cell_index)=averaging.Cell_To_Face(projection.p_grid,axis,p_face_index,cell_centered_levelset_multiple.phis(i));}
    heat_poisson.levelset_multiple->Project_Levelset(2);

    if(projection_dynamics.flame) Calculate_Velocity_Jump();
    heat_poisson.Find_Constant_beta_Multiphase(heat_poisson.levelset_multiple->phis);
}
//#####################################################################
// Function Setup_Boundary_Conditions
//#####################################################################
template<class TV> void IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<TV>::
Setup_Boundary_Conditions(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities)
{
    IMPLICIT_VISCOSITY_UNIFORM<TV>::Setup_Boundary_Conditions(face_velocities);
    // set neumann b.c. at zero viscosity faces
    POISSON_COLLIDABLE_UNIFORM<TV>& heat_poisson=dynamic_cast<POISSON_COLLIDABLE_UNIFORM<TV>&>(*heat_solver);
    for(FACE_ITERATOR<TV> iterator(face_grid);iterator.Valid();iterator.Next()){int face_axis=iterator.Axis();TV_INT face=iterator.Face_Index();
        if(!heat_poisson.beta_face(face_axis,face)) heat_poisson.psi_N(face_axis,face)=true;}
}
//#####################################################################
// Function Calculate_Velocity_Jump
//#####################################################################
template<class TV> void IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<TV>::
Calculate_Velocity_Jump()
{
    POISSON_COLLIDABLE_UNIFORM<TV>& heat_poisson=dynamic_cast<POISSON_COLLIDABLE_UNIFORM<TV>&>(*heat_solver);
    PROJECTION_DYNAMICS_UNIFORM<TV>& projection_dynamics=dynamic_cast<PROJECTION_DYNAMICS_UNIFORM<TV>&>(projection);
    heat_poisson.Set_Jump_Multiphase();
    const ARRAY<TRIPLE<T,T,T> ,VECTOR<int,2> >& flame_speed_constants=projection_dynamics.flame_speed_constants;
    TV_INT axis_offset=TV_INT::Axis_Vector(axis);
    for(FACE_ITERATOR<TV> iterator(face_grid);iterator.Valid();iterator.Next()){TV_INT face_index=iterator.Face_Index();TV_INT face_axis_offset=TV_INT::Axis_Vector(iterator.Axis());
        int region_1=heat_poisson.levelset_multiple->Inside_Region(iterator.First_Cell_Index()),region_2=heat_poisson.levelset_multiple->Inside_Region(iterator.Second_Cell_Index());
        // [Vn]=M*[1/density] with M=-density_fuel*flame_speed, [1/density]=(1/density_fuel-1/density_products), flame_speed_constant.z=(-density_fuel*[1/density])
        int fuel_region=region_1,product_region=region_2;
        if(densities(region_1)<densities(region_2)){fuel_region=region_2;product_region=region_1;}
        const TRIPLE<T,T,T>& constants=flame_speed_constants(fuel_region,product_region);if(constants.z==0)continue;
        const LEVELSET<TV>& levelset=*projection_dynamics.poisson_collidable->levelset_multiple->levelsets(fuel_region);
        T flame_speed=constants.x;
        if(constants.y){T face_curvature;TV_INT p_face_index=face_index;
            if(iterator.Axis()==axis)face_curvature=(*levelset.curvature)(p_face_index-axis_offset);
            else{face_curvature=(T).25*((*levelset.curvature)(p_face_index-axis_offset-face_axis_offset)+(*levelset.curvature)(p_face_index-face_axis_offset)+
                (*levelset.curvature)(p_face_index)+(*levelset.curvature)(p_face_index-axis_offset));}
            flame_speed+=constants.y*face_curvature;}
        T face_normal;TV_INT p_face_index=face_index;
        if(iterator.Axis()==axis)face_normal=(*levelset.normals)(p_face_index-axis_offset)[axis];
        else{face_normal=((*levelset.normals)(p_face_index-axis_offset-face_axis_offset)+(*levelset.normals)(p_face_index-face_axis_offset)+
            (*levelset.normals)(p_face_index)+(*levelset.normals)(p_face_index-axis_offset)).Normalized()[axis];}
        heat_poisson.u_jump_face.Component(iterator.Axis())(face_index)=LEVELSET_MULTIPLE<TV>::Sign(product_region,fuel_region)*constants.z*flame_speed*face_normal;}
}
//#####################################################################
// Debug_Write
//#####################################################################
template<class TV> void IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<TV>::
Debug_Write(const VIEWER_DIR& viewer_dir)
{
    POISSON_COLLIDABLE_UNIFORM<TV>& heat_poisson=dynamic_cast<POISSON_COLLIDABLE_UNIFORM<TV>&>(*heat_solver);
    static int frame[3]={0,0,0};
    // std::string viewer_dir.output_directory=output_directory_input;
    // if(mpi_grid) viewer_dir.output_directory+=LOG::sprintf("/processor%d",mpi_grid->rank);
    std::string output_directory_axis=LOG::sprintf("%s/%d",viewer_dir.current_directory,axis);
    Create_Directory(output_directory_axis);
    std::string f=Value_To_String(frame[axis]);
    Write_To_File<T>(output_directory_axis+"/grid",face_grid);
    Write_To_File<T>(output_directory_axis+"/psi_N."+f,heat_poisson.psi_N);
    Write_To_File<T>(output_directory_axis+"/psi_D."+f,heat_poisson.psi_D);
    Write_To_File<T>(output_directory_axis+"/colors."+f,heat_poisson.filled_region_colors);
    Write_To_File<T>(output_directory_axis+"/beta_face."+f,heat_poisson.beta_face);
    for(int i=0;i<densities.m;i++){
        std::string filename=LOG::sprintf("/levelset_%d.%s",i,f.c_str());
        Write_To_File<T>(output_directory_axis+filename,*heat_poisson.levelset_multiple->levelsets(i));}
}
//#####################################################################
namespace PhysBAM{
template class IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<VECTOR<float,1> >;
template class IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<VECTOR<float,2> >;
template class IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<VECTOR<float,3> >;
template class IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<VECTOR<double,1> >;
template class IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<VECTOR<double,2> >;
template class IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<VECTOR<double,3> >;
}
