//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Frank Losasso, Duc Nguyen, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Advection/ADVECTION.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Grids_Uniform_Computations/VORTICITY_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/maxabs.h>
#include <Tools/Math_Tools/min.h>
#include <Tools/Math_Tools/minmag.h>
#include <Tools/Math_Tools/sqr.h>
#include <Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <Dynamics/Incompressible_Flows/IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM.h>
#include <Dynamics/Incompressible_Flows/INCOMPRESSIBLE_MULTIPHASE_UNIFORM.h>
#include <Dynamics/Interpolation/FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM.h>
using namespace PhysBAM;
#ifdef _WIN32
#pragma warning(disable:4723)
#endif
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
INCOMPRESSIBLE_MULTIPHASE_UNIFORM(const T_GRID& grid_input,PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection_input)
    :INCOMPRESSIBLE_UNIFORM<T_GRID>(grid_input,projection_input),levelset_for_dirichlet_regions(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
~INCOMPRESSIBLE_MULTIPHASE_UNIFORM()
{
}
//#####################################################################
// Function Advance_One_Time_Step_Explicit_Part
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
Advance_One_Time_Step_Convection(const T dt,const T time,T_FACE_ARRAYS_SCALAR& advecting_face_velocities,T_FACE_ARRAYS_SCALAR& face_velocities_to_advect,const ARRAY<bool>* pseudo_dirichlet_regions,const int number_of_ghost_cells)
{
    // TODO: make efficient if advection velocities are same as advected velocities
    T_FACE_ARRAYS_SCALAR advection_face_velocities_ghost;advection_face_velocities_ghost.Resize(grid,number_of_ghost_cells,false);
    boundary->Fill_Ghost_Faces(grid,advecting_face_velocities,advection_face_velocities_ghost,time,number_of_ghost_cells);
    T_FACE_ARRAYS_SCALAR face_velocities_to_advect_ghost;face_velocities_to_advect_ghost.Resize(grid,number_of_ghost_cells,false);
    boundary->Fill_Ghost_Faces(grid,face_velocities_to_advect,face_velocities_to_advect_ghost,time,number_of_ghost_cells);

    // update convection
    if(pseudo_dirichlet_regions->Number_True()>0){
        T_FACE_ARRAYS_SCALAR face_velocities_liquid=face_velocities_to_advect;T_FACE_ARRAYS_SCALAR advection_face_velocities_ghost_extrapolated=advection_face_velocities_ghost;
        T_FACE_ARRAYS_SCALAR face_velocities_to_advect_ghost_extrapolated=face_velocities_to_advect_ghost;
        T_ARRAYS_SCALAR phi_for_pseudo_dirichlet_regions;T_GRID grid_temp(grid);
        LEVELSET<TV> levelset_for_pseudo_dirichlet_regions(grid_temp,phi_for_pseudo_dirichlet_regions);
        projection.poisson_collidable->levelset_multiple->Get_Single_Levelset(*pseudo_dirichlet_regions,levelset_for_pseudo_dirichlet_regions,false);
        Extrapolate_Velocity_Across_Interface(advection_face_velocities_ghost_extrapolated,phi_for_pseudo_dirichlet_regions);
        Extrapolate_Velocity_Across_Interface(face_velocities_to_advect_ghost_extrapolated,phi_for_pseudo_dirichlet_regions);
        advection->Update_Advection_Equation_Face(grid,face_velocities_liquid,face_velocities_to_advect_ghost_extrapolated,advection_face_velocities_ghost_extrapolated,*boundary,dt,time);
        advection->Update_Advection_Equation_Face(grid,face_velocities_to_advect,face_velocities_to_advect_ghost,advection_face_velocities_ghost,*boundary,dt,time);
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){const int axis=iterator.Axis();const TV_INT face_index=iterator.Face_Index();
            int region1=projection.poisson_collidable->levelset_multiple->Inside_Region(iterator.First_Cell_Index());
            int region2=projection.poisson_collidable->levelset_multiple->Inside_Region(iterator.Second_Cell_Index());
            if(!(*pseudo_dirichlet_regions)(region1)||!(*pseudo_dirichlet_regions)(region2))
                face_velocities_to_advect.Component(axis)(face_index)=face_velocities_liquid.Component(axis)(face_index);}}
    else advection->Update_Advection_Equation_Face(grid,face_velocities_to_advect,face_velocities_to_advect_ghost,advection_face_velocities_ghost,*boundary,dt,time);
}
//#####################################################################
// Function Advance_One_Time_Step_Explicit_Part
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
Advance_One_Time_Step_Forces(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const bool implicit_viscosity,const ARRAY<T_ARRAYS_SCALAR>* phi_ghost,const ARRAY<bool>* pseudo_dirichlet_regions,const int number_of_ghost_cells)
{
    T_FACE_ARRAYS_SCALAR face_velocities_ghost(grid.Domain_Indices(number_of_ghost_cells),false);
    boundary->Fill_Ghost_Faces(grid,face_velocities,face_velocities_ghost,time,number_of_ghost_cells);

    // update strain and apply elastic forces
    for(int i=0;i<strains.m;i++)if(strains(i)){
        // extrapolate the velocity across the interface to get better strain boundaries
        T_ARRAYS_SCALAR phi_ghost(grid.Domain_Indices(number_of_ghost_cells));projection.poisson_collidable->levelset_multiple->levelsets(i)->boundary->Fill_Ghost_Cells(grid,projection.poisson_collidable->levelset_multiple->phis(i),phi_ghost,dt,time,number_of_ghost_cells);
        T_FACE_ARRAYS_SCALAR face_velocities_temp=face_velocities_ghost;
        for(int axis=0;axis<T_GRID::dimension;axis++){
            T_GRID face_grid=grid.Get_Face_Grid(axis);T_ARRAYS_SCALAR phi_face(face_grid.Domain_Indices(),false);T_ARRAYS_BASE& face_velocity=face_velocities_temp.Component(axis);
            for(FACE_ITERATOR<TV> iterator(grid,0,T_GRID::WHOLE_REGION,-1,axis);iterator.Valid();iterator.Next())
                phi_face(iterator.Face_Index())=(T).5*(phi_ghost(iterator.First_Cell_Index())+phi_ghost(iterator.Second_Cell_Index()));
            int extrapolation_bandwidth=3;
            EXTRAPOLATION_UNIFORM<GRID<TV>,T>  extrapolate(face_grid,phi_face,face_velocity,number_of_ghost_cells);extrapolate.Set_Band_Width((T)extrapolation_bandwidth);
            extrapolate.Extrapolate();}
        strains(i)->Update_Strain_Equation_Multiphase(dt,time,projection.densities(i),face_velocities,face_velocities_temp,*projection.poisson_collidable->levelset_multiple,i,number_of_ghost_cells);}

    // update gravity
    if(gravity) for(int axis=0;axis<T_GRID::dimension;axis++) if(downward_direction[axis]) face_velocities.Component(axis)+=dt*gravity*downward_direction[axis];
    
    // update body force
    if(use_force) for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) 
        face_velocities.Component(iterator.Axis())(iterator.Face_Index())+=dt*force.Component(iterator.Axis())(iterator.Face_Index());

    if((viscosity || use_variable_viscosity) && (!implicit_viscosity || use_explicit_part_of_implicit_viscosity))
        Discretize_Explicit_Viscous_Terms(dt);

    if(!GFM && nonzero_surface_tension){
        T half_width=(T).5*number_of_interface_cells*grid.dX.Min();T dt_over_four=(T).25*dt;
        LEVELSET_MULTIPLE<T_GRID>& levelset_multiple=*projection.poisson_collidable->levelset_multiple;
        levelset_multiple.Compute_Normals();levelset_multiple.Compute_Curvature();
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT index=iterator.Face_Index();
            TV_INT cell_1=iterator.First_Cell_Index(),cell_2=iterator.Second_Cell_Index();
            int region_1,region_2;T phi_1,phi_2;levelset_multiple.Minimum_Regions(cell_1,cell_2,region_1,region_2,phi_1,phi_2);
            T sign=(region_1==region_2)?(T)1:(T)-1;
            T phi_face=(T).5*(phi_1+sign*phi_2);
            T twice_curvature=((*levelset_multiple.levelsets(region_1)->curvature)(cell_1)+sign*(*levelset_multiple.levelsets(region_2)->curvature)(cell_2));
            T twice_face_normal_component=((*levelset_multiple.levelsets(region_1)->normals)(cell_1)[axis]+sign*(*levelset_multiple.levelsets(region_2)->normals)(cell_2)[axis]);
            face_velocities.Component(axis)(index)+=dt_over_four*LEVELSET_UTILITIES<T>::Delta(phi_face,half_width)*surface_tensions(region_1,region_2)*
                twice_curvature*twice_face_normal_component/LEVELSET_UTILITIES<T>::Heaviside(phi_face,projection.densities(region_1),projection.densities(region_2),half_width);}}

    if(use_variable_vorticity_confinement){
        ARRAY<TV,TV_INT> F(grid.Cell_Indices(1),false);
        Compute_Vorticity_Confinement_Force(grid,face_velocities_ghost,F);
        F*=dt*(T).5;F*=variable_vorticity_confinement;
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
            face_velocities.Component(axis)(iterator.Face_Index())+=F(iterator.First_Cell_Index())[axis]+F(iterator.Second_Cell_Index())[axis];}}

    if(vorticity_confinements.Count_Matches(0)!=vorticity_confinements.m){
        ARRAY<TV,TV_INT> F(grid.Cell_Indices(1),false);
        Compute_Vorticity_Confinement_Force(grid,face_velocities_ghost,F);
        T half_dt=(T).5*dt;
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();
            if(pseudo_dirichlet_regions){
                int region1=projection.poisson_collidable->levelset_multiple->Inside_Region(iterator.First_Cell_Index());
                int region2=projection.poisson_collidable->levelset_multiple->Inside_Region(iterator.Second_Cell_Index());
                if(!(*pseudo_dirichlet_regions)(region1) || !(*pseudo_dirichlet_regions)(region2))
                    if((!(*pseudo_dirichlet_regions)(region1) && vorticity_confinements(region1)==0) || (!(*pseudo_dirichlet_regions)(region2) && vorticity_confinements(region2)==0))
                        continue;}
            int region=projection.poisson_collidable->levelset_multiple->Inside_Region_Face(iterator.Axis(),iterator.Face_Index());
            T face_F=F(iterator.First_Cell_Index())[axis]+F(iterator.Second_Cell_Index())[axis];
            face_velocities.Component(axis)(iterator.Face_Index())+=face_F*half_dt*vorticity_confinements(region);}}

    boundary->Apply_Boundary_Condition_Face(grid,face_velocities,time+dt);
}
//#####################################################################
// Function Advance_One_Time_Step_Implicit_Part
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
Advance_One_Time_Step_Implicit_Part(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const bool implicit_viscosity)
{
    // boundary conditions
    boundary->Apply_Boundary_Condition_Face(grid,face_velocities,time+dt);

    // set up poisson equation
    ARRAY<T> one_over_densities(projection.densities.m);for(int i=0;i<projection.densities.m;i++)one_over_densities(i)=1/projection.densities(i);
    projection.poisson->Set_Constant_beta(one_over_densities);
    if(!GFM){projection.poisson->Use_Delta_Function_Method(number_of_interface_cells);projection.poisson->Smear_One_Over_beta();}

    if((GFM&&nonzero_surface_tension)||projection.flame) projection.poisson_collidable->Set_Jump_Multiphase();
    
    if(GFM && nonzero_surface_tension){
        LEVELSET_MULTIPLE<T_GRID>& levelset_multiple=*projection.poisson_collidable->levelset_multiple;levelset_multiple.Compute_Curvature();
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT cell_1=iterator.First_Cell_Index(),cell_2=iterator.Second_Cell_Index();
            int region_1,region_2;T phi_1,phi_2;levelset_multiple.Minimum_Regions(cell_1,cell_2,region_1,region_2,phi_1,phi_2);
            if(region_1!=region_2){
                T sign=LEVELSET_MULTIPLE<T_GRID>::Sign(region_1,region_2);
                T curvature=LEVELSET_UTILITIES<T>::Average(phi_1,-sign*(*levelset_multiple.levelsets(region_1)->curvature)(cell_1),
                                                       -phi_2,sign*(*levelset_multiple.levelsets(region_2)->curvature)(cell_2));
                projection.poisson_collidable->u_jump_face.Component(iterator.Axis())(iterator.Face_Index())+=dt*surface_tensions(region_1,region_2)*curvature;}}}

    if(projection.flame) Calculate_Pressure_Jump(dt,time);

    // viscosity
    if(nonzero_viscosity && implicit_viscosity){
        projection.Make_Divergence_Free(face_velocities,dt,time);Implicit_Viscous_Update(face_velocities,dt,time);}

    // make divergence free
    projection.Make_Divergence_Free(face_velocities,dt,time);
}
//#####################################################################
// Function Calculate_Pressure_Jump
//#####################################################################
// flame_speed must be up to date
template<class T_GRID> void INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
Calculate_Pressure_Jump(const T dt,const T time)
{
    assert(projection.poisson_collidable->levelset_multiple);
    LEVELSET_MULTIPLE<T_GRID>& levelset_multiple=*projection.poisson_collidable->levelset_multiple;
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell1=iterator.First_Cell_Index(),cell2=iterator.Second_Cell_Index();
        int region1=levelset_multiple.Inside_Region(cell1),region2=levelset_multiple.Inside_Region(cell2);
        const TRIPLE<T,T,T>& constants=projection.flame_speed_constants(region1,region2);
        if(constants.z==0)continue;
        if(projection.densities(region1)<projection.densities(region2))exchange(region1,region2); //region1 is now the fuel region
        // [p]=dt*[1/density]*sqr(M) with M=-density_fuel*flame_speed, [1/density]=(1/density_fuel-1/density_products)
        // flame_speed_constant.z is (-density_fuel*[1/density])
        projection.poisson_collidable->u_jump_face.Component(iterator.Axis())(iterator.Face_Index())+=LEVELSET_MULTIPLE<T_GRID>::Sign(region1,region2)*
            dt*constants.z*projection.densities(region1)*sqr(projection.Flame_Speed_Face_Multiphase(iterator.Axis(),iterator.Face_Index(),region1,region2));}
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
CFL(T_FACE_ARRAYS_SCALAR& face_velocities,const bool inviscid,const bool viscous_only) const
{
    TV DX=grid.dX,sqr_DX=DX*DX,max_abs_V;
    // convection
    T dt_convection=0;
    if(!viscous_only){
        if(projection.flame){
            FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<T_GRID> face_velocities_fire(face_velocities,projection,projection.poisson_collidable->levelset_multiple);
            typename FIRE_INTERPOLATION_POLICY<T_GRID>::AVERAGING_FIRE_MULTIPHASE averaging;
            for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
                TV V=averaging.Face_To_Face_Vector(grid,iterator.Axis(),iterator.Face_Index(),face_velocities_fire);
                dt_convection=max(dt_convection,TV::Dot_Product(grid.one_over_dX,abs(V)));}}
        else{
            for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
                T local_V_norm=0;for(int axis=0;axis<T_GRID::dimension;axis++){
                    TV_INT first_face_index=grid.First_Face_Index_In_Cell(axis,cell),second_face_index=grid.Second_Face_Index_In_Cell(axis,cell);
                    local_V_norm+=grid.one_over_dX[axis]*maxabs(face_velocities(axis,first_face_index),face_velocities(axis,second_face_index));}
                dt_convection=max(dt_convection,local_V_norm);}}}
    // viscosity
    T dt_viscosity=0;
    if(!inviscid){
        T norm_2_over_sqr_DX=2*Inverse(sqr_DX).Sum_Abs();
        for(int i=0;i<viscosities.m;i++)dt_viscosity=max(dt_viscosity,viscosities(i)/projection.densities(i));
        dt_viscosity*=norm_2_over_sqr_DX;
        if(use_variable_viscosity) PHYSBAM_NOT_IMPLEMENTED();}
    // surface tension
    T dt_surface_tension=0;
    if(nonzero_surface_tension){
        LEVELSET_MULTIPLE<T_GRID>& levelset_multiple=*projection.poisson_collidable->levelset_multiple;levelset_multiple.Compute_Curvature();
        int ghost_cells=1;
        levelset_multiple.Fill_Ghost_Cells(levelset_multiple.phis,0,ghost_cells);T kappa_cfl=0;
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT cell_1=iterator.First_Cell_Index(),cell_2=iterator.Second_Cell_Index();
            int region_1,region_2;T phi_1,phi_2;levelset_multiple.Minimum_Regions(cell_1,cell_2,region_1,region_2,phi_1,phi_2);
            if(region_1!=region_2 && surface_tensions(region_1,region_2)){
                T curvature=LEVELSET_UTILITIES<T>::Average(phi_1,LEVELSET_MULTIPLE<T_GRID>::Sign(region_2,region_1)*(*levelset_multiple.levelsets(region_1)->curvature)(cell_1),
                    -phi_2,LEVELSET_MULTIPLE<T_GRID>::Sign(region_1,region_2)*(*levelset_multiple.levelsets(region_2)->curvature)(cell_2));
                kappa_cfl=max(kappa_cfl,abs(curvature*surface_tensions(region_1,region_2)/
                          LEVELSET_UTILITIES<T>::Heaviside((T).5*(phi_1-phi_2),projection.densities(region_1),projection.densities(region_2))));}}
        dt_surface_tension=sqrt(kappa_cfl)/grid.dX.Min();}
    TV max_force;
    if(use_force) max_force=force.Max_Abs();
    T dt_force=0;
    if(use_force) dt_force=(max_force/DX).Sum_Abs();
    if(gravity) dt_force+=abs(gravity)*(downward_direction/DX).Sum_Abs();
    T strain_cfl=FLT_MAX;
    if(strain){
        for(int i=0;i<strains.m;i++)if(strains(i))
            strain_cfl=min(strain_cfl,strains(i)->CFL(projection.densities(i)));}
    dt_force+=1/strain_cfl;
    T dt_overall=(dt_convection+dt_viscosity+sqrt(sqr(dt_convection+dt_viscosity)+4*dt_force+4*sqr(dt_surface_tension)))/2; 
    return 1/max(dt_overall,1/max_time_step);
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
Set_Dirichlet_Boundary_Conditions(ARRAY<T_ARRAYS_SCALAR>& phis,const ARRAY<bool>& dirichlet_regions,const ARRAY<T>* pressures)
{
    LEVELSET_MULTIPLE<T_GRID> levelset_multiple(grid,phis);
    if(dirichlet_regions.Number_True()>0){
        if(!pressures){for(CELL_ITERATOR<TV> iterator(projection.p_grid);iterator.Valid();iterator.Next()) if(dirichlet_regions(levelset_multiple.Inside_Region(iterator.Cell_Index()))){
            projection.elliptic_solver->psi_D(iterator.Cell_Index())=true;projection.p(iterator.Cell_Index())=0;}}
        else{for(CELL_ITERATOR<TV> iterator(projection.p_grid);iterator.Valid();iterator.Next()){
            int region=levelset_multiple.Inside_Region(iterator.Cell_Index());
            if(dirichlet_regions(levelset_multiple.Inside_Region(iterator.Cell_Index()))){
                projection.elliptic_solver->psi_D(iterator.Cell_Index())=true;projection.p(iterator.Cell_Index())=(*pressures)(region);}}}}
    if(mpi_grid){
        mpi_grid->Exchange_Boundary_Cell_Data(projection.elliptic_solver->psi_D,1,nonzero_viscosity); // need to exchange corners only in the case of implicit viscosity
        mpi_grid->Exchange_Boundary_Cell_Data(projection.p,1,false);}
}
//#####################################################################
// Function Add_Surface_Tension
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
Add_Surface_Tension(LEVELSET<TV>& levelset,const T time)
{
    LAPLACE_UNIFORM<T_GRID>& elliptic_solver=*projection.elliptic_solver;T_GRID& p_grid=elliptic_solver.grid;
    LAPLACE_COLLIDABLE<T_GRID>& collidable_solver=*projection.collidable_solver;
    LEVELSET_MULTIPLE<T_GRID>& levelset_multiple=*projection.poisson_collidable->levelset_multiple;
    T_ARRAYS_SCALAR& phi=levelset.phi; 

   if(collidable_solver.second_order_cut_cell_method) for(FACE_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()){
        TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index(),face_index=iterator.Face_Index();int axis=iterator.Axis();
        if(!projection.elliptic_solver->psi_N.Component(axis)(face_index) && LEVELSET_UTILITIES<T>::Interface(phi(first_cell_index),phi(second_cell_index))){
            int region_1,region_2;T phi_1,phi_2;levelset_multiple.Minimum_Regions(first_cell_index,second_cell_index,region_1,region_2,phi_1,phi_2);
            if(!surface_tensions(region_1,region_2))continue;
            T theta=LEVELSET_UTILITIES<T>::Theta(phi(first_cell_index),phi(second_cell_index));
            TV location=theta*(grid.Center(second_cell_index)-grid.Center(first_cell_index))+grid.Center(first_cell_index);
            T curvature_at_interface=levelset.Compute_Curvature(location);
            collidable_solver.u_interface.Component(axis)(face_index)=-surface_tensions(region_1,region_2)*curvature_at_interface;}}
   else{
       for(CELL_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()) if(elliptic_solver.psi_D(iterator.Cell_Index()) && phi(iterator.Cell_Index()) < 5*grid.dX.Max()){
           int minimum_region,second_minimum_region;T minimum_phi,second_minimum_phi;
           levelset_multiple.Two_Minimum_Regions(iterator.Cell_Index(),minimum_region,second_minimum_region,minimum_phi,second_minimum_phi);
           projection.p(iterator.Cell_Index())=-surface_tensions(minimum_region,second_minimum_region)*levelset.Compute_Curvature(phi,iterator.Cell_Index());}}
}
//#####################################################################
// Function Implicit_Viscous_Update
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
Implicit_Viscous_Update(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    for(int axis=0;axis<T_GRID::dimension;axis++){
        IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<T_GRID> implicit_viscosity(projection,variable_viscosity,projection.densities,viscosities,mpi_grid,axis,use_variable_viscosity);
        implicit_viscosity.Viscous_Update(grid,face_velocities,face_velocities,dt,time,maximum_implicit_viscosity_iterations);}
}
//#####################################################################
// Function Compute_Vorticity_Confinement_Force
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
Compute_Vorticity_Confinement_Force(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,ARRAY<TV,TV_INT>& F)
{
    T_ARRAYS_SPIN vorticity(grid.Cell_Indices(2),false);T_ARRAYS_SCALAR vorticity_magnitude(grid.Cell_Indices(2));
    if(projection.flame){FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<T_GRID> face_velocities_lookup(face_velocities_ghost,projection,projection.poisson_collidable->levelset_multiple);
        VORTICITY_UNIFORM<TV>::Vorticity(grid,face_velocities_lookup,vorticity,vorticity_magnitude);}
    else VORTICITY_UNIFORM<TV>::Vorticity(grid,FACE_LOOKUP_UNIFORM<T_GRID>(face_velocities_ghost),vorticity,vorticity_magnitude);
    for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()){
        TV vortex_normal_vector=LEVELSET<TV>::Normal_At_Node(grid,vorticity_magnitude,iterator.Cell_Index());
        F(iterator.Cell_Index())=TV::Cross_Product(vortex_normal_vector,vorticity(iterator.Cell_Index()));}
}
//#####################################################################
namespace PhysBAM{
template class INCOMPRESSIBLE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,1> > >;
template class INCOMPRESSIBLE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > >;
template class INCOMPRESSIBLE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >;
template class INCOMPRESSIBLE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,1> > >;
template class INCOMPRESSIBLE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > >;
template class INCOMPRESSIBLE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >;
}
