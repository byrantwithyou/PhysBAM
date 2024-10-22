//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Grid_PDE/Boundaries/BOUNDARY_MPI.h>
#include <Grid_PDE/Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <Dynamics/Incompressible_Flows/DETONATION_SHOCK_DYNAMICS.h>
#include <Dynamics/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
#include <Dynamics/Interpolation/FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PROJECTION_DYNAMICS_UNIFORM<TV>::
PROJECTION_DYNAMICS_UNIFORM(const GRID<TV>& mac_grid,const bool flame_input,const bool multiphase,const bool use_variable_beta,const bool use_poisson)
    :PROJECTION_COLLIDABLE_UNIFORM<TV>(mac_grid,multiphase,flame_input || multiphase || use_variable_beta || use_poisson,use_variable_beta),PROJECTION_DYNAMICS<T>(flame_input),
    use_flame_speed_multiplier(0),dsd(0),use_divergence_multiplier_save_for_sph(false),use_non_zero_divergence_save_for_sph(false),
    p_save_for_sph(0),divergence_save_for_sph(0),divergence_multiplier_save_for_sph(0),face_velocities_save_for_sph(0),elliptic_solver_save_for_sph(0),laplace_save_for_sph(0),poisson_save_for_sph(0),
    collidable_solver_save_for_sph(0)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PROJECTION_DYNAMICS_UNIFORM<TV>::
PROJECTION_DYNAMICS_UNIFORM(const GRID<TV>& mac_grid,LEVELSET<TV>& levelset_input)
    :PROJECTION_COLLIDABLE_UNIFORM<TV>(mac_grid,levelset_input),PROJECTION_DYNAMICS<T>(false),use_flame_speed_multiplier(0),dsd(0),collidable_solver_save_for_sph(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PROJECTION_DYNAMICS_UNIFORM<TV>::
~PROJECTION_DYNAMICS_UNIFORM()
{
    delete dsd;
}
//#####################################################################
// Function Initialize_Grid
//#####################################################################
template<class TV> void PROJECTION_DYNAMICS_UNIFORM<TV>::
Initialize_Grid(const GRID<TV>& mac_grid)
{
    BASE::Initialize_Grid(mac_grid);
    if(dsd)dsd->Initialize_Grid();
}
//#####################################################################
// Function Initialize_Dsd
//#####################################################################
template<class TV> void PROJECTION_DYNAMICS_UNIFORM<TV>::
Initialize_Dsd(const LEVELSET_MULTIPLE<TV>& levelset_multiple,const ARRAY<bool>& fuel_region)
{
    int region=0;if(!fuel_region.Find(true,region)) PHYSBAM_FATAL_ERROR();//TODO: multiple fuel regions
    delete dsd;dsd=new DETONATION_SHOCK_DYNAMICS<TV>(p_grid,*levelset_multiple.levelsets(region));
    if(elliptic_solver->mpi_grid)
        dsd->Set_Custom_Boundary(new BOUNDARY_MPI<TV>(elliptic_solver->mpi_grid,*(dsd->boundary)),
            new BOUNDARY_MPI<TV,TV>(elliptic_solver->mpi_grid,*(dsd->boundary_vector)));
}
//#####################################################################
// Function Initialize_Dsd
//#####################################################################
template<class TV> void PROJECTION_DYNAMICS_UNIFORM<TV>::
Initialize_Dsd(const LEVELSET<TV>& levelset,const ARRAY<bool>& fuel_region)
{
    int region=0;if(!fuel_region.Find(true,region)) PHYSBAM_FATAL_ERROR();//TODO: multiple fuel regions
    delete dsd;dsd=new DETONATION_SHOCK_DYNAMICS<TV>(p_grid,levelset);
    if(elliptic_solver->mpi_grid)
        dsd->Set_Custom_Boundary(new BOUNDARY_MPI<TV>(elliptic_solver->mpi_grid,*dsd->boundary),
            new BOUNDARY_MPI<TV,TV>(elliptic_solver->mpi_grid,*dsd->boundary_vector));
}
//#####################################################################
// Function Make_Divergence_Free
//#####################################################################
template<class TV> void PROJECTION_DYNAMICS_UNIFORM<TV>::
Make_Divergence_Free(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time)
{
    // find f - divergence of the velocity
    if(flame) Compute_Divergence(T_FACE_LOOKUP_FIRE_MULTIPHASE(face_velocities,*this,poisson_collidable->levelset_multiple),elliptic_solver);
    else Compute_Divergence(T_FACE_LOOKUP(face_velocities),elliptic_solver);

    // find the pressure
    elliptic_solver->Find_Solution_Regions(); // flood fill
    elliptic_solver->Compute_beta_And_Add_Jumps_To_b(dt,time); // only does something for poisson solver
    if(elliptic_solver->solve_neumann_regions) Enforce_Velocity_Compatibility(face_velocities); // make all the right hand sides compatible
    elliptic_solver->Solve(time,true); // solve all regions

    Apply_Pressure(face_velocities,dt,time);
}
//#####################################################################
// Function Compute_Divergence
//#####################################################################
template<class TV> void PROJECTION_DYNAMICS_UNIFORM<TV>::
Compute_Divergence(const T_FACE_LOOKUP_FIRE_MULTIPHASE& face_lookup,LAPLACE_UNIFORM<TV>* solver)
{
    TV one_over_dx=p_grid.one_over_dX;
    for(CELL_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()){
        const typename T_FACE_LOOKUP_FIRE_MULTIPHASE::LOOKUP& lookup=face_lookup.Starting_Point_Cell(iterator.Cell_Index());T divergence=0;
        for(int axis=0;axis<TV::m;axis++)divergence+=(lookup(axis,iterator.Second_Face_Index(axis))-lookup(axis,iterator.First_Face_Index(axis)))*one_over_dx[axis];
        solver->f(iterator.Cell_Index())=divergence;}
    
    if(use_non_zero_divergence) for(CELL_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next())
        solver->f(iterator.Cell_Index())-=divergence(iterator.Cell_Index());
    if(use_divergence_multiplier) for(CELL_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next())
        solver->f(iterator.Cell_Index())*=divergence_multiplier(iterator.Cell_Index());
}
//#####################################################################
// Function Apply_Pressure
//#####################################################################
template<class TV> void PROJECTION_DYNAMICS_UNIFORM<TV>::
Apply_Pressure(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time,bool scale_by_dt)
{
    BASE::Apply_Pressure(face_velocities,dt,time,scale_by_dt);

    // fix the jump in pressure - interior only
    if(poisson && poisson->u_jumps){
        ARRAY<bool,TV_INT>& psi_D=elliptic_solver->psi_D;
        ARRAY<bool,FACE_INDEX<TV::m> >& psi_N=elliptic_solver->psi_N;
        TV dx=p_grid.dX,one_over_dx=Inverse(dx);
        int ghost_cells=1;
        if(poisson->multiphase){
            ARRAY<ARRAY<T,TV_INT>> phis_ghost;phis_ghost.Resize(poisson_collidable->levelset_multiple->levelsets.m);
            for(int i=0;i<poisson_collidable->levelset_multiple->levelsets.m;i++){
                phis_ghost(i).Resize(p_grid.Domain_Indices(ghost_cells),no_init);
                poisson_collidable->levelset_multiple->levelsets(i)->boundary->Fill_Ghost_Cells(p_grid,poisson_collidable->levelset_multiple->levelsets(i)->phi,phis_ghost(i),dt,time,ghost_cells);}
            LEVELSET_MULTIPLE<TV> levelset_multiple(p_grid,phis_ghost);
            for(FACE_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()){
                int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index(),first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
                if(levelset_multiple.Interface(second_cell,first_cell) && !psi_N.Component(axis)(face_index) && !(psi_D(second_cell)&&psi_D(first_cell))){
                    int region_1,region_2;T phi_1,phi_2;levelset_multiple.Minimum_Regions(second_cell,first_cell,region_1,region_2,phi_1,phi_2);
                    face_velocities.Component(axis)(face_index)+=poisson->beta_face.Component(axis)(face_index)*one_over_dx[axis]*
                        LEVELSET_MULTIPLE<TV>::Sign(region_1,region_2)*poisson_collidable->u_jump_face.Component(axis)(face_index);}}}
        else{
            ARRAY<T,TV_INT> phi_ghost(p_grid.Domain_Indices(ghost_cells));poisson_collidable->levelset->boundary->Fill_Ghost_Cells(p_grid,poisson_collidable->levelset->phi,phi_ghost,dt,time,ghost_cells);
            for(FACE_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()){
                int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index(),first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
                if(LEVELSET_UTILITIES<T>::Interface(phi_ghost(second_cell),phi_ghost(first_cell)) && !psi_N.Component(axis)(face_index) && !(psi_D(second_cell)&&psi_D(first_cell))){
                    face_velocities.Component(axis)(face_index)+=poisson->beta_face.Component(axis)(face_index)*one_over_dx[axis]*LEVELSET_UTILITIES<T>::Sign(phi_ghost(second_cell))*
                        LEVELSET_UTILITIES<T>::Average(phi_ghost(second_cell),poisson_collidable->u_jump(second_cell),phi_ghost(first_cell),poisson_collidable->u_jump(first_cell));}}}}
}
//#####################################################################
// Function Set_Up_For_SPH
//#####################################################################
template<class TV> void PROJECTION_DYNAMICS_UNIFORM<TV>::
Set_Up_For_SPH(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const bool use_variable_density_solve,const bool use_one_way_coupling)
{
    if(use_variable_density_solve){
        POISSON_COLLIDABLE_UNIFORM<TV>* poisson_for_sph=new POISSON_COLLIDABLE_UNIFORM<TV>(p_grid,p,true,false,true);
        laplace_save_for_sph=laplace_collidable;elliptic_solver_save_for_sph=elliptic_solver;poisson_save_for_sph=poisson_collidable;
        collidable_solver_save_for_sph=collidable_solver;
        elliptic_solver=poisson=poisson_collidable=poisson_for_sph;laplace=laplace_collidable=0;
        collidable_solver=poisson_collidable;
        poisson->Solve_Neumann_Regions(elliptic_solver_save_for_sph->solve_neumann_regions);
        poisson->Set_Relative_Tolerance(elliptic_solver_save_for_sph->relative_tolerance);
        poisson->pcg.Set_Maximum_Iterations(elliptic_solver_save_for_sph->pcg.maximum_iterations);
        poisson->pcg.Show_Results();
        poisson->psi_N=elliptic_solver_save_for_sph->psi_N;poisson_for_sph->psi_D=elliptic_solver_save_for_sph->psi_D;
        poisson->mpi_grid=elliptic_solver_save_for_sph->mpi_grid;
        poisson_collidable->levelset=collidable_solver_save_for_sph->levelset;
        if(use_one_way_coupling){
            use_non_zero_divergence_save_for_sph=use_non_zero_divergence;
            use_divergence_multiplier_save_for_sph=use_divergence_multiplier;
            for(FACE_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()){
                TV_INT cell_1=iterator.First_Cell_Index(),cell_2=iterator.Second_Cell_Index();
                if(!elliptic_solver_save_for_sph->psi_D(cell_1) || !elliptic_solver_save_for_sph->psi_D(cell_2)) elliptic_solver->psi_N(iterator.Axis(),iterator.Face_Index())=true;}}}
    else if(use_one_way_coupling){
        face_velocities_save_for_sph=new ARRAY<T,FACE_INDEX<TV::m> >(face_velocities);
        p_save_for_sph=new ARRAY<T,TV_INT>(p);
        divergence_save_for_sph=new ARRAY<T,TV_INT>(divergence);
        divergence_multiplier_save_for_sph=new ARRAY<T,TV_INT>(divergence_multiplier);
        use_divergence_multiplier_save_for_sph=use_divergence_multiplier;
        use_non_zero_divergence_save_for_sph=use_non_zero_divergence;
        elliptic_solver->psi_D_save_for_sph=new ARRAY<bool,TV_INT>(elliptic_solver->psi_D);
        elliptic_solver->psi_N_save_for_sph=new ARRAY<bool,FACE_INDEX<TV::m> >(elliptic_solver->psi_N);
        for(FACE_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()){
            TV_INT cell_1=iterator.First_Cell_Index(),cell_2=iterator.Second_Cell_Index();
            if(!(*elliptic_solver->psi_D_save_for_sph)(cell_1) || !(*elliptic_solver->psi_D_save_for_sph)(cell_2)) elliptic_solver->psi_N(iterator.Axis(),iterator.Face_Index())=true;}}

    Use_Divergence_Multiplier(true);Use_Non_Zero_Divergence(true);
}
//#####################################################################
// Function Restore_After_SPH
//#####################################################################
template<class TV> void PROJECTION_DYNAMICS_UNIFORM<TV>::
Restore_After_SPH(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const bool use_variable_density_solve,const bool use_one_way_coupling)
{
    if(use_variable_density_solve){
        delete poisson;poisson=poisson_collidable=poisson_save_for_sph;
        laplace=laplace_collidable=laplace_save_for_sph;elliptic_solver=elliptic_solver_save_for_sph;collidable_solver=collidable_solver_save_for_sph;
        if(use_one_way_coupling){Use_Divergence_Multiplier(use_divergence_multiplier_save_for_sph);Use_Non_Zero_Divergence(use_non_zero_divergence_save_for_sph);}}
    else if(use_one_way_coupling){
        face_velocities=*face_velocities_save_for_sph;
        delete face_velocities_save_for_sph;face_velocities_save_for_sph=0;
        p=*p_save_for_sph;
        delete p_save_for_sph;p_save_for_sph=0;
        divergence=*divergence_save_for_sph;
        delete divergence_save_for_sph;divergence_save_for_sph=0;
        divergence_multiplier=*divergence_multiplier_save_for_sph;
        delete divergence_multiplier_save_for_sph;divergence_multiplier_save_for_sph=0;
        Use_Divergence_Multiplier(use_divergence_multiplier_save_for_sph);Use_Non_Zero_Divergence(use_non_zero_divergence_save_for_sph);
        elliptic_solver->psi_D=*elliptic_solver->psi_D_save_for_sph;
        delete elliptic_solver->psi_D_save_for_sph;elliptic_solver->psi_D_save_for_sph=0;
        elliptic_solver->psi_N=*elliptic_solver->psi_N_save_for_sph;
        delete elliptic_solver->psi_N_save_for_sph;elliptic_solver->psi_N_save_for_sph=0;}
}
//#####################################################################
// Function Update_Phi_And_Move_Velocity_Discontinuity
//#####################################################################
template<class TV> void PROJECTION_DYNAMICS_UNIFORM<TV>::
Update_Phi_And_Move_Velocity_Discontinuity(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,LEVELSET_MULTIPLE<TV>& levelset_multiple,const T time,const bool update_phi_only)
{
    assert(flame);
    int ghost_cells=3;
    levelset_multiple.Fill_Ghost_Cells(levelset_multiple.phis,time,ghost_cells);
    if(!update_phi_only){
        LEVELSET_MULTIPLE<TV>& levelset_multiple_old=*poisson_collidable->levelset_multiple;
        for(FACE_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face=iterator.Face_Index();
            int region_old=levelset_multiple_old.Inside_Region_Face(axis,face),region_new=levelset_multiple.Inside_Region_Face(axis,face);
            if(region_old!=region_new) face_velocities.Component(axis)(face)-=Face_Jump_Multiphase(axis,face,region_new,region_old);}}
    poisson_collidable->Update_Internal_Level_Set(levelset_multiple);poisson_collidable->levelset_multiple->Compute_Normals();
}
//#####################################################################
// Function Update_Phi_And_Move_Velocity_Discontinuity
//#####################################################################
template<class TV> typename TV::SCALAR PROJECTION_DYNAMICS_UNIFORM<TV>::
Flame_Speed_Face_Multiphase(const int axis,const TV_INT& face_index,const int fuel_region,const int product_region) const
{
    T multiplier=1;if(use_flame_speed_multiplier) multiplier=flame_speed_multiplier.Component(axis)(face_index);
    if(dsd) return multiplier*dsd->Normal_Flame_Speed(axis,face_index);
    TV_INT offset=TV_INT::Axis_Vector(axis);const TRIPLE<T,T,T>& constants=flame_speed_constants(fuel_region,product_region);
    const T normal_flame_speed=constants.x;const T curvature_flame_speed=constants.y;
    if(!curvature_flame_speed) return multiplier*normal_flame_speed;
    const LEVELSET<TV>* levelset=poisson_collidable->levelset_multiple->levelsets(fuel_region);
    T face_curvature=(T).5*((*levelset->curvature)(face_index)+(*levelset->curvature)(face_index-offset));
    return multiplier*(normal_flame_speed+curvature_flame_speed*face_curvature);
}
//#####################################################################
// Function Use_Flame_Speed_Multiplier
//#####################################################################
template<class TV> void PROJECTION_DYNAMICS_UNIFORM<TV>::
Use_Flame_Speed_Multiplier(const bool use_flame_speed_multiplier_input)
{
    use_flame_speed_multiplier=use_flame_speed_multiplier_input;
    if(use_flame_speed_multiplier) flame_speed_multiplier.Resize(p_grid,3);
    else flame_speed_multiplier.Clean_Memory();
}
//#####################################################################
// Function Face_Jump_Multiphase
//#####################################################################
template<class TV> typename TV::SCALAR PROJECTION_DYNAMICS_UNIFORM<TV>::
Face_Jump_Multiphase(const int axis,const TV_INT& face_index,const int current_region,const int face_region) const
{
    const TRIPLE<T,T,T>& constants=flame_speed_constants(current_region,face_region);
    if(constants.z==0) return 0; // getting the diagonal if the two regions are the same. i.e. .z will be zero and the jump will return as nothing
    TV_INT offset=TV_INT::Axis_Vector(axis);
    const LEVELSET<TV>* levelset=poisson_collidable->levelset_multiple->levelsets(current_region);
    T face_normal=(levelset->phi(face_index)-levelset->phi(face_index-offset))*p_grid.one_over_dX[axis];
    return constants.z*Flame_Speed_Face_Multiphase(axis,face_index,current_region,face_region)*face_normal;
}
//#####################################################################
namespace PhysBAM{
template class PROJECTION_DYNAMICS_UNIFORM<VECTOR<float,1> >;
template class PROJECTION_DYNAMICS_UNIFORM<VECTOR<float,2> >;
template class PROJECTION_DYNAMICS_UNIFORM<VECTOR<float,3> >;
template class PROJECTION_DYNAMICS_UNIFORM<VECTOR<double,1> >;
template class PROJECTION_DYNAMICS_UNIFORM<VECTOR<double,2> >;
template class PROJECTION_DYNAMICS_UNIFORM<VECTOR<double,3> >;
}
